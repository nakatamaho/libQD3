#!/usr/bin/env python3
"""Generate tests/c_api_test.cpp from the public C API headers.

Every function declared in include/qd/c_{dd,td,qd,edd}.h is classified from
its name (operation, operand types, result type) and exercised against a
quad-double reference with an epsilon-relative tolerance of the result type.
The generator fails if a declared function cannot be classified, so a new C
API entry point cannot be added without test coverage.

Usage: python3 qa/gen_c_api_test.py [--check] [SRC_DIR]

With --check, the generated text is compared with tests/c_api_test.cpp and
the script fails if the checked-in file is stale.
"""

import os
import re
import sys

CHECK = '--check' in sys.argv[1:]
_ARGS = [a for a in sys.argv[1:] if a != '--check']
SRC = _ARGS[0] if _ARGS else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TYPES = ['dd', 'td', 'qd', 'edd']
LIMB = {'dd': 'double', 'td': 'double', 'qd': 'double', 'edd': '_Float64x'}
NLIMB = {'dd': 2, 'td': 3, 'qd': 4, 'edd': 2}
UNARY = {
    # name: (reference expression on qd_real x, domain, tolerance in eps)
    'sqrt': ('sqrt(x)', 'positive', 4), 'sqr': ('sqr(x)', 'any', 4),
    'abs': ('abs(x)', 'any', 0), 'neg': ('-x', 'any', 0),
    'exp': ('exp(x)', 'any', 64), 'log': ('log(x)', 'positive', 64),
    'log10': ('log10(x)', 'positive', 64),
    'sin': ('sin(x)', 'any', 64), 'cos': ('cos(x)', 'any', 64),
    'tan': ('tan(x)', 'any', 64), 'asin': ('asin(x)', 'unit', 64),
    'acos': ('acos(x)', 'unit', 64), 'atan': ('atan(x)', 'any', 64),
    'sinh': ('sinh(x)', 'any', 64), 'cosh': ('cosh(x)', 'any', 64),
    'tanh': ('tanh(x)', 'any', 64), 'asinh': ('asinh(x)', 'any', 64),
    'acosh': ('acosh(x)', 'above1', 64), 'atanh': ('atanh(x)', 'unit', 64),
    'nint': ('nint(x)', 'wide', 0), 'aint': ('aint(x)', 'wide', 0),
    'floor': ('floor(x)', 'wide', 0), 'ceil': ('ceil(x)', 'wide', 0),
}
OPS = {'add': '+', 'sub': '-', 'mul': '*', 'div': '/'}


def declarations():
    decls = []
    for t in TYPES:
        text = open(os.path.join(SRC, 'include', 'qd', f'c_{t}.h')).read()
        for m in re.finditer(r'^\s*(\w+)\s+(c_' + t + r'_\w+)\s*\(([^)]*)\);', text, re.M):
            decls.append((t, m.group(2), m.group(3)))
    return decls


def arg(kind, name):
    """C argument expression for an operand of the given kind."""
    return name if kind == 'd' else f'{name}.limbs'


def operand(kind, var):
    return f'Operand<{kind}_tag> {var}' if kind != 'd' else f'double {var}'


def make(kind, var, domain='any'):
    if kind == 'd':
        return f'  double {var} = random_double(rng, Domain::{domain});\n'
    return f'  Operand<{kind}_tag> {var}(random_qd(rng, Domain::{domain}));\n'


def classify(t, name):
    rest = name[len(f'c_{t}_'):]
    m = re.fullmatch(r'(add|sub|mul|div)(?:_(d|dd|td|qd|edd)_(d|dd|td|qd|edd))?', rest)
    if m:
        a, b = (m.group(2) or t), (m.group(3) or t)
        return ('binary', m.group(1), a, b)
    m = re.fullmatch(r'self(add|sub|mul|div)(?:_(d|dd|td|qd|edd))?', rest)
    if m:
        return ('self', m.group(1), m.group(2) or t)
    m = re.fullmatch(r'copy(?:_(d|dd|td|qd|edd))?', rest)
    if m:
        return ('copy', m.group(1) or t)
    m = re.fullmatch(r'comp(?:_(d|dd|td|qd|edd)_(d|dd|td|qd|edd))?', rest)
    if m:
        return ('comp', m.group(1) or t, m.group(2) or t)
    if rest in UNARY:
        return ('unary', rest)
    if rest in ('sincos', 'sincosh'):
        return ('pair', rest)
    if rest in ('npwr', 'nroot', 'atan2', 'read', 'swrite', 'write', 'rand',
                'pi', '2pi', 'epsilon'):
        return (rest,)
    return None


def emit_case(t, name, cls):
    T = f'{t}_tag'
    out = f'void test_{name}(Rng &rng, Report &report) {{\n'
    out += '  for (int iter = 0; iter < kIterations; ++iter) {\n'
    body = ''
    kind = cls[0]
    if kind == 'binary':
        op, a, b = cls[1], cls[2], cls[3]
        body += make(a, 'a') + make(b, 'b', 'nonzero' if op == 'div' else 'any')
        body += f'  Operand<{T}> c;\n  {name}({arg(a, "a")}, {arg(b, "b")}, c.limbs);\n'
        body += f'  report.close("{name}", c.value(), ref_of(a) {OPS[op]} ref_of(b), Operand<{T}>::eps(), 4);\n'
    elif kind == 'self':
        op, a = cls[1], cls[2]
        body += make(a, 'a', 'nonzero' if op == 'div' else 'any') + make(t, 'b')
        body += f'  const qd_real expected = b.value() {OPS[op]} ref_of(a);\n'
        body += f'  {name}({arg(a, "a")}, b.limbs);\n'
        body += f'  report.close("{name}", b.value(), expected, Operand<{T}>::eps(), 4);\n'
    elif kind == 'copy':
        a = cls[1]
        body += make(a, 'a') + f'  Operand<{T}> b;\n  {name}({arg(a, "a")}, b.limbs);\n'
        body += f'  report.close("{name}", b.value(), ref_of(a), Operand<{T}>::eps(), 1);\n'
    elif kind == 'comp':
        a, b = cls[1], cls[2]
        body += make(a, 'a')
        # Second operand equal to, or within a few ulps of the narrower
        # type around, the first, so mixed-precision rounding is detected.
        if b == 'd':
            body += '  double b = near_double(ref_of(a), iter);\n'
        else:
            body += f'  Operand<{b}_tag> b = near_operand<{b}_tag>(ref_of(a), iter);\n'
        body += f'  int result = 2;\n  {name}({arg(a, "a")}, {arg(b, "b")}, &result);\n'
        body += '  const qd_real diff = ref_of(a) - ref_of(b);\n'
        body += '  const int expected = diff < 0.0 ? -1 : (diff > 0.0 ? 1 : 0);\n'
        body += f'  report.check("{name}", result == expected);\n'
    elif kind == 'unary':
        fn = cls[1]
        expr, domain, tol = UNARY[fn]
        body += make(t, 'a', domain) + f'  Operand<{T}> b;\n  {name}(a.limbs, b.limbs);\n'
        body += f'  const qd_real x = a.value();\n'
        body += f'  report.close("{name}", b.value(), {expr}, Operand<{T}>::eps(), {tol});\n'
    elif kind == 'pair':
        fn = cls[1]
        s_ref, c_ref = ('sin(x)', 'cos(x)') if fn == 'sincos' else ('sinh(x)', 'cosh(x)')
        body += make(t, 'a') + f'  Operand<{T}> s, c;\n  {name}(a.limbs, s.limbs, c.limbs);\n'
        body += '  const qd_real x = a.value();\n'
        body += f'  report.close("{name} sin", s.value(), {s_ref}, Operand<{T}>::eps(), 64);\n'
        body += f'  report.close("{name} cos", c.value(), {c_ref}, Operand<{T}>::eps(), 64);\n'
    elif kind == 'npwr':
        body += make(t, 'a', 'nonzero') + f'  Operand<{T}> b;\n  const int n = iter % 7 - 3;\n'
        body += f'  {name}(a.limbs, n, b.limbs);\n'
        body += f'  report.close("{name}", b.value(), npwr(a.value(), n), Operand<{T}>::eps(), 16);\n'
    elif kind == 'nroot':
        body += make(t, 'a', 'positive') + f'  Operand<{T}> b;\n  const int n = 2 + iter % 4;\n'
        body += f'  {name}(a.limbs, n, b.limbs);\n'
        body += f'  report.close("{name}", b.value(), nroot(a.value(), n), Operand<{T}>::eps(), 16);\n'
    elif kind == 'atan2':
        body += make(t, 'a') + make(t, 'b', 'nonzero') + f'  Operand<{T}> c;\n'
        body += f'  {name}(a.limbs, b.limbs, c.limbs);\n'
        body += f'  report.close("{name}", c.value(), atan2(a.value(), b.value()), Operand<{T}>::eps(), 64);\n'
    elif kind == 'read':
        body += f'  Operand<{T}> a;\n  {name}("1.5e-3", a.limbs);\n'
        body += f'  report.close("{name}", a.value(), qd_real(3) / qd_real(2000), Operand<{T}>::eps(), 1);\n'
    elif kind == 'swrite':
        body += make(t, 'a') + '  char text[128];\n'
        body += f'  {name}(a.limbs, Operand<{T}>::digits(), text, sizeof(text));\n'
        body += '  report.close("' + name + '", qd_real(text), a.value(), Operand<' + T + '>::eps(), 8);\n'
    elif kind == 'write':
        body += f'  if (iter == 0) {{\n    Operand<{T}> a(qd_real(0.5));\n    {name}(a.limbs);\n  }}\n'
    elif kind == 'rand':
        body += f'  Operand<{T}> a;\n  {name}(a.limbs);\n'
        body += f'  report.check("{name}", a.value() >= 0.0 && a.value() < 1.0);\n'
    elif kind in ('pi', '2pi'):
        body += f'  Operand<{T}> a;\n  {name}(a.limbs);\n'
        ref = 'qd_real::_pi' if kind == 'pi' else 'qd_real::_2pi'
        body += f'  report.close("{name}", a.value(), {ref}, Operand<{T}>::eps(), 1);\n'
    elif kind == 'epsilon':
        body += f'  report.check("{name}", qd_real(static_cast<double>({name}())) == qd_real(static_cast<double>(Operand<{T}>::eps())));\n'
    out += ''.join('  ' + line + '\n' for line in body.rstrip('\n').split('\n'))
    out += '  }\n}\n'
    return out


def main():
    decls = declarations()
    unknown = [n for t, n, _ in decls if classify(t, n) is None]
    if unknown:
        sys.exit('unclassified C API functions: ' + ', '.join(unknown))
    cases = {t: [] for t in TYPES}
    for t, n, _ in decls:
        cases[t].append((n, emit_case(t, n, classify(t, n))))
    tmpl = open(os.path.join(SRC, 'tests', 'c_api_test.cpp.in')).read()
    blocks = []
    runs = []
    for t in TYPES:
        guard = t == 'edd'
        code = ''.join(c for _, c in cases[t])
        call = ''.join(f'  test_{n}(rng, report);\n' for n, _ in cases[t])
        if guard:
            code = '#ifdef QD_HAVE_EDD_REAL\n' + code + '#endif\n'
            call = '#ifdef QD_HAVE_EDD_REAL\n' + call + '#endif\n'
        blocks.append(code)
        runs.append(call)
    text = tmpl.replace('@CASES@', '\n'.join(blocks))
    text = text.replace('@RUNS_DOUBLE@', ''.join(runs[:3])).replace('@RUNS_EDD@', runs[3])
    text = text.replace('@COUNT@', str(len(decls)))
    target = os.path.join(SRC, 'tests', 'c_api_test.cpp')
    if CHECK:
        current = open(target).read() if os.path.exists(target) else ''
        if current != text:
            sys.exit('tests/c_api_test.cpp is stale; run python3 qa/gen_c_api_test.py')
        print(f'tests/c_api_test.cpp is up to date ({len(decls)} C API functions)')
        return
    open(target, 'w').write(text)
    print(f'generated tests/c_api_test.cpp covering {len(decls)} C API functions')


if __name__ == '__main__':
    main()
