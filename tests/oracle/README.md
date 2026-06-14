# MPFR Oracle Tests

This directory contains the optional MPFR-backed oracle test suite. The
suite is not part of the default build; it is enabled only when MPFR is
requested explicitly.

## Build

Autotools:

```sh
./configure --enable-mpfr-tests
make check
```

CMake:

```sh
cmake -S . -B build-mpfr -DQD3_ENABLE_MPFR_TESTS=ON
cmake --build build-mpfr -j
ctest --test-dir build-mpfr --output-on-failure
```

## Run

The oracle programs accept the same
precision selection flags as `qd_test`:

```sh
tests/oracle/test_arith -dd
tests/oracle/test_arith -td
tests/oracle/test_arith -qd
tests/oracle/test_arith -edd
tests/oracle/test_arith -all --seed=12345

tests/oracle/test_unary -dd
tests/oracle/test_unary -td
tests/oracle/test_unary -qd
tests/oracle/test_unary -edd
tests/oracle/test_unary -all --seed=12345

tests/oracle/test_binary -all --seed=12345
tests/oracle/test_special -all
tests/oracle/test_identities -all --seed=12345
tests/oracle/test_io -all --seed=12345
tests/oracle/test_capi -all --seed=12345

tests/oracle/test_rounding_corners -all --seed=12345
```

`test_arith` checks `+`, `-`, `*`, `/`, `sqr`, compound assignment,
mixed `T op double`, and comparisons against exact MPFR images.
`test_unary` checks `sqrt`, `sqr`, `exp`, `log`, `log10`, `sin`, `cos`,
`tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, and `tanh` over bounded
random domains. Trigonometric tests are split into stable inputs, where
both `|sin(x)|` and `|cos(x)|` are at least `0.25`, and conditioned
inputs, where the allowed error is scaled by the relative condition
number. The TAP comments for trig cases include the worst input limbs,
MPFR `sin`, `cos`, `tan`, the condition number, and the `tan` error
obtained from the library's `sin(x) / cos(x)` path.
`test_binary` checks two-input functions and parameterized operations:
`pow(T,int)`, `pow(T,T)` where supported, `nroot`, `atan2`, `ldexp`,
`fmod`, and `drem`.
`test_special` is deterministic and does not accept `--seed`; it checks
NaN/Inf constants, documented domain-NaN behavior,
selected endpoint values, `nint` half cases, exact power-of-two scaling,
and `floor`/`ceil`/`aint` where those APIs exist.
`test_identities` checks composed identities such as `exp(log(x))`,
`log(exp(x))`, `sin^2(x)+cos^2(x)`, `sqrt(x)^2`, `nroot(x,3)^3`,
`tan(x)` versus `sin(x)/cos(x)`, and inverse trig principal-branch
round trips over bounded domains.
`test_io` checks finite full-precision bounded string, stream, and
`write` round trips, reduced-precision bounded round trips, format-grid invariants,
and the type-specific special-string parse contract.
`test_capi` checks the C wrappers against the C++ API with bit-identical
limb comparisons for arithmetic, unary, transcendental, trigonometric,
hyperbolic, comparison, constants, read, and swrite shims.

## Oracle Coverage Audit

This audit compares the public C++ and C API surfaces in `include/qd/*.h`
and `include/qd/c_*.h` with the current MPFR oracle programs.

Legend: Covered means there is a direct MPFR oracle, exact expansion oracle,
or bit-identical C shim oracle for the listed surface. Partial means the
surface is sampled or represented, but not every overload, boundary, or
parameter family is covered directly. TODO means the public surface still
needs an explicit oracle row or special-value grid.

### Registry and API Matrix

| Surface | Current oracle | Status | Notes |
| --- | --- | --- | --- |
| `+`, `-`, `*`, `/`, `sqr` | `test_arith`, `test_rounding_corners` | Covered | Same-type random MPFR checks plus deterministic exactness, tie, and variant cases. |
| Compound assignment | `test_arith` | Partial | Same-type compound operations are checked against the corresponding expression; mixed compound overloads still need a dedicated matrix. |
| Comparisons | `test_arith`, `test_capi` | Partial | Same-type comparisons are randomized against MPFR ordering; mixed dd/td/qd comparison overloads are not exhaustively enumerated. |
| `sqrt`, `sqr`, `exp`, `log`, `log10` | `fn_registry`, `test_unary` | Covered | Random domains are table-driven and report bound justification. |
| `sin`, `cos`, `tan` | `fn_registry`, `test_unary` | Covered | Stable and conditioned domains are split; trig diagnostics include condition number, MPFR components, and td reduction sectors. |
| `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh` | `fn_registry`, `test_unary` | Covered | Random bounded domains are covered for all enabled real types. |
| `pow(T,int)` | `test_binary` | Covered | Covered for dd, td, qd, and edd. |
| `pow(T,T)` | `test_binary` | Covered where available | Covered for dd, td, and qd; edd does not expose this overload. |
| `nroot` | `test_binary`, `test_special`, `test_identities` | Partial | `n=3` is covered, including a negative cube-root edge; other integer roots remain TODO. |
| `atan2` | `test_binary`, `test_special` | Covered | Random MPFR checks plus a quadrant-II endpoint case. |
| `ldexp` | `test_binary`, `test_special` | Covered | Exact MPFR scaling with zero allowed error. |
| `fmod`, `drem` | `test_binary` | Covered where available | Covered for dd and qd, matching the public C++ API. |
| `abs`, `fabs` | `test_capi`, `test_special` | Partial | C shim equivalence and `abs(-inf)` are checked; direct MPFR absolute-value rows are still TODO. |
| `nint`, `floor`, `ceil`, `aint`, `quick_nint` | `test_special`, `test_capi` | Partial | Half cases and C shim checks exist for exposed APIs; a direct MPFR rounding grid and `quick_nint` coverage remain TODO. |
| `sincos`, `sincosh` | `test_identities` indirectly | TODO | `tan` versus `sin/cos` is checked, but the public pair-output APIs need direct bit or MPFR checks. |
| `asinh`, `acosh`, `atanh` | none | TODO | Public C++ and C headers expose these for dd/td/qd; add direct MPFR rows. |
| String construction, stream read/write, `to_string`, `to_digits`, `write` | `test_io`, `test_capi` | Partial | Round trips and format invariants are covered; exact decimal boundary and locale/width grids remain TODO. |
| C API wrappers | `test_capi` | Partial | Core shims are compared to C++ API bit-identically; mixed C shim overloads, `sincos`, `sincosh`, inverse hyperbolic, and write/rand variants need more rows. |

### Special-Value Matrix

| Special class | Current coverage | Status | Remaining work |
| --- | --- | --- | --- |
| `NaN`, `+inf`, `-inf` constants | `test_special`, `test_io` | Partial | Constants, formatting, and selected domain-NaN contracts are covered; arithmetic propagation rows remain TODO. |
| Signed zero | `test_rounding_corners`, `test_special` | Partial | Exact zero identities and `atan2(0,0)` domain behavior are covered; signed-zero propagation by operation remains TODO. |
| Subnormal, min-normal, max, safe-max | none systematic | TODO | Add generator rows for `_min_normalized`, `_max`, `_safe_max`, nearby values, and subnormal-like limb patterns where representable. |
| Powers of two | `test_binary`, `test_rounding_corners`, `test_special` | Covered | `ldexp` and exact scaling cases exercise power-of-two behavior. |
| Domain endpoints | `test_special` | Partial | `asin(+1)`, `acos(-1)`, `atan2` quadrant II, `nroot(-8,3)`, `sqrt(-1)`, and `log(-1)` are covered; `log(0)`, `pow(0,0)`, and out-of-domain asin/acos need rows. |
| Trig zeros and poles | `test_unary` | Partial | Conditioned random tests cover near-zero and near-pole cases statistically; deterministic grids around multiples of pi/2 remain TODO. |
| Overflow and underflow | legacy `huge` plus constants | Partial | Add MPFR oracle rows for overflow, underflow, and safe max arithmetic once policy expectations are fixed. |


## Sanitizers, Coverage, and Matrix QA

CMake sanitizer run:

```sh
cmake -S . -B build-sanitize \
  -DQD3_ENABLE_MPFR_TESTS=ON \
  -DQD3_ENABLE_ASAN=ON \
  -DQD3_ENABLE_UBSAN=ON \
  -DBUILD_TESTING=ON
cmake --build build-sanitize -j
ctest --test-dir build-sanitize --output-on-failure
```

Autotools sanitizer run:

```sh
./configure --enable-mpfr-tests --enable-sanitizers=address,undefined
make check
```

CMake coverage run:

```sh
cmake -S . -B build-coverage \
  -DQD3_ENABLE_MPFR_TESTS=ON \
  -DQD3_ENABLE_COVERAGE=ON \
  -DBUILD_TESTING=ON
cmake --build build-coverage --target coverage -j
```

The coverage target writes `coverage.filtered.info` and
`coverage-html/index.html` under the build directory. The latest local
coverage run with GCC 15.2.0, MPFR 4.2.2, and the MPFR oracle suite
enabled passed 41/41 tests. Filtered lcov
coverage over `src/*.cpp`, `src/td_trig_reduce.h`, and
`include/qd/*_inline.h` was 69.7% line coverage and 75.2% function
coverage. The follow-up target
is to raise exercised arithmetic/transcendental function-body coverage to
at least 90% as the remaining corner-case tests are added.


## Unary Worst-Case Snapshot (seed 12345)

The following table is a representative comparison of the maximum observed
`relerr` per function+domain across `dd`, `td`, `qd`, and `edd`.

Command used:

```sh
build-mpfr/tests/oracle/test_unary -dd --seed=12345 --worst-report=/tmp/oracle_unary_dd_worst.csv
build-mpfr/tests/oracle/test_unary -td --seed=12345 --worst-report=/tmp/oracle_unary_td_worst.csv
build-mpfr/tests/oracle/test_unary -qd --seed=12345 --worst-report=/tmp/oracle_unary_qd_worst.csv
build-mpfr/tests/oracle/test_unary -edd --seed=12345 --worst-report=/tmp/oracle_unary_edd_worst.csv
python3 qa/compare_unary_worst_reports.py \
  --report dd=/tmp/oracle_unary_dd_worst.csv \
  --report td=/tmp/oracle_unary_td_worst.csv \
  --report qd=/tmp/oracle_unary_qd_worst.csv \
  --report edd=/tmp/oracle_unary_edd_worst.csv
```

| function | operation | domain | dd | td | qd | edd | best_precision | best_relerr | worst_relerr |
|---|---|---|---|---|---|---|---|---|
| cos_conditioned | cos | trig_conditioned | 234.5158837015159 | 0.1859678578838277 | 3.606068198292264 | 0.5842733118180256 | td:0.1859678578838277 | 234.5158837015159 |
| tan_conditioned | tan | trig_conditioned | 100.2354127653254 | 0.4293062542108887 | 5.337407515500874 | 0.7403564850562804 | td:0.4293062542108887 | 100.2354127653254 |
| tan_stable | tan | trig_stable | 6.894465606911766 | 0.6904400371696184 | 0.1798864565349701 | 0.5534120362369035 | qd:0.1798864565349701 | 6.894465606911766 |
| cos_stable | cos | trig_stable | 6.627038246121096 | 0.1448689713104189 | 0.1563370472712637 | 0.3846866974570015 | td:0.1448689713104189 | 6.627038246121096 |
| sin_stable | sin | trig_stable | 5.558650157962114 | 0.1646533043055813 | 0.1356551055455925 | 0.4341243022129826 | qd:0.1356551055455925 | 5.558650157962114 |
| sin_conditioned | sin | trig_conditioned | 4.061869175123492 | 0.2586139619497284 | 0.484758983433543 | 0.7634140708511088 | td:0.2586139619497284 | 4.061869175123492 |
| log10 | log10 | positive | 3.382881796619329 | 1.451569741775147 | 0.5365655341154324 | 0.2169708880475239 | edd:0.2169708880475239 | 3.382881796619329 |
| log | log | positive | 3.188037618468893 | 0.1733815365736716 | 0.02628300239941649 | 0.22259460670651 | qd:0.02628300239941649 | 3.188037618468893 |
| tanh | tanh | hyperbolic_moderate | 3.040973951650837 | 1.2649942359957 | 0.2710157193589349 | 0.9181692205563492 | qd:0.2710157193589349 | 3.040973951650837 |
| sinh | sinh | hyperbolic_moderate | 2.437613433479258 | 0.4052122187515949 | 0.1673047257728243 | 0.6210456354884193 | qd:0.1673047257728243 | 2.437613433479258 |
| exp | exp | exp_moderate | 1.789884216404387 | 0.7173980805162762 | 0.04589685531423539 | 0.7639450291224287 | qd:0.04589685531423539 | 1.789884216404387 |
| asin | asin | unit | 1.078179715832877 | 0.2747186105239964 | 0.0912821036830738 | 0.5494034213667686 | qd:0.0912821036830738 | 1.078179715832877 |
| atan | atan | all_moderate | 0.9421272497469996 | 0.2171871948503973 | 0.6420795173406116 | 0.9853352078363772 | td:0.2171871948503973 | 0.9853352078363772 |
| sqrt | sqrt | nonnegative | 0.6002603878042156 | 0.06519135199933343 | 0.2240958605086439 | 0.9717311669399569 | td:0.06519135199933343 | 0.9717311669399569 |
| cosh | cosh | hyperbolic_moderate | 0.856867811548415 | 0.4996760014112317 | 0.08926190204519797 | 0.7348719673112949 | qd:0.08926190204519797 | 0.856867811548415 |
| sqr | sqr | all_moderate | 0.5849260506966909 | 0.04975007369238089 | 0.05073075808077032 | 0.4151285873497249 | td:0.04975007369238089 | 0.5849260506966909 |
| acos | acos | unit | 0.2374165876842437 | 0.1021943776724964 | 0.03530422376174445 | 0.1447475598100737 | qd:0.03530422376174445 | 0.2374165876842437 |


Build-matrix QA:

```sh
qa/check_16_builds_cmake.sh
qa/check_oracle_matrix_cmake.sh
```

`qa/check_oracle_matrix_cmake.sh` runs the same CMake matrix with
`QD3_ENABLE_MPFR_TESTS=ON` in every configuration. For Autotools-generated
trees, `qa/check_oracle_matrix.sh` runs the MPFR oracle matrix through
`./configure --enable-mpfr-tests`.

For seeded oracle programs, seed precedence is:

1. `QD_TEST_SEED`
2. `--seed=N`
3. the fixed default `0x9E3779B97F4A7C15`

Every oracle program prints TAP version 13 output. By default, passing cases stay
compact and detailed diagnostics are emitted only for failures. Pass `-v` or
`-verbose` to emit TAP YAML diagnostics for passing MPFR comparison cases too.
For seeded oracle programs, those verbose blocks include the active seed and
replay command. For deterministic `test_special`, they include a replay command
without seed. Verbose blocks also include input limbs in
C99 hex, MPFR-formatted input values, a digit ruler above the MPFR reference,
the MPFR reference, the result converted back to MPFR, result limbs, absolute
MPFR error, measured error in eps, an `ulp_error_estimate` alias for the
eps-scaled error, the allowed bound, and the bound justification.

Reproducibility and governance notes:

- Seeded oracle programs print the active seed in TAP output and failure diagnostics.
- `test_special` is deterministic: it does not accept `--seed`, does not read
  `QD_TEST_SEED`, and prints replay commands without seed.
- For seeded programs, `QD_TEST_SEED` has highest precedence, then
  `--seed=<N>`, and the fallback default is `0x9E3779B97F4A7C15`.
- The command to replay any failure is printed in diagnostics. Seeded programs
  include `--seed=<N>`; deterministic `test_special` does not.
- Seed handling and MPFR-suite commands are aligned across Autotools and
  CMake so CI and local runs are reproducible.

`test_rounding_corners.cpp` is implemented as a dedicated corner-case oracle
and is wired into both Autotools and CMake oracle test matrices.

