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
tests/oracle/test_special -all --seed=12345
tests/oracle/test_identities -all --seed=12345
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
`test_special` checks NaN/Inf constants, documented domain-NaN behavior,
selected endpoint values, `nint` half cases, exact power-of-two scaling,
and `floor`/`ceil`/`aint` where those APIs exist.
`test_identities` checks composed identities such as `exp(log(x))`,
`log(exp(x))`, `sin^2(x)+cos^2(x)`, `sqrt(x)^2`, `nroot(x,3)^3`,
`tan(x)` versus `sin(x)/cos(x)`, and inverse trig principal-branch
round trips over bounded domains.

Seed precedence is:

1. `QD_TEST_SEED`
2. `--seed=N`
3. the fixed default `0x9E3779B97F4A7C15`

Every program prints TAP version 13 output. On failure, the diagnostic
block includes the active seed, replay command, input limbs in C99 hex,
the MPFR reference, result limbs, measured error in eps, allowed bound,
and the bound justification.
