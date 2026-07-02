# libQD3 QA Notes

## QA Philosophy

libQD3 QA treats correctness as layered numerical evidence, not as a single
golden run. The default tests stay dependency-light so normal builds catch API,
ABI, smoke, and regression failures quickly. The optional oracle tests provide
stronger numerical evidence by comparing libQD3 against MPFR and MPC reference
computations at higher precision and documented rounding modes.

A result is accepted only when the relevant public precision types are covered:
`dd_real`, `td_real`, `qd_real`, `edd_real` when available, and the matching
complex types. Real oracle references are computed independently with MPFR;
complex oracle references are computed independently with MPC and MPFR; endpoint
tests use fixed mathematical contracts for NaN, infinity, branch, and rounding
behavior.

The tests make the reason for acceptance explicit. Input domains are named and
separated into stable regions, conditioned regions, and deterministic endpoint
cases. Error budgets live in the oracle source, are expressed in library eps,
and are adjusted for condition number where the mathematical problem itself is
ill-conditioned. Randomized cases are reproducible through `QD_TEST_SEED` or
`--seed`, and failures print the operands, reference value, library value,
absolute error, relative error, ULP estimate, and replay command.

This is not a formal proof of the algorithms. It is a practical correctness
argument: independent high-precision oracles check the direct numerical result,
identity tests check consistency between related functions, corner tests check
documented contracts, smoke tests check the public API used by consumers, and
matrix scripts check that those guarantees survive supported build options.

This file summarizes how the QA suite covers real and complex arithmetic and
special-function behavior. Source links are included for each test path.

## Program Index

| Program | Purpose |
| --- | --- |
| [tests/qd_test.cpp](tests/qd_test.cpp) | Default legacy real-number smoke test. |
| [tests/complex_test.cpp](tests/complex_test.cpp) | Default complex smoke test and real C++ math compatibility smoke test. |
| [tests/oracle/fn_registry.h](tests/oracle/fn_registry.h) | Shared oracle operation registry, input domains, and error bounds. |
| [tests/oracle/test_arith.cpp](tests/oracle/test_arith.cpp) | MPFR-backed real arithmetic oracle. |
| [tests/oracle/test_unary.cpp](tests/oracle/test_unary.cpp) | MPFR-backed real unary special-function oracle. |
| [tests/oracle/test_binary.cpp](tests/oracle/test_binary.cpp) | MPFR-backed real binary and parameterized special-function oracle. |
| [tests/oracle/test_special.cpp](tests/oracle/test_special.cpp) | Deterministic real special-value and endpoint oracle. |
| [tests/oracle/test_identities.cpp](tests/oracle/test_identities.cpp) | Real composed-identity oracle for elementary functions. |
| [tests/oracle/test_rounding_corners.cpp](tests/oracle/test_rounding_corners.cpp) | Deterministic real arithmetic corner and build-variant oracle. |
| [tests/oracle/mpc_complex_oracle.h](tests/oracle/mpc_complex_oracle.h) | Shared MPC/MPFR complex oracle helpers and diagnostics. |
| [tests/oracle/test_complex_arith.cpp](tests/oracle/test_complex_arith.cpp) | MPC-backed complex arithmetic oracle. |
| [tests/oracle/test_complex_unary.cpp](tests/oracle/test_complex_unary.cpp) | MPC-backed complex unary special-function oracle. |
| [tests/oracle/test_complex_binary.cpp](tests/oracle/test_complex_binary.cpp) | MPC-backed complex `pow` and `ldexp` oracle. |
| [tests/oracle/test_complex_special.cpp](tests/oracle/test_complex_special.cpp) | MPC/MPFR-backed complex helper-function oracle. |

## Real Arithmetic

Real arithmetic covers `dd_real`, `td_real`, `qd_real`, and `edd_real` when
`edd_real` is available.

| Operation | Primary oracle test | Input strategy | Reference and pass criterion | Extra coverage |
| --- | --- | --- | --- | --- |
| `a + b` | [tests/oracle/test_arith.cpp](tests/oracle/test_arith.cpp) / `oracle_test_arith_*` | 64 random same-sign or same-negated-sign expansion pairs per type, exponents in `[-220, 220]` | Convert both operands to MPFR, compute `mpfr_add(..., MPFR_RNDN)`, then require finite relative error in library eps within `add_sub_bound()` from [tests/oracle/fn_registry.h](tests/oracle/fn_registry.h) | Exact add smoke with `1.25 + 2.5`; `+=` must be bit-identical to `a + b`; `T + double` is checked against MPFR; [tests/oracle/test_rounding_corners.cpp](tests/oracle/test_rounding_corners.cpp) adds exact non-overlapping powers-of-two sum, signed-zero identities, tie-adjacent add, and add-variant cases across the CMake matrix. |
| `a - b` | [tests/oracle/test_arith.cpp](tests/oracle/test_arith.cpp) / `oracle_test_arith_*` | 64 random ordered positive pairs per type, then `b` is scaled by `1/4`; both operands may be negated together | Convert to MPFR, compute `mpfr_sub(..., MPFR_RNDN)`, then check relative error against `add_sub_bound()` | `-=` must be bit-identical to `a - b`; [tests/oracle/test_rounding_corners.cpp](tests/oracle/test_rounding_corners.cpp) adds a Sterbenz exact subtraction case and signed-zero `a - a` behavior. |
| `a * b` | [tests/oracle/test_arith.cpp](tests/oracle/test_arith.cpp) / `oracle_test_arith_*` | 64 random expansion pairs per type, exponents in `[-140, 140]` | Convert to MPFR, compute `mpfr_mul(..., MPFR_RNDN)`, then check relative error against `mul_bound()` | `*=` must be bit-identical to `a * b`; `T * double` is checked against MPFR; [tests/oracle/test_rounding_corners.cpp](tests/oracle/test_rounding_corners.cpp) adds exact power-of-two scaling, tie-adjacent multiplication, and qd sloppy/accurate multiplication variant cases. |
| `a / b` | [tests/oracle/test_arith.cpp](tests/oracle/test_arith.cpp) / `oracle_test_arith_*` | 64 random numerators with exponents in `[-220, 220]` and denominators in `[-80, 80]` per type | Convert to MPFR, compute `mpfr_div(..., MPFR_RNDN)`, then check relative error against `div_bound()` | `/=` must be bit-identical to `a / b`; default [tests/qd_test.cpp](tests/qd_test.cpp) also contains division sanity and round-trip checks for the precision-specific smoke paths. |

The legacy [tests/qd_test.cpp](tests/qd_test.cpp) executable is still part of the
default CTest suite. It is not the main oracle for the table above, but it
exercises arithmetic through identity formulas, fixed regressions, mixed-mode
paths, and precision-specific sanity checks.

## Real Special Functions

Real special-function oracle coverage also covers `dd_real`, `td_real`,
`qd_real`, and available `edd_real`.

| Function group | Tested functions | Source | Input strategy | Reference and pass criterion | Extra coverage / diagnostics |
| --- | --- | --- | --- | --- | --- |
| Unary algebraic, exponential, logarithmic, and rounding functions | `sqrt`, `sqr`, `exp`, `exp2`, `expm1`, `log`, `log10`, `log1p`, `log2`, `cbrt`, `trunc`, `round` | [tests/oracle/test_unary.cpp](tests/oracle/test_unary.cpp), driven by [tests/oracle/fn_registry.h](tests/oracle/fn_registry.h) | 48 random samples per registry row. Domains include nonnegative, positive, moderate, exp-moderate, and `log1p`-specific ranges. | Calls the matching MPFR function, for example `mpfr_sqrt`, `mpfr_exp2`, `mpfr_log1p`, `mpfr_cbrt`, `mpfr_trunc`, or `mpfr_round`, then checks relative error in library eps against the registry bound. | Failure diagnostics include seed, replay command, input limbs/value, MPFR reference, library limbs/value, absolute error, relative error, ULP estimate, and bound justification. |
| Trigonometric functions | `sin`, `cos`, `tan` | [tests/oracle/test_unary.cpp](tests/oracle/test_unary.cpp), bounds in [tests/oracle/fn_registry.h](tests/oracle/fn_registry.h) | 48 samples for both stable and conditioned domains. Stable cases require the absolute values of `sin(x)` and `cos(x)` to be at least `0.25`; conditioned cases deliberately include near-zero or near-pole behavior. | Uses MPFR `sin`, `cos`, and `tan`. Conditioned cases scale the allowed error by the computed relative condition number. | Diagnostics additionally record MPFR sin/cos/tan components, `tan` via `sin/cos`, and TD argument-reduction sectors and reduced arguments. |
| Inverse trig and hyperbolic unary functions | `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh` | [tests/oracle/test_unary.cpp](tests/oracle/test_unary.cpp), registry in [tests/oracle/fn_registry.h](tests/oracle/fn_registry.h) | 48 random samples per function. `asin`/`acos` use unit-domain inputs; hyperbolic functions use moderate inputs. | Uses MPFR `asin`, `acos`, `atan`, `sinh`, `cosh`, and `tanh`, with transcendental registry bounds. | [tests/oracle/test_identities.cpp](tests/oracle/test_identities.cpp) adds principal-branch checks for `asin(sin(x))` and `atan(tan(x))`. |
| Binary or parameterized special functions | `pow(T,int)`, `pow(T,T)`, `nroot(x,3)`, `atan2`, `ldexp`, `hypot`, `fmod`, `drem` where exposed | [tests/oracle/test_binary.cpp](tests/oracle/test_binary.cpp) | 48 random samples per case, with positive bases where real powers require them and bounded exponents/scales. | Uses MPFR `pow_si`, `pow`, `rootn_ui`, `atan2`, `mul_2si`, `hypot`, `fmod`, and `remainder`. `ldexp` has a zero-error exactness bound. | DD/QD cover `drem`; DD/TD/QD/EDD cover `fmod` and `hypot` where the public API exposes them. |
| Deterministic endpoints and domain contracts | NaN/Inf constants, `sqrt(-1)`, `log(-1)`, `atan2(0,0)`, `asin(+1)`, `acos(-1)`, `atan2(1,-1)`, `nroot(-8,3)`, `nint`, `ldexp`, selected `floor`/`ceil`/`aint` | [tests/oracle/test_special.cpp](tests/oracle/test_special.cpp) | Fixed endpoint and domain samples rather than random sampling. | Uses MPFR endpoint references where applicable and exact equality for documented constants/rounding contracts. | Suppresses expected domain error messages while checking the NaN-domain contract. |
| Composed identities | `exp(log(x))`, `log(exp(x))`, `sin^2+cos^2`, `sqrt(x)^2`, `sqr(x) == x*x`, `nroot(x,3)^3`, `tan(x) == sin(x)/cos(x)` | [tests/oracle/test_identities.cpp](tests/oracle/test_identities.cpp) | 48 random samples per identity over positive, moderate, principal, or trig-stable domains. | Compares the composed library result against the expected MPFR image or independently-composed reference under `identity_bound()`. | Complements direct oracle rows by catching inconsistencies between related functions. |
| C++ compatibility smoke | `log2`, `exp2`, `expm1`, `log1p`, `hypot`, `cbrt`, `trunc`, `round`, `fmod`, `pow`, `asinh`, `acosh`, `atanh` | [tests/complex_test.cpp](tests/complex_test.cpp) | Fixed small values run for each real value type through `run_real_math_type`. | Checks simple identities and exact/simple expected values with unqualified calls after `using std::...`. | This is default-build smoke coverage, not an MPFR oracle row. |

## Complex Arithmetic

Complex arithmetic covers `dd_complex`, `td_complex`, `qd_complex`, and
`edd_complex` when `edd_real` is available. The default smoke test is
[tests/complex_test.cpp](tests/complex_test.cpp); the optional MPC oracle tests
are enabled with `QD3_ENABLE_MPC_TESTS=ON`.

| Operation | Default smoke test | MPC oracle test | Input strategy | Reference and pass criterion |
| --- | --- | --- | --- | --- |
| `z + w` | [tests/complex_test.cpp](tests/complex_test.cpp) checks `(1+2i) + (3-4i) == 4-2i`, compound `+=`, and real/complex mixed addition | [tests/oracle/test_complex_arith.cpp](tests/oracle/test_complex_arith.cpp) / `oracle_test_complex_arith_*` uses `mpc_add(..., MPC_RNDNN)` | 24 random complex pairs per type from `random_complex<Real>(-2, 2)` | Convert real and imaginary parts through MPFR into MPC, compare the library result component-wise, and require `max(real_err, imag_err)` in eps to fit `complex_arithmetic_bound()` (`256` eps). |
| `z - w` | [tests/complex_test.cpp](tests/complex_test.cpp) checks `(1+2i) - (3-4i) == -2+6i`, compound `-=`, and real/complex mixed subtraction | [tests/oracle/test_complex_arith.cpp](tests/oracle/test_complex_arith.cpp) / `oracle_test_complex_arith_*` uses `mpc_sub(..., MPC_RNDNN)` | 24 random complex pairs per type from `random_complex<Real>(-2, 2)` | Same component-wise MPC comparison and `256` eps arithmetic bound. |
| `z * w` | [tests/complex_test.cpp](tests/complex_test.cpp) checks `(1+2i) * (3-4i) == 11+2i`, compound `*=`, and real/complex mixed multiplication | [tests/oracle/test_complex_arith.cpp](tests/oracle/test_complex_arith.cpp) / `oracle_test_complex_arith_*` uses `mpc_mul(..., MPC_RNDNN)` | 24 random complex pairs per type from `random_complex<Real>(-2, 2)` | Same component-wise MPC comparison and `256` eps arithmetic bound. |
| `z / w` | [tests/complex_test.cpp](tests/complex_test.cpp) checks `(1+2i) / (3-4i)` against `-1/5 + 2/5 i`, compound `/=`, and real/complex mixed division | [tests/oracle/test_complex_arith.cpp](tests/oracle/test_complex_arith.cpp) / `oracle_test_complex_arith_*` uses `mpc_div(..., MPC_RNDNN)` | 24 random complex pairs per type; for division, the denominator is shifted by `1.25 - 0.75i` to avoid near-zero random denominators | Same component-wise MPC comparison and `256` eps arithmetic bound. |

Complex diagnostics include the seed, replay command, input limbs, MPFR/MPC
reference values, library output limbs, absolute error, relative error in eps,
and the active MPC rounding mode.

## Complex Special Functions

Complex special-function coverage is split between default smoke tests and the
optional MPC oracle suite.

| Function group | Tested functions | Source | Input strategy | Reference and pass criterion | Extra coverage / diagnostics |
| --- | --- | --- | --- | --- | --- |
| Default complex smoke and ADL compatibility | Construction/conversion, `real`, `imag`, `conj`, `abs`, `norm`, `arg`, `polar`, `sqrt`, `exp`, `log`, `log10`, `log2`, `sin`, `cos`, `tan`, inverse trig, hyperbolic and inverse hyperbolic functions, `pow`, `ldexp`, `ceil`, `proj` | [tests/complex_test.cpp](tests/complex_test.cpp) | Fixed representative complex values, mixed real/complex operands, and `using std::...; f(z)` calls. | Checks exact simple arithmetic, finite outputs, and identities such as `sqrt(z)^2`, `exp(log(z))`, `sin^2+cos^2`, `tan == sin/cos`, `tanh == sinh/cosh`, inverse-function round trips, and `proj` signed-zero behavior. | Runs in the default CTest suite without MPC. |
| MPC unary elementary functions | `sqrt`, `exp`, `log`, `log10`, `log2`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh` | [tests/oracle/test_complex_unary.cpp](tests/oracle/test_complex_unary.cpp), helpers in [tests/oracle/mpc_complex_oracle.h](tests/oracle/mpc_complex_oracle.h) | 20 random complex samples per function and type from `random_complex<Real>(-2, 1)`; log-family inputs are shifted away from the branch point by `1.5 + 0.25i`. | Uses the matching MPC functions. `log2` is computed as `mpc_log(z) / log(2)`. The result is compared component-wise in library eps against `complex_transcendental_bound()` (`4096` eps). | Diagnostics include MPC reference components, library limbs, component-wise absolute error, relative error in eps, seed, replay command, and MPC rounding mode. |
| MPC binary and parameterized functions | `pow(z,int)`, `pow(z,Real)`, `pow(Real,z)`, `pow(z,z)`, `ldexp(z,int)` | [tests/oracle/test_complex_binary.cpp](tests/oracle/test_complex_binary.cpp) | 20 random cases per function and type; bases are shifted by `1.25 + 0.125i`, exponents are bounded, and `ldexp` uses integer shifts. | Uses `mpc_pow`; integer powers use an explicit MPC exponentiation loop; `ldexp` applies `mpfr_mul_2si` to real and imaginary components. `ldexp` has an exact bound, while `pow` uses the complex transcendental bound. | Covers all public complex power overload shapes. |
| MPC/MPFR helper functions | `abs`, `norm`, `arg`, `polar`, component-wise `ceil`, `conj`, `proj` | [tests/oracle/test_complex_special.cpp](tests/oracle/test_complex_special.cpp) | Deterministic representative inputs per type. | Uses `mpc_abs`, MPFR `re^2+im^2`, MPFR `atan2`, MPFR sin/cos construction for `polar`, MPFR `ceil`, exact conjugation, and explicit finite/infinite projection checks. | `proj` verifies the infinite projection and signed zero on the imaginary component. |

## Common Commands

```sh
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure

cmake -S . -B build-mpfr -DBUILD_TESTING=ON -DQD3_ENABLE_MPFR_TESTS=ON
cmake --build build-mpfr -j
ctest --test-dir build-mpfr --output-on-failure

cmake -S . -B build-mpc -DBUILD_TESTING=ON -DQD3_ENABLE_MPC_TESTS=ON
cmake --build build-mpc -j
ctest --test-dir build-mpc --output-on-failure
```

For wider build-option coverage, use:

- [qa/check_16_builds_cmake.sh](qa/check_16_builds_cmake.sh)
- [qa/check_oracle_matrix_cmake.sh](qa/check_oracle_matrix_cmake.sh)
- [qa/check_complex_oracle_matrix_cmake.sh](qa/check_complex_oracle_matrix_cmake.sh)
