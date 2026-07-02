# libQD3 QA Notes

This file summarizes how the arithmetic QA covers real and complex addition,
subtraction, multiplication, and division.

## Real Arithmetic

Real arithmetic covers `dd_real`, `td_real`, `qd_real`, and `edd_real` when
`edd_real` is available.

| Operation | Primary oracle test | Input strategy | Reference and pass criterion | Extra coverage |
| --- | --- | --- | --- | --- |
| `a + b` | `tests/oracle/test_arith.cpp` / `oracle_test_arith_*` | 64 random same-sign or same-negated-sign expansion pairs per type, exponents in `[-220, 220]` | Convert both operands to MPFR, compute `mpfr_add(..., MPFR_RNDN)`, then require finite relative error in library eps within `add_sub_bound()` from `tests/oracle/fn_registry.h` | Exact add smoke with `1.25 + 2.5`; `+=` must be bit-identical to `a + b`; `T + double` is checked against MPFR; `test_rounding_corners.cpp` adds exact non-overlapping powers-of-two sum, signed-zero identities, tie-adjacent add, and add-variant cases across the CMake matrix. |
| `a - b` | `tests/oracle/test_arith.cpp` / `oracle_test_arith_*` | 64 random ordered positive pairs per type, then `b` is scaled by `1/4`; both operands may be negated together | Convert to MPFR, compute `mpfr_sub(..., MPFR_RNDN)`, then check relative error against `add_sub_bound()` | `-=` must be bit-identical to `a - b`; `test_rounding_corners.cpp` adds a Sterbenz exact subtraction case and signed-zero `a - a` behavior. |
| `a * b` | `tests/oracle/test_arith.cpp` / `oracle_test_arith_*` | 64 random expansion pairs per type, exponents in `[-140, 140]` | Convert to MPFR, compute `mpfr_mul(..., MPFR_RNDN)`, then check relative error against `mul_bound()` | `*=` must be bit-identical to `a * b`; `T * double` is checked against MPFR; `test_rounding_corners.cpp` adds exact power-of-two scaling, tie-adjacent multiplication, and qd sloppy/accurate multiplication variant cases. |
| `a / b` | `tests/oracle/test_arith.cpp` / `oracle_test_arith_*` | 64 random numerators with exponents in `[-220, 220]` and denominators in `[-80, 80]` per type | Convert to MPFR, compute `mpfr_div(..., MPFR_RNDN)`, then check relative error against `div_bound()` | `/=` must be bit-identical to `a / b`; default `qd_test` also contains division sanity and round-trip checks for the precision-specific smoke paths. |

The legacy `qd_test` executable is still part of the default CTest suite. It is
not the main oracle for the table above, but it exercises arithmetic through
identity formulas, fixed regressions, mixed-mode paths, and precision-specific
sanity checks.

## Complex Arithmetic

Complex arithmetic covers `dd_complex`, `td_complex`, `qd_complex`, and
`edd_complex` when `edd_real` is available. The default smoke test is
`tests/complex_test.cpp`; the optional MPC oracle tests are enabled with
`QD3_ENABLE_MPC_TESTS=ON`.

| Operation | Default smoke test | MPC oracle test | Input strategy | Reference and pass criterion |
| --- | --- | --- | --- | --- |
| `z + w` | `complex_test` checks `(1+2i) + (3-4i) == 4-2i`, compound `+=`, and real/complex mixed addition | `tests/oracle/test_complex_arith.cpp` / `oracle_test_complex_arith_*` uses `mpc_add(..., MPC_RNDNN)` | 24 random complex pairs per type from `random_complex<Real>(-2, 2)` | Convert real and imaginary parts through MPFR into MPC, compare the library result component-wise, and require `max(real_err, imag_err)` in eps to fit `complex_arithmetic_bound()` (`256` eps). |
| `z - w` | `complex_test` checks `(1+2i) - (3-4i) == -2+6i`, compound `-=`, and real/complex mixed subtraction | `oracle_test_complex_arith_*` uses `mpc_sub(..., MPC_RNDNN)` | 24 random complex pairs per type from `random_complex<Real>(-2, 2)` | Same component-wise MPC comparison and `256` eps arithmetic bound. |
| `z * w` | `complex_test` checks `(1+2i) * (3-4i) == 11+2i`, compound `*=`, and real/complex mixed multiplication | `oracle_test_complex_arith_*` uses `mpc_mul(..., MPC_RNDNN)` | 24 random complex pairs per type from `random_complex<Real>(-2, 2)` | Same component-wise MPC comparison and `256` eps arithmetic bound. |
| `z / w` | `complex_test` checks `(1+2i) / (3-4i)` against `-1/5 + 2/5 i`, compound `/=`, and real/complex mixed division | `oracle_test_complex_arith_*` uses `mpc_div(..., MPC_RNDNN)` | 24 random complex pairs per type; for division, the denominator is shifted by `1.25 - 0.75i` to avoid near-zero random denominators | Same component-wise MPC comparison and `256` eps arithmetic bound. |

Complex diagnostics include the seed, replay command, input limbs, MPFR/MPC
reference values, library output limbs, absolute error, relative error in eps,
and the active MPC rounding mode.

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

```sh
qa/check_16_builds_cmake.sh
qa/check_oracle_matrix_cmake.sh
qa/check_complex_oracle_matrix_cmake.sh
```
