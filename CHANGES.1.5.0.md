# Changes for libQD3 1.5.0

This release adds binary32-based float expansion arithmetic and extends the
C++ and Fortran QA and packaging coverage.

## Float expansion types

- Added `ds_real`, `ts_real`, and `qs_real`, using two, three, and four
  binary32 limbs respectively.
- Added public C++ headers and support for construction, conversion,
  arithmetic, comparisons, elementary functions, special values,
  parsing/formatting, random generation, and `numeric_limits` for the new
  types.
- Added Fortran modules `dsmodule`, `tsmodule`, and `qsmodule`, including the
  `ds_real`, `ts_real`, and `qs_real` types, their complex counterparts,
  generic arithmetic and conversion interfaces, elementary functions, I/O,
  comparisons, and random-number interfaces.
- Added the corresponding low-level Fortran declarations and C++ wrappers,
  and integrated installation of the generated Fortran modules with CMake.

## Tests and build integration

- Added dependency-free `regression_smoke` coverage for special-value,
  overflow, division, and square-root regression paths.
- Added `single_smoke` and `single_test` coverage for the new float expansion
  types, including arithmetic, elementary functions, special values,
  formatting, and state behavior.
- Added seeded MPFR oracle coverage for DS, TS, and QS, with replay
  diagnostics and dedicated CTest entries:
  `oracle_test_single_ds`, `oracle_test_single_ts`, and
  `oracle_test_single_qs`.
- Extended the Fortran test program to exercise DS, TS, and QS arithmetic,
  elementary functions, rounding, comparisons, metadata, complex helpers,
  and random-number interfaces.

## Verification

- Default CMake/Fortran CTest suite: 14/14 tests passed.
- MPFR-enabled CTest suite: 52/52 tests passed.
- CMake installation of all six new Fortran modules was verified.
- Legacy Autotools entry points remain absent; the branch continues to use
  the CMake-only build flow.
