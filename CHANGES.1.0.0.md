# Changes for libQD3 1.0.0

Commit range summarized:

- `98f6c05289382fdd1433c596f8854cdf026f19fe`
- through `HEAD`

This release starts libQD3 as an independently versioned fork of QD and
adds triple-double arithmetic as a first-class precision type alongside
`dd_real` and `qd_real`.

## Project Rename and Versioning

- Renamed the package from `qd` to `qd3`.
- Set the package version to `1.0.0`.
- Reset `PATCH_VERSION` to `0`.
- Updated package metadata and maintainer contact for independent libQD3
  development.

## New `td_real` Triple-Double Type

- Added native `td_real`, a three-double floating-point type with normalized
  three-word storage.
- Added public headers, source files, constants, build integration, and test
  support for `td_real`.
- Implemented core `td_real` arithmetic:
  addition, subtraction, multiplication, division, unary minus, `sqr`, `sqrt`,
  `abs`, `fabs`, comparisons, stream I/O, parsing, and conversion helpers.
- Added `td_real` constants including:
  `_nan`, `_inf`, `_eps`, `_min_normalized`, `_max`, `_safe_max`, `_ndigits`,
  `_2pi`, `_pi`, `_3pi4`, `_pi2`, `_pi4`, `_e`, `_log2`, and `_log10`.
- Moved `td_real` constants into `src/td_const.cpp`, matching the existing
  `dd_const.cpp` and `qd_const.cpp` organization.
- Added `td_inline.h` and reorganized the triple-double inline layer to mirror
  the structure used by `dd_inline.h` and `qd_inline.h`.

## Native Triple-Double Math

- Added practical `td_real` transcendental support and then replaced the
  temporary qd-backed fallback implementation with native triple-double
  algorithms.
- Native `td_real` functions now include:
  `exp`, `log`, `log10`, `sin`, `cos`, `sincos`, `tan`, `atan`, `atan2`,
  `asin`, `acos`, `sinh`, `cosh`, `sincosh`, `tanh`, `asinh`, `acosh`,
  `atanh`, and `pow`.
- Added native `td_real` `nroot` and `nint`.
- Aligned `td_real` division with the qd-style long-division design:
  accurate division extracts an additional guard digit, while sloppy division
  uses fewer quotient digits but still updates residuals in full triple-double
  arithmetic.
- Added sloppy add and sloppy division variants for `td_real`.

## Mixed-Precision C++ API

- Expanded mixed-mode arithmetic and comparison coverage among:
  `double`, `dd_real`, `td_real`, and `qd_real`.
- Added missing mixed-mode self operators and non-member operators involving
  `td_real`.
- Added conversions between `dd_real`, `td_real`, and `qd_real`.
- Added normalized comparison semantics for `dd_real` and `qd_real`:
  values are compared by exact normalized differences rather than direct
  component layout, and relational operators follow the same semantics.
  NaN and infinity cases are handled before subtraction.

## C API

- Added and completed C wrapper coverage for `td_real`.
- Added cross-type C wrapper support across `dd`, `td`, and `qd` for:
  arithmetic, self operations, copy helpers, and comparison helpers.
- Extended C tests to cover the new mixed-type wrapper paths.

## Fortran API

- Added native Fortran support for `td_real`.
- Added Fortran wrapper and module files for triple-double support:
  `fortran/f_td.cpp`, `fortran/tdext.f`, and `fortran/tdmod.f`.
- Added Fortran `tdmodule` functionality for:
  constructors, constants, inquiry helpers, arithmetic, comparisons,
  elementary functions, rounding functions, input/output helpers, and
  conversions among `dd_real`, `td_real`, and `qd_real`.
- Added `td_complex` to the Fortran module.
- Added `td_complex` arithmetic, constructors, conversions, I/O helpers, and
  utilities such as `conjg`, `aimag`, `abs`, `exp`, and `log`.
- Added Fortran cross-type complex conversions among:
  `td_complex`, `dd_complex`, and `qd_complex`.

## Documentation

- Added `docs/td.tex` and integrated triple-double documentation into the
  documentation build.
- Added `td.pdf` targets to `docs/Makefile` and top-level documentation rules.
- Added binary64 flop-count analysis for DD, TD, and QD arithmetic.
- Updated TD division documentation to match the implemented long-division
  design and cost model.
- Added 2-clause BSD headers to newly introduced triple-double source files.

## Tests and QA

- Extended `tests/qd_test.cpp` with `-td` support and a broad set of
  triple-double tests.
- Added tests for `td_real` normalization, cancellation-sensitive arithmetic,
  mixed-magnitude operations, division, square roots, parsing, formatting,
  mixed-mode arithmetic, conversions, special values, and transcendental
  functions.
- Extended `tests/huge.cpp` to cover `td_real` alongside `dd_real` and
  `qd_real`.
- Promoted `huge` into the default autotools `make check` test suite.
- Added `qa/check_16_builds.sh`, which validates the release matrix across
  all 16 combinations of:
  `ieee_add`, `sloppy_mul`, `sloppy_div`, and `fma`, plus the default
  configuration.
- The QA script performs out-of-tree builds, runs `make check` and
  `tests/huge`, and writes per-build logs and summary files.
- The QA script copies `docs/` into each out-of-tree build directory so the
  top-level documentation targets work when `make` builds `docs/qd.pdf` and
  `docs/td.pdf`.

## Test Maintenance

- Fixed `tests/qd_test.cpp` test numbering for the td-specific tests.
- Replaced several raw dynamic arrays in tests with RAII containers.
- Added an RAII guard for temporary suppression of td error messages.
- Fixed a loop counter type and cleaned up misleading local variable names in
  tests.
- Stripped trailing whitespace in touched test code.

