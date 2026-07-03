# Changes for libQD3 1.2.0

This file summarizes all commits after tag `1.1.0` through `HEAD` for the
1.2.0 release.

## User-Facing Enhancements

### Extended-Precision Complex C++ API

- Added first-class public complex headers:
  - `include/qd/detail/complex_impl.h`
  - `include/qd/dd_complex.h`
  - `include/qd/td_complex.h`
  - `include/qd/qd_complex.h`
  - `include/qd/edd_complex.h`
  - `include/qd/complex.h`
- Added `dd_complex`, `td_complex`, `qd_complex`, and `edd_complex` when
  `edd_real` is available.
- Implemented the complex types with one shared `qd3_complex<Real>` template
  body instead of separate copied implementations.
- Added construction, real/imag accessors, arithmetic, mixed real/complex
  operations, conversion from `std::complex<double>`, and conversion back to
  `std::complex<double>`.
- Added unqualified-call and ADL-friendly overloads for the standard complex
  function set without adding overloads to namespace `std`.
- Complex functions now cover `conj`, `abs`, `norm`, `arg`, `polar`, `proj`,
  `sqrt`, `exp`, `log`, `log10`, `log2`, `sin`, `cos`, `tan`, `asin`, `acos`,
  `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `pow`, `ldexp`,
  and component-wise `ceil`.

### Real C++ Math API Coverage

- Added missing practical C++17-style real math overloads for `dd_real`,
  `td_real`, `qd_real`, and `edd_real`:
  `log2`, `exp2`, `expm1`, `log1p`, `hypot`, `cbrt`, `trunc`, and `round`.
- Closed related type gaps, including TD/EDD `fmod` coverage, `edd_real`
  `pow(real, real)`, and missing EDD inverse-hyperbolic C++ declarations and
  definitions.
- Added exact `std::int64_t` and `std::uint64_t` constructors and assignment
  operators for `dd_real`, `td_real`, `qd_real`, and `edd_real`, avoiding silent
  rounding through `double` for values such as `UINT64_MAX` and `2^53 + 1`.

### Accuracy, Portability, and Compatibility Fixes

- Fixed an `O(eps^3)` cross-term accumulation bug in mixed
  `qd_real * dd_real` multiplication.
- Fixed fixed-format zero-precision output so requested field width is
  preserved.
- Reworked `td_real` trigonometric argument reduction to use a native
  four-limb residual helper shared by production code and oracle diagnostics.
- Added deterministic TD trigonometric worst-case regressions for the improved
  reduction path.
- Improved `edd_real` transcendental and trigonometric paths, including large
  finite argument handling and backend-specific low-level math wrappers.
- Enabled `edd_real` on binary80 `long double` targets when `_Float64x`
  builtins are unavailable, matching MinGW-style binary80 environments.
- Defined EDD constants directly as normalized binary80 limbs to avoid static
  initialization order dependencies on `qd_real` constants.
- Fixed the EDD C API oracle bridge.

### Build and Packaging

- Added a first-class CMake build system with install rules, package config,
  pkg-config metadata, CMake presets, and build-option normalization.
- Made this branch CMake-only and removed legacy Autotools entry points and
  metadata from the release path.
- Added documented CMake build knobs for inline mode, arithmetic variants,
  FMA selection, static/shared libraries, Fortran support, MPFR tests, MPC
  tests, sanitizer builds, and coverage builds.
- Added consumer smoke tests for installed CMake package and pkg-config usage.

### Optional Branch-Free Arithmetic

- Added opt-in branch-free addition and multiplication paths for `dd_real`,
  `td_real`, and `qd_real`, based on the Kouya/Zhang-Aiken algorithms.
- The feature is off by default and can be enabled with `QD_BF`, or more
  selectively with `QD_BF_ADD` and `QD_BF_MUL`.
- CMake options are `-DQD_ENABLE_BF=ON`, `-DQD_ENABLE_BF_ADD=ON`, and
  `-DQD_ENABLE_BF_MUL=ON`.
- The default arithmetic behavior remains unchanged when branch-free options
  are not enabled.

### Documentation and Examples

- Added `examples/edd_euler_gamma.cpp`.
- Refreshed README and release documentation for EDD support, CMake usage,
  complex support, new real math overloads, and optional oracle testing.
- Moved the TODO list to Markdown.

## QA Enhancements

### MPFR Real Oracle Tests

- Added an optional MPFR-backed oracle suite under `tests/oracle`.
- Real oracle coverage includes arithmetic, unary functions, binary and
  parameterized functions, composed identities, deterministic endpoint and
  special-value contracts, I/O, and C API checks.
- Added deterministic rounding-corner tests and matrix comparisons for
  arithmetic build variants.
- Added coverage for the new real math overloads and for the TD/EDD API gaps
  closed in this release.

### MPC Complex Oracle Tests

- Added an optional MPC + MPFR complex oracle suite controlled by
  `QD3_ENABLE_MPC_TESTS`.
- MPC tests are independent from MPFR-only real oracle tests; enabling the
  existing MPFR real oracle path does not require MPC.
- Complex oracle coverage includes arithmetic, unary elementary functions,
  binary/parameterized functions, and helpers such as `abs`, `norm`, `arg`,
  `polar`, `ceil`, `conj`, and `proj`.
- Made the MPC complex oracle portable across MinGW by using
  `pkg-config --define-prefix` during cross builds and by linking with
  pkg-config library names, include directories, and link directories.
- Changed complex oracle error measurement to a norm-wise comparison to avoid
  false failures from cancellation in a single component.

### Reproducible Diagnostics

- Added seeded TAP diagnostics with reproducible replay commands.
- Oracle seed resolution is `QD_TEST_SEED`, then `--seed=<N>`, then the fixed
  default `0x9E3779B97F4A7C15`.
- Added verbose MPFR/MPC diagnostics for passing cases, including input limbs,
  reference values, library values, absolute error, relative error in eps, ULP
  estimates, bounds, and replay commands.
- Added a TAP legacy log reformatter and deterministic special-test diagnostic
  behavior.
- Added worst-case unary reports, TD trigonometric reduction sector
  diagnostics, and comparison tooling for precision snapshots.

### Build Matrix and Release Gate QA

- Added CMake matrix scripts for the default smoke matrix, MPFR oracle matrix,
  and MPC complex oracle gate:
  - `qa/check_16_builds_cmake.sh`
  - `qa/check_oracle_matrix_cmake.sh`
  - `qa/check_complex_oracle_matrix_cmake.sh`
- Made the CMake FMA matrix host-aware. The matrix now probes whether an
  x86_64/amd64 host can actually run x86-64-v3/FMA code before selecting that
  configuration, and non-x86 hosts use an appropriate generic configuration.
- Added sanitizer and coverage workflows for the oracle suite.
- Refreshed coverage filtering to include compiled `src/*.cpp`,
  `src/td_trig_reduce.h`, and `include/qd/*_inline.h`.
- Added `QA.md`, with linked tables describing how real arithmetic, real
  special functions, complex arithmetic, and complex special functions are
  tested and why the oracle strategy is a practical correctness argument.

## Release Test Matrix

Release-oriented testing for 1.2.0 is CMake-only and MPC-inclusive. Each
release candidate should be tested from a clean tree on these platforms:

| Platform | Required release gate |
| --- | --- |
| Linux, Ubuntu 26.04, amd64 | Clean CMake smoke matrix, MPFR oracle matrix, and MPC complex oracle gate. |
| macOS, arm64 | Clean CMake smoke matrix, MPFR oracle matrix, and MPC complex oracle gate. |
| macOS, amd64 | Clean CMake smoke matrix, MPFR oracle matrix, and MPC complex oracle gate. |
| Windows, MinGW | Clean CMake/MinGW build and CTest run with MPC complex oracle tests enabled. |

The release gate is:

```sh
git status --short --branch
test ! -e configure && test ! -e configure.ac && test ! -e Makefile.am && test ! -e autogen.sh
git clean -x -f -d

bash qa/check_16_builds_cmake.sh
bash qa/check_oracle_matrix_cmake.sh
bash qa/check_complex_oracle_matrix_cmake.sh
```

For Windows MinGW, the same source tree is configured with the MinGW C and C++
compilers, MPFR/MPC/GMP from the MinGW prefix, and `QD3_ENABLE_MPC_TESTS=ON`;
the resulting test executables are run through the MinGW runtime or Wine.

## Scope Notes

- No BLAS/LAPACK complex wrappers were added.
- No complex C API was added.
- No Fortran complex bindings were added.
- No FFT support was added.
- No `std::complex<dd_real>` specialization was added.
