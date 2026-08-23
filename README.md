# libQD3

libQD3 is a fork of the [QD library](https://github.com/BL-highprecision/QD).
It keeps the original double-double and quad-double functionality and adds
native extended-double and triple-double support.

libQD3 provides extended-precision floating-point types for C++, C, and
Fortran:

- `dd_real`: double-double precision, about 106 mantissa bits
- `edd_real`: extended-double precision based on two binary80 `_Float64x`
  limbs
- `td_real`: triple-double precision, about 159 mantissa bits
- `qd_real`: quad-double precision, about 212 mantissa bits

libQD3 also provides matching C++ complex aliases through `<qd/complex.h>`:
`dd_complex`, `td_complex`, `qd_complex`, and `edd_complex` when `edd_real` is
available. These overload arithmetic, standard complex elementary functions,
`proj`, `ldexp`, and component-wise `ceil` through normal unqualified lookup
plus ADL; no overloads are added to namespace `std`.

## News
2026-08-23 libQD3 1.4.0 was released.  This numerical robustness
release hardens division rescaling and IEEE special-value handling for DD, TD,
QD, and EDD, adds square-root special-value handling, fixes the EDD maximum
constant on targets without a `_Float64x` `numeric_limits` specialization, and
adds focused arithmetic smoke coverage.  See
[CHANGES.1.4.0.md](CHANGES.1.4.0.md) for the release notes.

2026-08-10 libQD3 1.3.1 was released.  This focused portability release fixes
x87/i386 FPU-mode separation for EDD: EDD arithmetic stays in verified 80-bit
mode, while DD/TD/QD oracle operations and the internal `qd_real` path used by
EDD trigonometric argument reduction run in QD round-to-double mode.  See
[CHANGES.1.3.1.md](CHANGES.1.3.1.md) for the 1.3.1 release notes.

2026-08-07 libQD3 1.3.0 was released.  This maintenance release folds
MPLAPACK downstream portability patches into upstream: integer mixed-mode
operators for DD/TD/QD/EDD real types, std::complex<double> interop for the
C++ complex wrappers, CMake uninstall support, stricter FMA probing, x87
80-bit FPU handling for EDD tests, shared-library SOVERSION 2, and MinGW/Wine
EDD notes.  See [CHANGES.1.3.0.md](CHANGES.1.3.0.md) for the 1.3.0 release
notes.

2026-07-03 libQD3 1.2.0 was released.  This release adds first-class
extended-precision complex C++ types, a CMake-only build, additional real C++
math overloads, optional branch-free arithmetic, improved TD/EDD
trigonometric paths, and MPFR/MPC oracle QA.  See
[CHANGES.1.2.0.md](CHANGES.1.2.0.md) for the 1.2.0 release notes.

2026-04-22 libQD3 1.1.0 was released.  This release adds
binary80-based `edd_real` extended-double support, including core
arithmetic, constants, conversions, transcendental functions, a small C API,
documentation, and AGM examples.  See
[CHANGES.1.1.0.md](CHANGES.1.1.0.md) for the 1.1.0 release notes.

2026-04-20 libQD3 1.0.0 was released.  This first release adds
TD/triple-double support through the new `td_real` type.  See
[CHANGES.1.0.0.md](CHANGES.1.0.0.md) for the 1.0.0 release notes.

## Citation
Y. Hida, X. S. Li and D. H. Bailey, "Algorithms for quad-double precision floating point arithmetic," Proceedings 15th IEEE Symposium on Computer Arithmetic. ARITH-15 2001, Vail, CO, USA, 2001, pp. 155-162, doi: 10.1109/ARITH.2001.930115.

## End-user documentation

For historical reasons, the original user-focused documentation is
located in the `README` file, sans the `.md` extension. In-depth
technical documentation is also provided in the `docs` directory.
Release tarballs include the generated PDFs.

To build the PDFs from LaTeX sources, run

```
$ make -C docs qd.pdf
$ make -C docs edd.pdf
$ make -C docs td.pdf
```

after installing the necessary LaTeX bits on your system.

Release-specific notes for libQD3 1.4.0 are in
[CHANGES.1.4.0.md](CHANGES.1.4.0.md).

## Tips for developers

libQD3 uses CMake as its build system. A normal developer build is:

```
$ cmake -S . -B build
$ cmake --build build -j
$ ctest --test-dir build --output-on-failure
```

Build knobs are passed at configure time with `-DNAME=VALUE`, for example:

```
$ cmake -S . -B build-bf-mul -DQD_ENABLE_BF_MUL=ON
$ cmake --build build-bf-mul -j
```

Common build knobs:

| Option | Default | Purpose |
| --- | --- | --- |
| `QD_ENABLE_INLINE` | `ON` | Install inline definitions for commonly used functions. |
| `QD_ENABLE_IEEE_ADD` | `OFF` | Use IEEE-style addition instead of the default sloppy addition. |
| `QD_ENABLE_SLOPPY_MUL` | `ON` | Use the existing faster sloppy multiplication path where applicable. |
| `QD_ENABLE_SLOPPY_DIV` | `ON` | Use the existing faster sloppy division path where applicable. |
| `QD_ENABLE_BF` | `OFF` | Enable both branch-free addition and multiplication. |
| `QD_ENABLE_BF_ADD` | `OFF` | Enable only branch-free addition. |
| `QD_ENABLE_BF_MUL` | `OFF` | Enable only branch-free multiplication. |
| `QD_FMA` | `auto` | Select FMA implementation: `auto`, `yes`, `no`, `gnu`, `c99`, `ibm`, `ia64`, or `compiler`. |
| `QD_ARCH` | `generic` | Select code generation target: `generic`, `x86-64-v3`, or `native`. |
| `QD_BUILD_FORTRAN` | `AUTO` | Build Fortran interfaces: `AUTO`, `ON`, or `OFF`. |
| `QD_ENABLE_EDD_REAL` | `AUTO` | Build binary80 `_Float64x` `edd_real` sources: `AUTO`, `ON`, or `OFF`. |
| `BUILD_SHARED_LIBS` | `ON` | Build shared libraries. |
| `QD_BUILD_STATIC` | `ON` | Build static libraries. |
| `QD_PROPAGATE_FP_CONTRACT_FLAG` | `OFF` | Export the selected FP-contraction flag to downstream C++ consumers. |
| `QD3_ENABLE_MPFR_TESTS` | `OFF` | Build the optional MPFR oracle test suite. |
| `QD3_ENABLE_MPC_TESTS` | `OFF` | Build the optional MPC complex oracle test suite. |
| `QD3_ENABLE_ASAN` | `OFF` | Build QA targets with AddressSanitizer. |
| `QD3_ENABLE_UBSAN` | `OFF` | Build QA targets with UndefinedBehaviorSanitizer. |
| `QD3_ENABLE_COVERAGE` | `OFF` | Build QA targets with coverage instrumentation and add the `coverage` target. |

`QD_ENABLE_BF=ON` turns on both `QD_ENABLE_BF_ADD` and `QD_ENABLE_BF_MUL`. At
least one of `BUILD_SHARED_LIBS` or `QD_BUILD_STATIC` must be `ON`.

### Optional branch-free arithmetic

Branch-free addition and multiplication for `dd_real`, `td_real`, and `qd_real`
can be enabled at compile time. Use `-DQD_ENABLE_BF_MUL=ON` for the scalar
TD/QD multiply speedup; use `-DQD_ENABLE_BF=ON` to enable both BF add and BF
mul for SIMD-oriented or deterministic dataflow builds. BF is off by default.

### Running the test suite

The libQD3 library comes with an automated CTest suite that should always pass:

```
$ cmake -S . -B build
$ cmake --build build -j
$ ctest --test-dir build --output-on-failure
```

The default suite also builds `complex_test`, which covers the public
`dd_complex`, `td_complex`, `qd_complex`, and optional `edd_complex` headers.

### Known issue: MinGW/Wine `edd_real` trigonometry

`edd_real` is implemented on native binary80 limbs.  On MinGW targets executed
under Wine, the binary80 `long double` trigonometric and inverse-trigonometric
runtime path may not satisfy libQD3's EDD test tolerances, especially for
large-argument reduction.  In observed `x86_64-w64-mingw32-g++` builds under
Wine, the EDD-only parts of the default suite failed in:

- `qd_test -edd`: Test 15, EDD trig identities and qd oracle agreement
- `qd_test -edd`: Test 16, EDD inverse trig consistency
- `qd_test -edd`: Test 17, EDD large-argument reduction and hyperbolic checks
- `complex_test`: `edd_complex`

The `dd_real`, `td_real`, and `qd_real` tests passed in the same environment.
This is tracked as a libQD3 platform/runtime issue, not as a downstream
consumer issue.  Release builds that only require DD/TD/QD functionality may
disable EDD on MinGW with:

```
$ cmake -S . -B build-mingw -DQD_ENABLE_EDD_REAL=OFF
```

The optional MPC complex oracle suite can be enabled independently from the
MPFR real oracle suite:

```
cmake -S . -B build-mpc -DBUILD_TESTING=ON -DQD3_ENABLE_MPC_TESTS=ON
cmake --build build-mpc -j
ctest --test-dir build-mpc --output-on-failure
```

The main numeric test can also be run by precision after building:

```
$ ./build/tests/qd_test -dd
$ ./build/tests/qd_test -edd
$ ./build/tests/qd_test -td
$ ./build/tests/qd_test -qd
```

For release-oriented configuration coverage, run:

```
$ bash qa/check_16_builds_cmake.sh
```

Additional QA-only builds are opt-in and do not affect the default test path:

```
$ cmake -S . -B build-qa -DQD3_ENABLE_ASAN=ON -DQD3_ENABLE_UBSAN=ON
$ cmake --build build-qa -j
$ ctest --test-dir build-qa --output-on-failure

$ cmake -S . -B build-coverage -DQD3_ENABLE_COVERAGE=ON
$ cmake --build build-coverage --target coverage
```

The optional MPFR oracle suite is guarded by `-DQD3_ENABLE_MPFR_TESTS=ON`.
When enabled, configuration fails immediately if MPFR/GMP cannot be found. The
oracle build matrix runner is:

```
$ bash qa/check_oracle_matrix_cmake.sh
$ ENABLE_MPFR_ORACLE=ON bash qa/check_16_builds_cmake.sh
```

Both commands export a deterministic `QD_TEST_SEED` by default and store logs
under `_build_matrix_cmake/` unless `BUILD_ROOT` is set.

### Making a release

There are several steps that need to be performed when making a new release:

1. Ensure that all important user-facing changes are mentioned in the `NEWS`
   file and the appropriate `CHANGES.*.md` file.

2. Update the package version number in `CMakeLists.txt`.

3. Build and test from a clean tree:

```
$ git clean -x -f -d
$ cmake -S . -B build
$ cmake --build build -j
$ ctest --test-dir build --output-on-failure
$ bash qa/check_16_builds_cmake.sh
```

4. Create the source archive from the tagged commit. The CMake build provides a
   `dist-cmake` convenience target based on `git archive`.

5. Tag the commit that corresponds to the release with `git tag -s
   <version-number>`.

6. Push everything to GitHub and upload the generated archive to the GitHub
   release page.
