# libQD3 CMake Migration Notes

## Project Summary

libQD3 is the qd3 1.1.0 fork of the QD extended-precision library. The CMake build added in this phase builds the C++ and C interfaces, optional Fortran interfaces, optional GNU x86 `_Float64x` `edd_real` sources, CTest tests, install rules, pkg-config metadata, `qd-config`, and an exported `QD3` CMake package.

## Build-System Coexistence Policy

Autotools remains authoritative during this phase. CMake is implemented side-by-side and no Autotools files are removed. The CMake build rejects in-source builds and is intended for Unix Makefiles and Ninja.

## CMake Option Parity

| Autotools | CMake | Default |
| --- | --- | --- |
| `--enable-inline` | `QD_ENABLE_INLINE` | `ON` |
| `--enable-ieee-add` | `QD_ENABLE_IEEE_ADD` | `OFF` |
| `--enable-sloppy-mul` | `QD_ENABLE_SLOPPY_MUL` | `ON` |
| `--enable-sloppy-div` | `QD_ENABLE_SLOPPY_DIV` | `ON` |
| `--enable-fma=...` | `QD_FMA=auto/yes/no/gnu/c99/ibm/ia64/compiler` | `auto` |
| `--with-arch=...` | `QD_ARCH=generic/x86-64-v3/native` | `generic` |
| `--enable-fortran` | `QD_BUILD_FORTRAN=AUTO/ON/OFF` | `AUTO` |
| implicit probe | `QD_ENABLE_EDD_REAL=AUTO/ON/OFF` | `AUTO` |
| `--enable-shared` | `BUILD_SHARED_LIBS` | `ON` |
| `--enable-static` | `QD_BUILD_STATIC` | `ON` |
| none | `QD_PROPAGATE_FP_CONTRACT_FLAG` | `OFF` |

## Configure-Probe Parity

| Autotools probe | CMake handling | Outcome |
| --- | --- | --- |
| `AC_CHECK_HEADERS([ieeefp.h])` | `check_include_file_cxx(ieeefp.h)` | private `HAVE_IEEEFP_H` |
| x86 host case | `CMAKE_SYSTEM_PROCESSOR` check | private `X86` and gated FPU probe |
| `AC_CHECK_HEADERS([fpu_control.h])` | `check_include_file_cxx(fpu_control.h)` | private `HAVE_FPU_CONTROL_H` |
| `AC_CHECK_FUNCS([gettimeofday])` | `check_symbol_exists(gettimeofday sys/time.h)` | private `HAVE_GETTIMEOFDAY` |
| `AX_CXX_CLOCK_GETTIME` / `AC_SEARCH_LIBS` | libc symbol check, then `rt` library check | private `HAVE_CLOCK_GETTIME` and optional `rt` link |
| `AC_CHECK_LIB(m,sqrt)` | `check_library_exists(m sqrt ...)` | optional `m` link |
| `AX_CXX_FMA` | run-test candidate checks | public/private `QD_FMA` and `QD_FMS` macros |
| `_Float64x` `edd_real` probe | compile test matching Autotools intent | optional `QD_HAVE_EDD_REAL` macros and sources |
| `AC_FC_WRAPPERS` | `FortranCInterface` | generated `FC_FUNC` / `FC_FUNC_` bridge |
| `AC_HEADER_STDBOOL` | documented as obsolete | no consumed generated macro in current sources |
| `AC_STRUCT_TM` | documented as obsolete | no consumed generated macro in current sources |

The old `AX_CXX_ISNAN`, `AX_CXX_ISINF`, `AX_CXX_ISFINITE`, and `AX_CXX_COPYSIGN` probes are simplified in CMake: C++ code uses `std::isnan`, `std::isinf`, `std::isfinite`, and `std::copysign`; C consumers get the C `<math.h>` forms.

## FMA Detection Behavior

CMake uses `check_cxx_source_runs()` for FMA candidates. `QD_FMA=auto` uses the Autotools host preference order: PowerPC tries `ibm gnu c99`, IA-64 tries `ia64 gnu c99`, and other targets try `gnu c99`. `QD_FMA=yes` uses the historical default candidate list `ibm gnu c99 compiler` and fails if none work. Explicit modes test only the requested candidate and fail on mismatch.

Cross-compiling without `CMAKE_CROSSCOMPILING_EMULATOR` disables `auto` FMA detection with a warning and rejects explicit FMA modes. Manual `QD_FMA_EXPR` and `QD_FMS_EXPR` cache variables are available for toolchains where the maintainer must assert the macro bodies.

## Cross-Compilation Behavior

CMake does not run FMA probes while cross-compiling unless `CMAKE_CROSSCOMPILING_EMULATOR` is configured. `QD_ARCH=native` is rejected when cross-compiling. `QD_ARCH=x86-64-v3` is validated against `CMAKE_SYSTEM_PROCESSOR` and by a compiler flag probe. Manual FMA macro overrides are the supported escape hatch for cross toolchains that cannot execute probe binaries.

## FP-Contraction Policy

The CMake build probes `-fp-model strict` first, then `-ffp-contract=off`, and fails if neither is accepted. The selected flag is applied to in-tree C++ library and test targets. CMake does not add `-mfma` and does not add architecture flags unless `QD_ARCH` requests them.

`QD_PROPAGATE_FP_CONTRACT_FLAG=ON` adds the selected flag to exported C++ target usage requirements. The default is `OFF` for Autotools parity, but propagation may be safer for consumers when `QD_INLINE=1` causes EFT kernels to be compiled in downstream translation units.

## edd_real Probe Behavior

`QD_ENABLE_EDD_REAL=AUTO` enables `edd_real` only when the GNU C++ x86 `_Float64x` binary80 compile probe succeeds. `ON` fails if the probe fails, and `OFF` disables without probing. Headers are still installed in all cases, matching Autotools header installation behavior.

## Generated Headers

CMake generates a private build-tree `config.h` and a public build-tree `include/qd/qd_config.h` from separate CMake templates. Only the public `qd_config.h` is installed. The private header carries build-only macros such as `HAVE_FORTRAN`, `HAVE_CLOCK_GETTIME`, `HAVE_GETTIMEOFDAY`, `HAVE_FPU_CONTROL_H`, `X86`, and generated Fortran symbol wrappers.

## Fortran Support

Fortran is controlled by `QD_BUILD_FORTRAN=AUTO/ON/OFF`. When enabled, CMake uses `FortranCInterface` to verify C++/Fortran compatibility and to generate mangling macros used by `FC_FUNC` and `FC_FUNC_`. The `.f` module sources are marked free-form, and shared/static Fortran libraries use separate module output directories to avoid parallel-build races. Generated `.mod` files from the primary Fortran library are installed to `${CMAKE_INSTALL_INCLUDEDIR}/qd`.

CMake builds and tests the required `f_test`. The historical Fortran demo/timer programs are not wired into the CMake `demo` or `time` targets in this phase because several demo sources generate the same module names and require per-target module directories to avoid races. `fortran/second.f` is generated when Fortran is enabled so that adding those programs later has the Autotools substitution point available.

## Known Install Differences

CMake installs `QD3Config.cmake`, `QD3ConfigVersion.cmake`, and `QD3Targets.cmake` under `${CMAKE_INSTALL_LIBDIR}/cmake/QD3`. Autotools installs libtool `.la` files; CMake does not. CMake installs Fortran `.mod` files under `include/qd` as required by this migration; the current Autotools `pkginclude` path installs them under `include/qd3`. The CMake-generated `qd-config` omits build-tree-only Autotools modes such as `--src`, `--libs-la`, `--build-flags`, `--build-libs`, and `--configure-args`.

## QA Commands

```sh
autoreconf -fi
./configure
make -j"$(nproc)"
make check

cmake -S . -B build
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure

cmake -S . -B build-ninja -G Ninja
cmake --build build-ninja -j"$(nproc)"
ctest --test-dir build-ninja --output-on-failure

bash qa/check_16_builds_cmake.sh
bash qa/compare_install_trees.sh
bash qa/check_version_sync.sh
```

Downstream smoke tests after installing to a prefix:

```sh
bash qa/consumer_smoke/run_cmake.sh /path/to/prefix
PKG_CONFIG_PATH=/path/to/prefix/lib/pkgconfig bash qa/consumer_smoke/run_pkg_config.sh /path/to/prefix
```

## Known Non-Goals

This phase does not remove Autotools, rewrite numerical algorithms, modernize the source, add `-ffast-math`, silently enable hardware ISA flags, implement MSVC/Windows support, or replace `make dist` completely. The `dist-cmake` target is only a `git archive` convenience target.

## Future Autotools Removal Checklist

1. Run the full Autotools and CMake QA matrices on GCC, Clang, Intel oneAPI, Ninja, and Unix Makefiles.
2. Verify install tree compatibility on lib and lib64 platforms.
3. Verify downstream CMake, pkg-config, and `qd-config` consumers.
4. Decide whether `QD_PROPAGATE_FP_CONTRACT_FLAG` should become the default.
5. Add race-free CMake wiring for historical Fortran demos/timers if they remain supported.
6. Only then remove Autotools files in a separate migration phase.
