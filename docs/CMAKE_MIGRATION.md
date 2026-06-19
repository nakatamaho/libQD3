# libQD3 CMake Build Notes

## Project Summary

libQD3 uses CMake as its build system. The CMake build covers the C++ and C
interfaces, optional Fortran interfaces, optional GNU x86 `_Float64x`
`edd_real` sources, CTest tests, install rules, pkg-config metadata,
`qd-config`, and an exported `QD3` CMake package.

## Core Options

| CMake option | Default | Notes |
| --- | --- | --- |
| `QD_ENABLE_INLINE` | `ON` | Install inline arithmetic definitions. |
| `QD_ENABLE_IEEE_ADD` | `OFF` | Select IEEE-style addition instead of the default sloppy add. |
| `QD_ENABLE_SLOPPY_MUL` | `ON` | Select sloppy QD multiplication where applicable. |
| `QD_ENABLE_SLOPPY_DIV` | `ON` | Select sloppy division where applicable. |
| `QD_ENABLE_BF` | `OFF` | Enable BF addition and multiplication. |
| `QD_ENABLE_BF_ADD` | `OFF` | Enable BF addition only. |
| `QD_ENABLE_BF_MUL` | `OFF` | Enable BF multiplication only. |
| `QD_FMA` | `auto` | `auto`, `yes`, `no`, `gnu`, `c99`, `ibm`, `ia64`, or `compiler`. |
| `QD_ARCH` | `generic` | `generic`, `x86-64-v3`, or `native`. |
| `QD_BUILD_FORTRAN` | `AUTO` | `AUTO`, `ON`, or `OFF`. |
| `QD_ENABLE_EDD_REAL` | `AUTO` | `AUTO`, `ON`, or `OFF`. |
| `BUILD_SHARED_LIBS` | `ON` | Build shared libraries. |
| `QD_BUILD_STATIC` | `ON` | Build static libraries. |
| `QD_PROPAGATE_FP_CONTRACT_FLAG` | `OFF` | Export the selected FP-contraction flag to consumers. |

## Probes

CMake generates a private build-tree `config.h` and a public build-tree
`include/qd/qd_config.h`. Only the public `qd_config.h` is installed. The
private header carries build-only macros such as `HAVE_FORTRAN`,
`HAVE_CLOCK_GETTIME`, `HAVE_GETTIMEOFDAY`, `HAVE_FPU_CONTROL_H`, `X86`, and
generated Fortran symbol wrappers.

CMake uses `check_cxx_source_runs()` for FMA candidates. `QD_FMA=auto` uses the
historical host preference order: PowerPC tries `ibm gnu c99`, IA-64 tries
`ia64 gnu c99`, and other targets try `gnu c99`. Explicit modes test only the
requested candidate and fail on mismatch.

Cross-compiling without `CMAKE_CROSSCOMPILING_EMULATOR` disables `auto` FMA
detection with a warning and rejects explicit FMA modes. Manual `QD_FMA_EXPR`
and `QD_FMS_EXPR` cache variables are available for toolchains where the
maintainer must assert the macro bodies.

The build probes `-fp-model strict` first, then `-ffp-contract=off`, and fails
if neither is accepted. The selected flag is applied to in-tree C++ library and
test targets. CMake does not add `-mfma` and does not add architecture flags
unless `QD_ARCH` requests them.

## edd_real Probe Behavior

`QD_ENABLE_EDD_REAL=AUTO` enables `edd_real` only when the GNU C++ x86
`_Float64x` binary80 compile probe succeeds. `ON` fails if the probe fails, and
`OFF` disables without probing. Headers are still installed in all cases.

## Fortran Support

Fortran is controlled by `QD_BUILD_FORTRAN=AUTO/ON/OFF`. When enabled, CMake
uses `FortranCInterface` to verify C++/Fortran compatibility and to generate
mangling macros used by `FC_FUNC` and `FC_FUNC_`. The `.f` module sources are
marked free-form, and shared/static Fortran libraries use separate module output
directories to avoid parallel-build races. Generated `.mod` files from the
primary Fortran library are installed to `${CMAKE_INSTALL_INCLUDEDIR}/qd`.

## QA Commands

```sh
cmake --preset default
cmake --build --preset default -j"$(nproc)"
ctest --preset default

cmake --preset ninja
cmake --build --preset ninja -j"$(nproc)"
ctest --preset ninja

bash qa/check_16_builds_cmake.sh
ctest -S qa/cmake_matrix.cmake -DSRC_DIR="$PWD" -DBUILD_ROOT="$PWD/_build_matrix_cmake"
bash qa/check_oracle_matrix_cmake.sh
```

Downstream smoke tests after installing to a prefix:

```sh
bash qa/consumer_smoke/run_cmake.sh /path/to/prefix
PKG_CONFIG_PATH=/path/to/prefix/lib/pkgconfig bash qa/consumer_smoke/run_pkg_config.sh /path/to/prefix
```

## Distribution

The `dist-cmake` target is a `git archive` convenience target for creating a
source archive from the current tree.
