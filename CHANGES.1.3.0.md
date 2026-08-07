# Changes for libQD3 1.3.0

This file summarizes the changes after tag `1.2.0` through the 1.3.0 release.
The release primarily upstreams the libQD3 patches carried by MPLAPACK.

## C++ API Compatibility

- Added `std::complex<double>` mixed arithmetic operators for libQD3 complex
  wrapper types.  This makes generic complex code work without relying on an
  ambiguous conversion path through `double`.
- Added integer mixed-mode arithmetic overloads for `dd_real`, `td_real`,
  `qd_real`, and `edd_real`: `+`, `-`, `*`, and `/` with supported integer
  operands on either side.
- Added integer comparison overloads for `dd_real`, `td_real`, `qd_real`, and
  `edd_real`: `==`, `!=`, `<`, `>`, `<=`, and `>=` with supported integer
  operands on either side.
- Added explicit integer conversion helpers such as `to_long`,
  `to_unsigned_long`, `to_long_long`, `to_unsigned_long_long`, `to_int64_t`,
  and `to_uint64_t`.
- Kept the existing integral constructor and assignment template work from
  post-1.2.0 `main`, which avoids routing supported integer inputs through
  `double`.

## Build And Packaging

- Added a CMake `uninstall` target based on the generated install manifest.
- Changed shared-library ABI names from SOVERSION 0 to SOVERSION 2 for `qd`,
  `qdmod`, and `qd_f_main`, avoiding accidental linkage against older local
  libQD3 installations.
- Changed the CMake release archive target to produce
  `libQD3-<version>.tar.gz` with a `libQD3-<version>/` archive prefix.
- Added a Branch-Free CMake QA gate on top of the default CMake matrix.

## Numerical Portability

- Strengthened CMake FMA auto-detection with a QD residual correctness probe.
  A compiler/runtime FMA candidate is selected only when it both runs and
  produces the residual expected by QD's split-product path.
- Added an x87 80-bit FPU mode probe and `fpu_fix_start_80bit`.
  DD/TD/QD tests continue to use the round-to-double x87 mode; EDD tests can
  switch to verified 80-bit extended mode where the platform supports it.
- Documented the known MinGW/Wine `edd_real` limitation: long-double
  trigonometric argument reduction may not satisfy the current EDD test
  tolerance under Wine, while DD/TD/QD tests pass in the same environment.

## Documentation

- Clarified unary worst-case eps-scaled errors.
- Updated README and NEWS for the 1.3.0 release.

## Scope Notes

- No new BLAS/LAPACK wrappers were added.
- No complex C API was added.
- No Fortran API expansion was added.
- The MinGW/Wine EDD trigonometry issue is documented but not fixed in this
  release.
