# Changes for libQD3 1.3.1

This file summarizes the changes after tag `1.3.0` through the 1.3.1 release.
This is a focused portability fix release for x87/i386 EDD testing and
trigonometric argument reduction.

## Numerical Portability

- Kept EDD arithmetic in verified 80-bit x87 mode while evaluating DD/TD/QD
  oracle expressions only inside scoped round-to-double FPU regions.
- Fixed EDD trigonometric argument reduction so its internal `qd_real`
  reduction work runs under the QD round-to-double FPU mode before converting
  the residual back to `edd_real`.
- Avoided using `-ffloat-store` as a global workaround. The fix is local to
  the places that need QD-mode evaluation and leaves EDD arithmetic in
  binary80 mode.

## Tests

- Split `complex_test` execution so `dd_complex`, `td_complex`, and
  `qd_complex` run under QD FPU mode, while `edd_complex` runs under 80-bit
  EDD FPU mode.
- Split `qd_test` EDD oracle checks so `qd_real` and `dd_real` construction,
  conversion, comparison, and transcendental oracle evaluation do not run while
  the process is in EDD 80-bit mode.
- Verified the Debian i386 CMake CTest suite with the MPLAPACK tier1 i386
  Docker image: 9/9 tests passed.
