# Changes for libQD3 1.6.0

libQD3 1.6.0 completes the C and Fortran interfaces for all six expansion
types, fixes numerical defects found by new epsilon-relative QA, and makes the
release gate cover Fortran.

C++ core:

- Added `<qd/ds_complex.h>`, `<qd/ts_complex.h>` and `<qd/qs_complex.h>`
  (`ds_complex`, `ts_complex`, `qs_complex`), included by `<qd/complex.h>`.
- Decimal parsing: dd/td/qd/edd multiplied by an inexact 10^-k, so exactly
  representable literals such as `td_real("3.0")` or `dd_real("0.0625")` were
  off by a fraction of an ulp; they now divide by the exact power.  The
  binary32 reader (ds/ts/qs) truncated to a fixed digit count and scaled in
  target precision (up to 6.7 eps error); it now converts with quad-double
  precision and rounds once (< 0.1 eps).
- `ds/ts/qs_real` built from `double`/`long double` `-0.0` lost the sign.
- `dd_real::_eps` and `qd_real::_eps` are now exactly 2^-104 and 2^-209.

C API:

- New `c_ds.h`, `c_ts.h`, `c_qs.h` with the same function set as the
  dd/td/qd C API.
- Mixed-precision comparisons in `c_dd`/`c_td` rounded the wider operand to
  the narrower type (1 and 1 + 1e-40 compared equal); they now compare in
  the wider type.
- Filled gaps: same-type `c_dd_self*`; `c_td` rounding functions, `nroot`,
  `rand` and dd comparisons; `c_qd` dd comparisons; `c_edd` mixed double
  arithmetic, self operators, powers/roots, rounding, inverse trigonometric
  and hyperbolic functions, `sincos`/`sincosh`, `comp_d_edd`, `rand`.

Fortran:

- `nint` truncated for dd/qd/ds/qs; `int()` and integer assignment used only
  the leading limb; `sign(a, +0)` returned `-|a|`.
- ds/ts/qs: integers above 2^24 were rounded to one binary32 limb in all
  mixed operations; `dble()` returned the leading binary32 limb; there was no
  `real*8`/`complex*16` interoperability (now added).
- ts/qs <-> dd/qd conversions copied limbs between binary32 and binary64
  arrays (`ts = dd` kept ~1e-8 accuracy); they now convert through C++.
- dsmodule and qsmodule now provide the same dd/qd conversions as tsmodule
  (real and complex, constructors and assignment in both directions).
- Literal parsing now uses the C++ readers (the old Fortran parsers were
  inexact even for `0.5`); `d`/`D` exponents are accepted.
- The qs `epsilon()` parameter was 1.25 * 2^-94.
- Added integer `-`, complex mixed-mode `+ - / == /=`, complex elementary
  functions (`sqrt`, trigonometric, hyperbolic and their inverses,
  `complex**complex`, `complex**real`), `floor`/`ceiling` for all modules,
  `precision`/`range` for dd/ds, `hypot`, `modulo` and `dim`.

QA:

- `f_suite`: type-generic Fortran precision suite with MPFR/MPC references
  and epsilon-relative tolerances for all six modules (`f_test` tolerances
  for ds/ts/qs were 1e-5 and qd results were not checked).
- `c_api_test`: generated conformance test for all 599 C API functions.
- The release-gate scripts build with `QD_BUILD_FORTRAN=ON`; set
  `QD3_QA_FORTRAN=OFF` to opt out explicitly.
- `regression_smoke` checks exact and long-literal parsing for all types.
- The MPC complex oracle also runs for `ds_complex`, `ts_complex` and
  `qs_complex` (seeded random inputs).

Compatibility:

- `libqdmod` SOVERSION is now 3: unused Fortran module procedures were
  removed and `dble()` of ds/ts/qs values now returns `real*8`.  `libqd`
  keeps SOVERSION 2 (symbols were only added).
