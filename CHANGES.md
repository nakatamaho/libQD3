# CHANGES.md

## 1.6.0

libQD3 1.6.0 completes the C and Fortran interfaces for all six expansion
types, fixes numerical defects found by new epsilon-relative QA, and makes the
release gate cover Fortran.  See CHANGES.1.6.0.md for the full notes.


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
- `log` near 1 kept only absolute accuracy (the Newton step
  `x + a*exp(-x) - 1` cancels): relative errors reached 2e11 eps for
  `dd_real`, 1e9 eps for `qd_real` and 5e5 eps for the binary32 types.  For
  |a - 1| < 1/8 it now evaluates `log1p(a - 1)` (`a - 1` is exact there);
  all seven types are within 0.5 eps.  `log2`/`log10` inherit the fix.
- `ds/ts/qs` `fmod` rounded `b * n` before subtracting (up to 17754 eps); it
  now forms the products exactly and rounds once, and corrects `n` when
  `a / b` rounded across an integer.
- `cbrt(dd_real)` reached 33 eps for tiny arguments; it is evaluated in
  quad-double.
- `dd_real::_eps` and `qd_real::_eps` are now exactly 2^-104 and 2^-209.

Random numbers:

- All random functions used `std::rand()`.  `ds/ts/qs_real::rand()` (and the
  C and Fortran wrappers) scaled the first draw by 2^-24 and returned values
  in [0, 6e-8) on every platform.  `ddrand`/`qdrand` assumed 31-bit
  `rand()`; where `RAND_MAX` is 32767 (MinGW, MSVC) they returned values in
  [0, 1.5e-5) with zero gaps between the draws.
- They now share a xoshiro256** generator (`<qd/qd_random.h>`: `qd_srand`,
  `qd_rand_u64`; Fortran: `use qdrandom`, `call qd_random_seed(seed)`) and
  fill every limb exactly, uniformly in [0, 1), with identical sequences on
  every platform.  Added `tdrand()` and `eddrand()`.  Sequences differ from
  1.5.0.

MinGW / Windows:

- `fpu_fix_end` restored the raw x87 control word through `_control87` on
  MinGW, which reinterprets the bits and switched rounding to toward-zero
  (and the SSE rounding mode on x86-64); every `edd_real` operation after
  the first `fpu_fix_start`/`fpu_fix_end` pair was wrong.  It now uses
  `fldcw` like `fpu_fix_start`.
- The pi/16 constant used by `edd_real` argument reduction was computed in
  a DLL static initializer before the x87 precision is set to 64 bits; it is
  now an exact literal.  With both fixes the documented MinGW `edd_real`
  failures are gone.
- `ds/ts/qs` did not compile on MinGW (`<math.h>` defines a `_nan()` macro),
  and their exact products used `fmaf`, which the MinGW runtime does not
  round correctly (errors of ~3000 eps); products are now formed exactly in
  binary64 without FMA.
- `qd_f_main` cannot be a DLL (it calls the user's `f_main`); on Windows it
  is a static archive (`qd_f_main_shared` when linking the shared
  libraries).
- MPC is found without pkg-config when `mpc.pc` is missing (Debian, Ubuntu).
- New GitHub Actions CI: Linux (GCC, gfortran) and Windows (MSYS2 UCRT64
  MinGW-w64, gfortran), both with the MPFR/MPC oracles.

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
- The MPFR oracle drew Windows inputs with the top 20 mantissa bits zero
  (`mpfr_set_ui` with a 52-bit value); all platforms now draw the same
  inputs.  Every oracle program passes seeds 1-100.
- The MPC complex oracle also runs for `ds_complex`, `ts_complex` and
  `qs_complex` (seeded random inputs).
- Random functions are checked for range coverage, mean, limb fill and
  seed reproducibility in `regression_smoke`, `c_api_test` and `f_suite`.

Compatibility:

- `libqdmod` SOVERSION is now 3: unused Fortran module procedures were
  removed and `dble()` of ds/ts/qs values now returns `real*8`.  `libqd`
  keeps SOVERSION 2 (symbols were only added).

## 1.4.0

libQD3 1.4.0 hardens IEEE-style special values and overflow handling for all
supported expansion types.

Highlights:

- Added overflow-safe division rescaling for `dd_real`, `td_real`, `qd_real`,
  and `edd_real`, including scalar and mixed QD/DD division paths.
- Return canonical signed infinities when a quotient overflows before or after
  expansion reconstruction, avoiding residual `Inf - Inf` NaNs.
- Added consistent NaN, infinity, zero, and signed-zero handling for division
  and square root across DD/TD/QD/EDD.
- Fixed the EDD maximum constant on targets where the standard library does
  not specialize `numeric_limits<_Float64x>`.
- Added the `arithmetic_smoke` CTest regression covering special values,
  large-operand rescaling, direct quotient overflow, and mixed division.
- Verified the default CTest suite, the 16-configuration CMake matrix, and
  the binary-float matrix.

## 1.3.1

libQD3 1.3.1 is a focused x87/i386 portability fix release. It separates the
FPU mode used by EDD arithmetic from the round-to-double mode required by
DD/TD/QD oracle calculations and by the `qd_real` reduction path inside EDD
trigonometric functions.

Highlights:

- Fixed EDD trigonometric argument reduction to evaluate its internal
  `qd_real` work under QD round-to-double FPU mode before converting back to
  `edd_real`.
- Split EDD tests so `qd_real` and `dd_real` oracle operations do not run in
  EDD 80-bit FPU mode.
- Split `complex_test` so DD/TD/QD complex tests and EDD complex tests run
  under their respective FPU modes.
- Verified Debian i386 CTest without adding a global `-ffloat-store`
  workaround.

## 1.3.0

libQD3 1.3.0 is a maintenance release that upstreams the downstream patches
used by MPLAPACK. It improves C++ overload coverage for generated LAPACK-style
code, hardens CMake packaging, and records the remaining MinGW/Wine EDD
trigonometry limitation.

Highlights:

- Added integer mixed-mode arithmetic and comparison overloads for `dd_real`,
  `td_real`, `qd_real`, and `edd_real`.
- Added explicit integer conversion helpers for all real types.
- Added `std::complex<double>` interop operators for the libQD3 complex wrapper
  types.
- Added a CMake `uninstall` target and changed shared-library SOVERSION to 2.
- Strengthened FMA auto-detection with a residual correctness probe.
- Added an x87 80-bit FPU mode probe and EDD test wrappers.
- Documented the known MinGW/Wine `edd_real` long-double trigonometry issue.

## 1.2.0

### Extended-precision complex types

libQD3 now installs first-class complex headers for `dd_complex`, `td_complex`,
`qd_complex`, and `edd_complex` when `edd_real` is available. The implementation
uses one shared `qd3_complex<Real>` template body and thin type aliases, with
ADL-discovered overloads for arithmetic, elementary functions, `ldexp`,
`proj`, and component-wise `ceil`. The complex overload set now covers the
standard-compatible `tan`, `tanh`, inverse trigonometric, and inverse
hyperbolic functions without adding overloads to namespace `std`.

A default `complex_test` target covers the public headers and generic
`using std::sqrt; sqrt(z)` style calls. Optional MPC-backed complex oracle tests
are available with `-DQD3_ENABLE_MPC_TESTS=ON` and remain independent from the
existing MPFR-only real oracle suite.

The real types also gained matching practical C++ math coverage for `log2`,
`exp2`, `expm1`, `log1p`, `hypot`, `cbrt`, `trunc`, and `round`, with TD/EDD
API gaps closed where needed so generic `using std::...; f(x)` code behaves
consistently across DD, TD, QD, and EDD builds.

### Optional branch-free arithmetic

libQD3 now has optional branch-free (BF) addition and multiplication for
`dd_real`, `td_real`, and `qd_real`, based on Kouya's transcription of Zhang and
Aiken's algorithms in arXiv:2603.14926v2. BF arithmetic is off by default. Enable
it with `QD_BF`, or separately with `QD_BF_ADD` and `QD_BF_MUL`. CMake
options are `-DQD_ENABLE_BF=ON`, `-DQD_ENABLE_BF_ADD=ON`, and
`-DQD_ENABLE_BF_MUL=ON`.

| Type | Op | Default (unchanged) | With BF macro | Effect, scalar |
|------|----|---------------------|---------------|----------------|
| DD | add | `sloppy_add` (~11 flop) | `bf_add`, Algorithm 6 (~20 flop, error about `2u^2`) | ~2x add cost, ~1 bit more accuracy |
| DD | mul | current `operator*` | `bf_mul`, Algorithm 8 | No-op: current code already uses the same dataflow |
| TD | add | `sloppy_add` (~41 flop) | `bf_add`, Algorithm 11 (~57 flop) | Slower scalar; branch-free and SIMD-ready |
| TD | mul | renorm-based (~128 ops) | `bf_mul`, Algorithm 12 (~39 ops) | Faster, branch-free, exactly commutative |
| QD | add | `sloppy_add` (26/58) | `bf_add`, Algorithm 13 (37/66) | Slower scalar; branch-free and SIMD-ready |
| QD | mul | `accurate_mul` (renorm) | `bf_mul`, Algorithm 14 | Faster, branch-free, exactly commutative |

Addition and multiplication are split deliberately. `QD_BF_MUL` is the scalar
performance option for TD/QD multiplication, with DD multiplication unchanged.
`QD_BF_ADD` trades scalar addition speed for the BF dataflow and, for DD, a
tighter addition bound. Scalar builds should usually try `QD_BF_MUL`; SIMD or
cross-backend determinism builds can use full `QD_BF`.

BF paths use `TwoSum` and `TwoProd` dataflows directly. Their special-value
behavior can differ from the default path: signed zero may normalize to `+0.0`,
and `Inf - Inf` follows the EFT operations to `NaN`. Default builds keep the
existing special-value behavior.

Local scalar timing snapshot from `tests/qd_timer.cpp` on this container (GCC 15.2, `-O2 -ffp-contract=off`; one run, DD/QD only because this timer has no TD section):

| Build | DD add | DD mul | QD add | QD mul |
|-------|--------|--------|--------|--------|
| default | 90.44 mop/s | 109.18 mop/s | 6.38 mop/s | 5.57 mop/s |
| `QD_BF_MUL` | 82.81 mop/s | 107.35 mop/s | 6.22 mop/s | 7.08 mop/s |
| `QD_BF` | 37.58 mop/s | 108.71 mop/s | 8.31 mop/s | 7.10 mop/s |

The DD multiply rate is effectively unchanged, DD BF addition is slower, and QD BF multiplication is faster in this scalar run. Treat these as local smoke measurements rather than release-grade benchmark data.
