# CHANGES.md

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
