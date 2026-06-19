# CHANGES.md

## Unreleased

### Optional branch-free arithmetic

libQD3 now has optional branch-free (BF) addition and multiplication for
`dd_real`, `td_real`, and `qd_real`, based on Kouya's transcription of Zhang and
Aiken's algorithms in arXiv:2603.14926v2. BF arithmetic is off by default. Enable
it with `QD_BF`, or separately with `QD_BF_ADD` and `QD_BF_MUL`; Autotools also
accepts `--enable-bf`, `--enable-bf-add`, and `--enable-bf-mul`.

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
