# Changes for libQD3 1.4.0

This release hardens division and special-value behavior across the DD, TD,
QD, and EDD expansion types.

## Numerical robustness

- Division rescales very large numerators before expansion arithmetic and now
  checks both the initial scalar quotient and the final rescaled result for
  overflow.
- Direct and reconstructed overflow returns canonical signed infinities,
  preventing an intermediate `Inf - Inf` from becoming a NaN.
- Division consistently handles NaN, infinity, zero, and signed zero for
  scalar, expansion, and QD/DD mixed operands.
- Square root now handles NaN, positive and negative infinity, and signed zero
  explicitly for DD, TD, QD, and EDD.
- Fixed `edd_real::_max` on platforms where the standard library does not
  provide `numeric_limits<_Float64x>`.

## Tests

- Added `arithmetic_smoke`, a CTest regression covering special values, large
  division and square root, direct quotient overflow, rescaled overflow,
  scalar division, and QD/DD mixed division.
- Verified the default CTest suite with EDD enabled: 11/11 tests passed.
- Verified all 16 configurations in `qa/check_16_builds_cmake.sh`.
- Verified all 5 configurations in `qa/check_bf_matrix_cmake.sh`.
- Verified the MPLAPACK QD SVD regression case containing the previous
  `Rgesvj` `Inf`/`NaN` failure.
