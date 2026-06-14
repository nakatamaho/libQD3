* update documentation
* oracle QA follow-up: implement optional JUnit output (`--junit=FILE` /
  `QD3_TEST_JUNIT`) or keep public docs strictly TAP-only until it exists
* oracle QA follow-up: wire Automake TAP-driver reporting for oracle subtests
  if per-subtest Automake summaries are required
* oracle QA follow-up: raise filtered lcov coverage from 69.7% line /
  75.2% function coverage toward the documented 90% function-body target
* oracle QA follow-up: add direct MPFR rows for `abs`/`fabs`, `asinh`,
  `acosh`, `atanh`, `sincos`, and `sincosh`
* oracle QA follow-up: extend rounding coverage for `nint`, `floor`, `ceil`,
  `aint`, and `quick_nint` with MPFR tie and large-value grids
* oracle QA follow-up: expand special-value coverage for signed zero,
  arithmetic NaN/Inf propagation, `_min_normalized`, `_max`, `_safe_max`,
  subnormal-like limb patterns, `log(0)`, `pow(0,0)`, and asin/acos
  out-of-domain inputs
* oracle QA follow-up: enumerate mixed dd/td/qd C++ overloads and mixed C API
  shims instead of relying on representative smoke coverage
* add integer format support
* support complex types.
* partial template specialization for complex divide.
* support x86 double-extended format (see Hladky's work)
* perhaps rewrite core code in C preprocessor, with C/C++ wrappers.
* complete numeric_limits
* sane handling of overflow / underflow / NaNs.
* handle more general streams, e.g., wide character streams
* use automake within the docs/ directory.
