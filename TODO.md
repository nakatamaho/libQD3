* update documentation
* oracle QA follow-up: implement optional JUnit output (`--junit=FILE` /
  `QD3_TEST_JUNIT`) or keep public docs strictly TAP-only until it exists
* oracle QA follow-up: wire Automake TAP-driver reporting for oracle subtests
  if per-subtest Automake summaries are required
* oracle QA follow-up: raise filtered lcov coverage from 69.7% line /
  75.2% function coverage toward the documented 90% function-body target
* oracle QA follow-up: audit the oracle registry and special-value matrix
  against every public unary, binary, rounding, I/O, and C API entry before release
* add integer format support
* support complex types.
* partial template specialization for complex divide.
* support x86 double-extended format (see Hladky's work)
* perhaps rewrite core code in C preprocessor, with C/C++ wrappers.
* complete numeric_limits
* sane handling of overflow / underflow / NaNs.
* handle more general streams, e.g., wide character streams
* use automake within the docs/ directory.
