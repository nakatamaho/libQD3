* update documentation
* replace td trig qd-assisted argument reduction with a native td reduction
  that avoids the qd round-trip while preserving sub-eps reduced-argument accuracy;
  simple staged subtraction of the fourth pi limb was tested and still left
  conditioned tan worst cases around 5-10 eps, so this likely needs true
  four-limb residual accumulation
* add integer format support
* support complex types.
* partial template specialization for complex divide.
* support x86 double-extended format (see Hladky's work)
* perhaps rewrite core code in C preprocessor, with C/C++ wrappers.
* complete numeric_limits
* sane handling of overflow / underflow / NaNs.
* handle more general streams, e.g., wide character streams
* use automake within the docs/ directory.
