# TODO

Open items after 1.6.0.  Completed items were removed (complex types,
x86 double-extended support via `edd_real`, Automake integration, which no
longer applies to the CMake-only build, and the mixed C API coverage now in
`c_api_test`).

## Interfaces

* complex C API: no `c_*` entry points exist for `dd_complex` ... `qs_complex`.
* `edd_real` Fortran module (the C and C++ interfaces exist).
* Fortran model inquiry functions for the expansion types: `exponent`,
  `fraction`, `scale`, `set_exponent`, `spacing`, `nearest`, `rrspacing`.
* `std::numeric_limits` is incomplete: members not overridden are inherited
  from `double`/`float` (for example `max_digits10` is 17 for `dd_real`), and
  `infinity()`, `quiet_NaN()`, `lowest()` and `denorm_min()` return the limb
  type rather than the expansion type.
* integer formatting support in the I/O routines.
* wide-character and other general stream support.

## Numerics

* overflow / underflow / NaN handling beyond division and square root (which
  were hardened in 1.4.0).
* partial template specialization for complex division.

## QA

* oracle: add direct MPFR rows for `abs`/`fabs`, `asinh`, `acosh`, `atanh`,
  `sincos`, and `sincosh`.
* oracle: MPFR rounding coverage for `nint`, `floor`, `ceil`, `aint`, and
  `quick_nint` with tie and large-value grids (`test_rounding_corners` does not
  exercise them).
* oracle: expand special-value coverage for signed zero, arithmetic NaN/Inf
  propagation, `_min_normalized`, `_max`, `_safe_max`, subnormal-like limb
  patterns, `log(0)`, `pow(0,0)`, and asin/acos out-of-domain inputs.
* oracle: enumerate the mixed dd/td/qd C++ operator overloads (the mixed C API
  shims are covered by `c_api_test`).
* raise filtered lcov coverage toward the 90% function-body target (last
  recorded: 69.7% line / 75.2% function, before 1.6.0).
* optional JUnit output (`--junit=FILE` / `QD3_TEST_JUNIT`), or keep public
  docs strictly TAP-only.
* CI runs the default and oracle suites only; the 16-configuration and BF
  matrices of the release gate are still run by hand.
* platforms not verified for 1.6.0: i386/x87, macOS, MSVC, and Fortran
  compilers other than gfortran.

## Documentation

* document `ds_real`/`ts_real`/`qs_real`, the C API and the Fortran modules
  in `docs/` (only `qd.tex`, `td.tex` and `edd.tex` exist).
* build the LaTeX documentation from CMake (`docs/Makefile` is standalone).

## Ideas

* rewrite the core code with the C preprocessor and thin C/C++ wrappers.
