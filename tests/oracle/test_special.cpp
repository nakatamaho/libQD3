#include "fn_registry.h"
#include "mpfr_oracle.h"
#include "tap.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <qd/fpu.h>

namespace {

struct Options {
  bool test_dd;
  bool test_td;
  bool test_qd;
  bool test_edd;
  bool verbose;
  Options()
      : test_dd(false), test_td(false), test_qd(false), test_edd(false),
        verbose(false) {}
};

std::string double_text(double value) {
  std::ostringstream os;
  os << value;
  return os.str();
}

void print_usage() {
  std::cout << "oracle_test_special [-dd] [-td] [-qd] [-edd] [-all] [-v]\n";
}

template <class T> struct SuppressErrors;

template <> struct SuppressErrors<dd_real> {
  bool old;
  SuppressErrors() : old(dd_suppress_error_messages) {
    dd_suppress_error_messages = true;
  }
  ~SuppressErrors() { dd_suppress_error_messages = old; }
};

template <> struct SuppressErrors<td_real> {
  bool old;
  SuppressErrors() : old(td_suppress_error_messages) {
    td_suppress_error_messages = true;
  }
  ~SuppressErrors() { td_suppress_error_messages = old; }
};

template <> struct SuppressErrors<qd_real> {
  bool old;
  SuppressErrors() : old(qd_suppress_error_messages) {
    qd_suppress_error_messages = true;
  }
  ~SuppressErrors() { qd_suppress_error_messages = old; }
};

#ifdef QD_HAVE_EDD_REAL
template <> struct SuppressErrors<edd_real> {
  bool old;
  bool old_qd;
  SuppressErrors()
      : old(edd_suppress_error_messages), old_qd(qd_suppress_error_messages) {
    edd_suppress_error_messages = true;
    qd_suppress_error_messages = true;
  }
  ~SuppressErrors() {
    edd_suppress_error_messages = old;
    qd_suppress_error_messages = old_qd;
  }
};
#endif

template <class T>
std::vector<qd_oracle::Tap::Diagnostic> base_diag(const char *case_name) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  std::ostringstream replay;
  replay << "tests/oracle/test_special -" << traits::name();
  diag.push_back(qd_oracle::Tap::Diagnostic("replay", replay.str()));
  diag.push_back(qd_oracle::Tap::Diagnostic("case", case_name));
  return diag;
}

template <class T>
bool report(qd_oracle::Tap &tap, const char *case_name, bool pass,
            const std::vector<qd_oracle::Tap::Diagnostic> &diag,
            bool verbose = false) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::string name = std::string(traits::name()) + " " + case_name;
  tap.ok(pass, name, diag, verbose);
  return pass;
}

template <class T>
bool check_close(qd_oracle::Tap &tap, const char *case_name, const T &got,
                 mpfr_t ref, double allowed_eps, bool verbose) {
  bool pass = true;
  double relerr = qd_oracle::relerr_in_eps(got, ref);
  if (!std::isfinite(relerr) || relerr > allowed_eps) pass = false;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = base_diag<T>(case_name);
    diag.push_back(qd_oracle::Tap::Diagnostic("mpfr_reference",
                                              qd_oracle::mpfr_to_string(ref)));
    diag.push_back(qd_oracle::Tap::Diagnostic("got_value",
                                              qd_oracle::value_to_mpfr_string(got)));
    diag.push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                              qd_oracle::limbs_hex(got)));
    diag.push_back(qd_oracle::Tap::Diagnostic("abs_error_mpfr",
                                              qd_oracle::abs_error_to_string(got, ref)));
    diag.push_back(qd_oracle::Tap::Diagnostic("relerr_eps",
                                              double_text(relerr)));
    diag.push_back(qd_oracle::Tap::Diagnostic("ulp_error_estimate",
                                              double_text(qd_oracle::ulp_error_estimate(got, ref))));
    diag.push_back(qd_oracle::Tap::Diagnostic("allowed_eps_multiplier",
                                              double_text(allowed_eps)));
  }
  return report<T>(tap, case_name, pass, diag, verbose);
}

template <class T>
bool run_nan_inf(qd_oracle::Tap &tap) {
  bool pass = true;
  pass &= isnan(T::_nan);
  pass &= isinf(T::_inf);
  pass &= isinf(abs(-T::_inf));
  pass &= !isnan(T(0));
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) diag = base_diag<T>("nan/inf constants");
  return report<T>(tap, "nan/inf constants", pass, diag);
}

template <class T>
bool run_domain_nan(qd_oracle::Tap &tap) {
  SuppressErrors<T> suppress;
  bool pass = true;
  pass &= isnan(sqrt(T(-1)));
  pass &= isnan(log(T(-1)));
  pass &= isnan(atan2(T(0), T(0)));
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) diag = base_diag<T>("domain nan contract");
  return report<T>(tap, "domain nan contract", pass, diag);
}

template <class T>
bool run_domain_edges(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  mpfr_t a, b, ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a, b, ref, (mpfr_ptr) 0);

  mpfr_set_ui(a, 1, MPFR_RNDN);
  mpfr_asin(ref, a, MPFR_RNDN);
  pass &= check_close<T>(tap, "asin(+1) edge", asin(T(1)), ref, 8.0, verbose);

  mpfr_set_si(a, -1, MPFR_RNDN);
  mpfr_acos(ref, a, MPFR_RNDN);
  pass &= check_close<T>(tap, "acos(-1) edge", acos(T(-1)), ref, 8.0, verbose);

  mpfr_set_si(a, 1, MPFR_RNDN);
  mpfr_set_si(b, -1, MPFR_RNDN);
  mpfr_atan2(ref, a, b, MPFR_RNDN);
  pass &= check_close<T>(tap, "atan2 quadrant II", atan2(T(1), T(-1)), ref,
                         16.0, verbose);

  mpfr_set_si(a, -8, MPFR_RNDN);
  mpfr_rootn_ui(ref, a, 3, MPFR_RNDN);
  pass &= check_close<T>(tap, "nroot(-8,3)", nroot(T(-8), 3), ref, 16.0, verbose);

  mpfr_clears(a, b, ref, (mpfr_ptr) 0);
  return pass;
}

template <class T>
bool run_rounding(qd_oracle::Tap &tap) {
  bool pass = true;
  pass &= (nint(T(1.5)) == T(2));
  pass &= (nint(T(-1.5)) == T(-1));
  pass &= (nint(T(2.0)) == T(2));
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) diag = base_diag<T>("nint half cases");
  return report<T>(tap, "nint half cases", pass, diag);
}

template <class T>
bool run_ldexp_power(qd_oracle::Tap &tap, bool verbose) {
  mpfr_t ref;
  mpfr_init2(ref, qd_oracle::ref_prec<T>());
  mpfr_set_ui(ref, 1, MPFR_RNDN);
  mpfr_mul_2si(ref, ref, 10, MPFR_RNDN);
  bool pass = check_close<T>(tap, "ldexp power-of-two", ldexp(T(1), 10), ref,
                             0.0, verbose);
  mpfr_clear(ref);
  return pass;
}

bool run_floor_ceil_aint(qd_oracle::Tap &tap, dd_real tag) {
  (void) tag;
  bool pass = true;
  pass &= (floor(dd_real(1.75)) == dd_real(1));
  pass &= (ceil(dd_real(-1.75)) == dd_real(-1));
  pass &= (aint(dd_real(-1.75)) == dd_real(-1));
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) diag = base_diag<dd_real>("floor/ceil/aint cases");
  return report<dd_real>(tap, "floor/ceil/aint cases", pass, diag);
}

bool run_floor_ceil_aint(qd_oracle::Tap &tap, qd_real tag) {
  (void) tag;
  bool pass = true;
  pass &= (floor(qd_real(1.75)) == qd_real(1));
  pass &= (ceil(qd_real(-1.75)) == qd_real(-1));
  pass &= (aint(qd_real(-1.75)) == qd_real(-1));
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) diag = base_diag<qd_real>("floor/ceil/aint cases");
  return report<qd_real>(tap, "floor/ceil/aint cases", pass, diag);
}

template <class T>
bool run_common_type(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  pass &= run_nan_inf<T>(tap);
  pass &= run_domain_nan<T>(tap);
  pass &= run_domain_edges<T>(tap, verbose);
  pass &= run_rounding<T>(tap);
  pass &= run_ldexp_power<T>(tap, verbose);
  return pass;
}

template <class T> bool run_type(qd_oracle::Tap &tap, bool verbose);

template <> bool run_type<dd_real>(qd_oracle::Tap &tap, bool verbose) {
  bool pass = run_common_type<dd_real>(tap, verbose);
  pass &= run_floor_ceil_aint(tap, dd_real(0));
  return pass;
}

template <> bool run_type<td_real>(qd_oracle::Tap &tap, bool verbose) {
  return run_common_type<td_real>(tap, verbose);
}

template <> bool run_type<qd_real>(qd_oracle::Tap &tap, bool verbose) {
  bool pass = run_common_type<qd_real>(tap, verbose);
  pass &= run_floor_ceil_aint(tap, qd_real(0));
  return pass;
}

#ifdef QD_HAVE_EDD_REAL
template <> bool run_type<edd_real>(qd_oracle::Tap &tap, bool verbose) {
  return run_common_type<edd_real>(tap, verbose);
}
#endif

int selected_count(const Options &options) {
  int count = 0;
  if (options.test_dd) count += 9;
  if (options.test_td) count += 8;
  if (options.test_qd) count += 9;
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) count += 8;
#endif
  return count;
}

void select_all(Options *options) {
  options->test_dd = true;
  options->test_td = true;
  options->test_qd = true;
#ifdef QD_HAVE_EDD_REAL
  options->test_edd = true;
#endif
}

bool parse_args(int argc, char **argv, Options *options) {
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-h") == 0 ||
        std::strcmp(argv[i], "-help") == 0) {
      print_usage();
      std::exit(0);
    } else if (std::strcmp(argv[i], "-dd") == 0) {
      options->test_dd = true;
    } else if (std::strcmp(argv[i], "-td") == 0) {
      options->test_td = true;
    } else if (std::strcmp(argv[i], "-qd") == 0) {
      options->test_qd = true;
    } else if (std::strcmp(argv[i], "-edd") == 0) {
#ifdef QD_HAVE_EDD_REAL
      options->test_edd = true;
#else
      std::cerr << "edd_real is not enabled in this build\n";
      return false;
#endif
    } else if (std::strcmp(argv[i], "-all") == 0) {
      select_all(options);
    } else if (std::strcmp(argv[i], "-v") == 0 ||
               std::strcmp(argv[i], "-verbose") == 0) {
      options->verbose = true;
    } else {
      std::cerr << "Unknown flag: " << argv[i] << "\n";
      return false;
    }
  }
  if (selected_count(*options) == 0) select_all(options);
  return true;
}

} // namespace

int main(int argc, char **argv) {
  Options options;
  if (!parse_args(argc, argv, &options)) {
    print_usage();
    return 2;
  }
  qd_oracle::Tap tap(selected_count(options));

  bool pass = true;
  unsigned int old_cw = 0;
  bool fpu_fixed = true;
  fpu_fix_start(&old_cw);
  if (options.test_dd) pass &= run_type<dd_real>(tap, options.verbose);
  if (options.test_td) pass &= run_type<td_real>(tap, options.verbose);
  if (options.test_qd) pass &= run_type<qd_real>(tap, options.verbose);
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) {
    fpu_fix_end(&old_cw);
    fpu_fixed = false;
    pass &= run_type<edd_real>(tap, options.verbose);
  }
#endif
  if (fpu_fixed) fpu_fix_end(&old_cw);
  return pass ? tap.exit_status() : 1;
}
