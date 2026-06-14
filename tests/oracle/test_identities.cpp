#include "fn_registry.h"
#include "mpfr_oracle.h"
#include "qd_rng.h"
#include "tap.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <qd/fpu.h>

namespace {

const int kSamples = 48;

struct Options {
  bool test_dd;
  bool test_td;
  bool test_qd;
  bool test_edd;
  bool verbose;
  bool has_seed;
  std::uint64_t seed;

  Options()
      : test_dd(false), test_td(false), test_qd(false), test_edd(false),
        verbose(false), has_seed(false), seed(0) {}
};

std::string double_text(double value) {
  std::ostringstream os;
  os << value;
  return os.str();
}

std::string int_text(int value) {
  std::ostringstream os;
  os << value;
  return os.str();
}

void print_usage() {
  std::cout << "oracle_test_identities [-dd] [-td] [-qd] [-edd] [-all] [-v]"
            << " [--seed=N]\n";
}

template <class T>
std::vector<qd_oracle::Tap::Diagnostic>
base_diag(const char *case_name, int iteration) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  std::ostringstream replay;
  replay << "tests/oracle/test_identities -" << traits::name()
         << " --seed=" << qd_oracle::rng::active_seed();
  diag.push_back(qd_oracle::Tap::Diagnostic(
      "seed", std::to_string(qd_oracle::rng::active_seed())));
  diag.push_back(qd_oracle::Tap::Diagnostic("replay", replay.str()));
  diag.push_back(qd_oracle::Tap::Diagnostic("case", case_name));
  diag.push_back(qd_oracle::Tap::Diagnostic("iteration", int_text(iteration)));
  return diag;
}

template <class T>
void append_diag(std::vector<qd_oracle::Tap::Diagnostic> *diag,
                 const T &input, const T &got, mpfr_t ref, double relerr,
                 const qd_oracle::ErrorBound &bound) {
  diag->push_back(qd_oracle::Tap::Diagnostic("input_limbs",
                                             qd_oracle::limbs_hex(input)));
  diag->push_back(qd_oracle::Tap::Diagnostic("input_value",
                                             qd_oracle::value_to_mpfr_string(input)));
  diag->push_back(qd_oracle::Tap::Diagnostic("mpfr_reference",
                                             qd_oracle::mpfr_to_string(ref)));
  diag->push_back(qd_oracle::Tap::Diagnostic("got_value",
                                             qd_oracle::value_to_mpfr_string(got)));
  diag->push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                             qd_oracle::limbs_hex(got)));
  diag->push_back(qd_oracle::Tap::Diagnostic("abs_error_mpfr",
                                             qd_oracle::abs_error_to_string(got, ref)));
  diag->push_back(qd_oracle::Tap::Diagnostic("relerr_eps",
                                             double_text(relerr)));
  diag->push_back(qd_oracle::Tap::Diagnostic("ulp_error_estimate",
                                             double_text(qd_oracle::ulp_error_estimate(got, ref))));
  diag->push_back(qd_oracle::Tap::Diagnostic(
      "allowed_eps_multiplier", double_text(bound.eps_multiplier)));
  diag->push_back(qd_oracle::Tap::Diagnostic("bound_justification",
                                             bound.justification));
}

template <class T>
bool report_case(qd_oracle::Tap &tap, const char *case_name, bool pass,
                 double worst, int worst_iteration, const T &worst_input,
                 const T &worst_got, mpfr_t worst_ref,
                 const qd_oracle::ErrorBound &bound, bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = base_diag<T>(case_name, worst_iteration);
    append_diag(&diag, worst_input, worst_got, worst_ref, worst, bound);
  }
  std::string name = std::string(traits::name()) + " " + case_name;
  tap.ok(pass, name, diag, verbose);
  if (verbose || !pass) {
    std::cout << "# " << traits::name() << " " << case_name
              << " worst_relerr_eps=" << worst
              << " allowed_eps=" << bound.eps_multiplier
              << " samples=" << kSamples << "\n";
  }
  return pass;
}

template <class T>
bool check_identity(qd_oracle::Tap &tap, const char *case_name,
                    const qd_oracle::ErrorBound &bound, bool verbose,
                    T (*make_input)(), T (*apply_identity)(const T &),
                    void (*make_reference)(mpfr_t, const T &)) {
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_input;
  T worst_got;
  mpfr_t ref;
  mpfr_t worst_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), ref, worst_ref, (mpfr_ptr) 0);

  for (int i = 0; i < kSamples; ++i) {
    T input = make_input();
    T got = apply_identity(input);
    make_reference(ref, input);
    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr;
      worst_iteration = i;
      worst_input = input;
      worst_got = got;
      mpfr_set(worst_ref, ref, MPFR_RNDN);
    }
    if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) {
      pass = false;
    }
  }

  bool result = report_case(tap, case_name, pass, worst, worst_iteration,
                            worst_input, worst_got, worst_ref, bound,
                            verbose);
  mpfr_clears(ref, worst_ref, (mpfr_ptr) 0);
  return result;
}

template <class T>
T input_positive_wide() {
  return qd_oracle::rng::positive_type<T>(-8, 8);
}

template <class T>
T input_positive_moderate() {
  return qd_oracle::rng::positive_type<T>(-6, 6);
}

template <class T>
T input_moderate() {
  return qd_oracle::rng::uniform_type<T>(-2, 2);
}

template <class T>
T input_principal() {
  return qd_oracle::rng::uniform_type<T>(-4, -1);
}

template <class T>
T input_trig_stable() {
  return qd_oracle::rng::uniform_type<T>(-3, -1);
}

template <class T>
T identity_exp_log(const T &x) {
  return exp(log(x));
}

template <class T>
T identity_log_exp(const T &x) {
  return log(exp(x));
}

template <class T>
T identity_sin2_cos2(const T &x) {
  T s = sin(x);
  T c = cos(x);
  return sqr(s) + sqr(c);
}

template <class T>
T identity_sqrt_sqr(const T &x) {
  T r = sqrt(x);
  return sqr(r);
}

template <class T>
T identity_sqr_mul(const T &x) {
  return sqr(x);
}

template <class T>
T identity_nroot3_cube(const T &x) {
  T r = nroot(x, 3);
  return sqr(r) * r;
}

template <class T>
T identity_tan_sincos(const T &x) {
  return tan(x);
}

template <class T>
T identity_asin_sin(const T &x) {
  return asin(sin(x));
}

template <class T>
T identity_atan_tan(const T &x) {
  return atan(tan(x));
}

template <class T>
void ref_input(mpfr_t ref, const T &x) {
  qd_oracle::to_mpfr(ref, x);
}

template <class T>
void ref_one(mpfr_t ref, const T &) {
  mpfr_set_ui(ref, 1, MPFR_RNDN);
}

template <class T>
void ref_square(mpfr_t ref, const T &x) {
  T rhs = x * x;
  qd_oracle::to_mpfr(ref, rhs);
}

template <class T>
void ref_tan_sincos(mpfr_t ref, const T &x) {
  T rhs = sin(x) / cos(x);
  qd_oracle::to_mpfr(ref, rhs);
}

template <class T>
bool run_type(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  pass &= check_identity<T>(
      tap, "exp(log(x)) identity", qd_oracle::identity_bound("exp_log"),
      verbose, input_positive_wide<T>, identity_exp_log<T>, ref_input<T>);
  pass &= check_identity<T>(
      tap, "log(exp(x)) identity", qd_oracle::identity_bound("log_exp"),
      verbose, input_moderate<T>, identity_log_exp<T>, ref_input<T>);
  pass &= check_identity<T>(
      tap, "sin^2+cos^2 identity", qd_oracle::identity_bound("sin2_cos2"),
      verbose, input_trig_stable<T>, identity_sin2_cos2<T>, ref_one<T>);
  pass &= check_identity<T>(
      tap, "sqrt(x)^2 identity", qd_oracle::identity_bound("sqrt_sqr"),
      verbose, input_positive_moderate<T>, identity_sqrt_sqr<T>,
      ref_input<T>);
  pass &= check_identity<T>(
      tap, "sqr(x) vs x*x identity", qd_oracle::identity_bound("sqr_mul"),
      verbose, input_moderate<T>, identity_sqr_mul<T>, ref_square<T>);
  pass &= check_identity<T>(
      tap, "nroot(x,3)^3 identity", qd_oracle::identity_bound("nroot3_cube"),
      verbose, input_moderate<T>, identity_nroot3_cube<T>, ref_input<T>);
  pass &= check_identity<T>(
      tap, "tan(x) vs sin(x)/cos(x) identity",
      qd_oracle::identity_bound("tan_sincos"), verbose, input_trig_stable<T>,
      identity_tan_sincos<T>, ref_tan_sincos<T>);
  pass &= check_identity<T>(
      tap, "asin(sin(x)) principal identity",
      qd_oracle::identity_bound("asin_sin"), verbose, input_principal<T>,
      identity_asin_sin<T>, ref_input<T>);
  pass &= check_identity<T>(
      tap, "atan(tan(x)) principal identity",
      qd_oracle::identity_bound("atan_tan"), verbose, input_principal<T>,
      identity_atan_tan<T>, ref_input<T>);
  return pass;
}

int type_plan() { return 9; }

int selected_plan(const Options &options) {
  int count = 0;
  if (options.test_dd) count += type_plan();
  if (options.test_td) count += type_plan();
  if (options.test_qd) count += type_plan();
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) count += type_plan();
#endif
  return count;
}

int selected_count(const Options &options) {
  int count = 0;
  if (options.test_dd) ++count;
  if (options.test_td) ++count;
  if (options.test_qd) ++count;
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) ++count;
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
      std::uint64_t seed = 0;
      if (qd_oracle::rng::parse_seed_arg(argv[i], &seed)) {
        options->has_seed = true;
        options->seed = seed;
      } else {
        std::cerr << "Unknown flag `" << argv[i] << "'.\n";
        return false;
      }
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
  qd_oracle::rng::configure(options.has_seed, options.seed);
  qd_oracle::Tap tap(selected_plan(options));
  std::cout << "# seed: " << qd_oracle::rng::active_seed() << "\n";

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
