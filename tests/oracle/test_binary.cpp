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
  std::cout << "oracle_test_binary [-dd] [-td] [-qd] [-edd] [-all] [-v]"
            << " [--seed=N]\n";
}

template <class T>
std::vector<qd_oracle::Tap::Diagnostic>
base_diag(const char *case_name, int iteration) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  std::ostringstream replay;
  replay << "tests/oracle/test_binary -" << traits::name()
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
                 const T &a, const T &b, const T &got, mpfr_t ref,
                 double relerr, const qd_oracle::ErrorBound &bound) {
  diag->push_back(qd_oracle::Tap::Diagnostic("input_a_limbs",
                                             qd_oracle::limbs_hex(a)));
  diag->push_back(qd_oracle::Tap::Diagnostic("input_b_limbs",
                                             qd_oracle::limbs_hex(b)));
  diag->push_back(qd_oracle::Tap::Diagnostic("mpfr_reference",
                                             qd_oracle::mpfr_to_string(ref)));
  diag->push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                             qd_oracle::limbs_hex(got)));
  diag->push_back(qd_oracle::Tap::Diagnostic("relerr_eps",
                                             double_text(relerr)));
  diag->push_back(qd_oracle::Tap::Diagnostic(
      "allowed_eps_multiplier", double_text(bound.eps_multiplier)));
  diag->push_back(qd_oracle::Tap::Diagnostic("bound_justification",
                                             bound.justification));
}

template <class T>
bool report_case(qd_oracle::Tap &tap, const char *case_name, bool pass,
                 double worst, int worst_iteration, const T &worst_a,
                 const T &worst_b, const T &worst_got, mpfr_t worst_ref,
                 const qd_oracle::ErrorBound &bound, bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) {
    diag = base_diag<T>(case_name, worst_iteration);
    append_diag(&diag, worst_a, worst_b, worst_got, worst_ref, worst, bound);
  }
  std::string name = std::string(traits::name()) + " " + case_name;
  tap.ok(pass, name, diag);
  if (verbose || pass) {
    std::cout << "# " << traits::name() << " " << case_name
              << " worst_relerr_eps=" << worst
              << " allowed_eps=" << bound.eps_multiplier
              << " samples=" << kSamples << "\n";
  }
  return pass;
}

template <class T>
bool run_pow_int(qd_oracle::Tap &tap, bool verbose) {
  qd_oracle::ErrorBound bound = qd_oracle::binary_algebraic_bound("pow_int");
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_a, worst_b, worst_got;
  mpfr_t a_mp, ref, worst_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, ref, worst_ref, (mpfr_ptr) 0);
  for (int i = 0; i < kSamples; ++i) {
    T a = qd_oracle::rng::uniform_type<T>(-3, 3);
    int n = qd_oracle::rng::uniform_int(-5, 5);
    if (n == 0) n = 3;
    T got = pow(a, n);
    qd_oracle::to_mpfr(a_mp, a);
    mpfr_pow_si(ref, a_mp, n, MPFR_RNDN);
    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr; worst_iteration = i; worst_a = a; worst_b = T(n);
      worst_got = got; mpfr_set(worst_ref, ref, MPFR_RNDN);
    }
    if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) pass = false;
  }
  bool result = report_case(tap, "pow_int oracle", pass, worst, worst_iteration,
                            worst_a, worst_b, worst_got, worst_ref, bound,
                            verbose);
  mpfr_clears(a_mp, ref, worst_ref, (mpfr_ptr) 0);
  return result;
}

template <class T>
bool run_nroot3(qd_oracle::Tap &tap, bool verbose) {
  qd_oracle::ErrorBound bound = qd_oracle::binary_algebraic_bound("nroot3");
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_a, worst_b, worst_got;
  mpfr_t a_mp, ref, worst_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, ref, worst_ref, (mpfr_ptr) 0);
  for (int i = 0; i < kSamples; ++i) {
    T a = qd_oracle::rng::uniform_type<T>(-18, 18);
    T got = nroot(a, 3);
    qd_oracle::to_mpfr(a_mp, a);
    mpfr_rootn_ui(ref, a_mp, 3, MPFR_RNDN);
    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr; worst_iteration = i; worst_a = a; worst_b = T(3);
      worst_got = got; mpfr_set(worst_ref, ref, MPFR_RNDN);
    }
    if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) pass = false;
  }
  bool result = report_case(tap, "nroot3 oracle", pass, worst, worst_iteration,
                            worst_a, worst_b, worst_got, worst_ref, bound,
                            verbose);
  mpfr_clears(a_mp, ref, worst_ref, (mpfr_ptr) 0);
  return result;
}

template <class T>
bool run_atan2(qd_oracle::Tap &tap, bool verbose) {
  qd_oracle::ErrorBound bound =
      qd_oracle::binary_transcendental_bound("atan2");
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_y, worst_x, worst_got;
  mpfr_t y_mp, x_mp, ref, worst_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), y_mp, x_mp, ref, worst_ref,
              (mpfr_ptr) 0);
  for (int i = 0; i < kSamples; ++i) {
    T y = qd_oracle::rng::uniform_type<T>(-8, 8);
    T x = qd_oracle::rng::uniform_type<T>(-8, 8);
    T got = atan2(y, x);
    qd_oracle::to_mpfr(y_mp, y);
    qd_oracle::to_mpfr(x_mp, x);
    mpfr_atan2(ref, y_mp, x_mp, MPFR_RNDN);
    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr; worst_iteration = i; worst_y = y; worst_x = x;
      worst_got = got; mpfr_set(worst_ref, ref, MPFR_RNDN);
    }
    if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) pass = false;
  }
  bool result = report_case(tap, "atan2 oracle", pass, worst, worst_iteration,
                            worst_y, worst_x, worst_got, worst_ref, bound,
                            verbose);
  mpfr_clears(y_mp, x_mp, ref, worst_ref, (mpfr_ptr) 0);
  return result;
}

template <class T>
bool run_ldexp(qd_oracle::Tap &tap, bool verbose) {
  qd_oracle::ErrorBound bound = qd_oracle::exact_ldexp_bound();
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_a, worst_b, worst_got;
  mpfr_t a_mp, ref, worst_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, ref, worst_ref, (mpfr_ptr) 0);
  for (int i = 0; i < kSamples; ++i) {
    T a = qd_oracle::rng::uniform_type<T>(-40, 40);
    int n = qd_oracle::rng::uniform_int(-20, 20);
    T got = ldexp(a, n);
    qd_oracle::to_mpfr(a_mp, a);
    mpfr_mul_2si(ref, a_mp, n, MPFR_RNDN);
    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr; worst_iteration = i; worst_a = a; worst_b = T(n);
      worst_got = got; mpfr_set(worst_ref, ref, MPFR_RNDN);
    }
    if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) pass = false;
  }
  bool result = report_case(tap, "ldexp exact oracle", pass, worst,
                            worst_iteration, worst_a, worst_b, worst_got,
                            worst_ref, bound, verbose);
  mpfr_clears(a_mp, ref, worst_ref, (mpfr_ptr) 0);
  return result;
}

template <class T>
bool run_pow_real_impl(qd_oracle::Tap &tap, bool verbose) {
  qd_oracle::ErrorBound bound =
      qd_oracle::binary_transcendental_bound("pow_real");
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_a, worst_b, worst_got;
  mpfr_t a_mp, b_mp, ref, worst_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, ref, worst_ref,
              (mpfr_ptr) 0);
  for (int i = 0; i < kSamples; ++i) {
    T a = qd_oracle::rng::positive_type<T>(-4, 4);
    T b = qd_oracle::rng::uniform_type<T>(-2, 2);
    T got = pow(a, b);
    qd_oracle::to_mpfr(a_mp, a);
    qd_oracle::to_mpfr(b_mp, b);
    mpfr_pow(ref, a_mp, b_mp, MPFR_RNDN);
    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr; worst_iteration = i; worst_a = a; worst_b = b;
      worst_got = got; mpfr_set(worst_ref, ref, MPFR_RNDN);
    }
    if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) pass = false;
  }
  bool result = report_case(tap, "pow_real oracle", pass, worst,
                            worst_iteration, worst_a, worst_b, worst_got,
                            worst_ref, bound, verbose);
  mpfr_clears(a_mp, b_mp, ref, worst_ref, (mpfr_ptr) 0);
  return result;
}

template <class T>
bool run_fmod_impl(qd_oracle::Tap &tap, bool verbose) {
  qd_oracle::ErrorBound bound = qd_oracle::binary_algebraic_bound("fmod");
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_a, worst_b, worst_got;
  mpfr_t a_mp, b_mp, ref, worst_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, ref, worst_ref,
              (mpfr_ptr) 0);
  for (int i = 0; i < kSamples; ++i) {
    double bd = std::ldexp(1.0 + static_cast<double>(i % 7) / 8.0,
                           qd_oracle::rng::uniform_int(-5, 5));
    double rd = bd * (0.25 + static_cast<double>((i % 5) + 1) / 16.0);
    int q = qd_oracle::rng::uniform_int(-12, 12);
    T b(bd);
    T a = T(q) * b + T(rd);
    if (qd_oracle::rng::coin()) a = -a;
    T got = fmod(a, b);
    qd_oracle::to_mpfr(a_mp, a);
    qd_oracle::to_mpfr(b_mp, b);
    mpfr_fmod(ref, a_mp, b_mp, MPFR_RNDN);
    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr; worst_iteration = i; worst_a = a; worst_b = b;
      worst_got = got; mpfr_set(worst_ref, ref, MPFR_RNDN);
    }
    if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) pass = false;
  }
  bool result = report_case(tap, "fmod oracle", pass, worst, worst_iteration,
                            worst_a, worst_b, worst_got, worst_ref, bound,
                            verbose);
  mpfr_clears(a_mp, b_mp, ref, worst_ref, (mpfr_ptr) 0);
  return result;
}

template <class T>
bool run_drem_impl(qd_oracle::Tap &tap, bool verbose) {
  qd_oracle::ErrorBound bound = qd_oracle::binary_algebraic_bound("drem");
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_a, worst_b, worst_got;
  mpfr_t a_mp, b_mp, ref, worst_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, ref, worst_ref,
              (mpfr_ptr) 0);
  for (int i = 0; i < kSamples; ++i) {
    double bd = std::ldexp(1.0 + static_cast<double>(i % 9) / 16.0,
                           qd_oracle::rng::uniform_int(-5, 5));
    double rd = bd * (0.125 + static_cast<double>((i % 3) + 1) / 16.0);
    int q = qd_oracle::rng::uniform_int(-12, 12);
    T b(bd);
    T a = T(q) * b + T(rd);
    if (qd_oracle::rng::coin()) a = -a;
    T got = drem(a, b);
    qd_oracle::to_mpfr(a_mp, a);
    qd_oracle::to_mpfr(b_mp, b);
    mpfr_remainder(ref, a_mp, b_mp, MPFR_RNDN);
    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr; worst_iteration = i; worst_a = a; worst_b = b;
      worst_got = got; mpfr_set(worst_ref, ref, MPFR_RNDN);
    }
    if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) pass = false;
  }
  bool result = report_case(tap, "drem oracle", pass, worst, worst_iteration,
                            worst_a, worst_b, worst_got, worst_ref, bound,
                            verbose);
  mpfr_clears(a_mp, b_mp, ref, worst_ref, (mpfr_ptr) 0);
  return result;
}

template <class T>
bool run_common_type(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  pass &= run_pow_int<T>(tap, verbose);
  pass &= run_nroot3<T>(tap, verbose);
  pass &= run_atan2<T>(tap, verbose);
  pass &= run_ldexp<T>(tap, verbose);
  return pass;
}

template <class T>
bool run_type(qd_oracle::Tap &tap, bool verbose);

template <>
bool run_type<dd_real>(qd_oracle::Tap &tap, bool verbose) {
  bool pass = run_common_type<dd_real>(tap, verbose);
  pass &= run_pow_real_impl<dd_real>(tap, verbose);
  pass &= run_fmod_impl<dd_real>(tap, verbose);
  pass &= run_drem_impl<dd_real>(tap, verbose);
  return pass;
}

template <>
bool run_type<td_real>(qd_oracle::Tap &tap, bool verbose) {
  bool pass = run_common_type<td_real>(tap, verbose);
  pass &= run_pow_real_impl<td_real>(tap, verbose);
  return pass;
}

template <>
bool run_type<qd_real>(qd_oracle::Tap &tap, bool verbose) {
  bool pass = run_common_type<qd_real>(tap, verbose);
  pass &= run_pow_real_impl<qd_real>(tap, verbose);
  pass &= run_fmod_impl<qd_real>(tap, verbose);
  pass &= run_drem_impl<qd_real>(tap, verbose);
  return pass;
}

#ifdef QD_HAVE_EDD_REAL
template <>
bool run_type<edd_real>(qd_oracle::Tap &tap, bool verbose) {
  return run_common_type<edd_real>(tap, verbose);
}
#endif

int type_plan_dd_qd() { return 7; }
int type_plan_td() { return 5; }
int type_plan_edd() { return 4; }

int selected_plan(const Options &options) {
  int count = 0;
  if (options.test_dd) count += type_plan_dd_qd();
  if (options.test_td) count += type_plan_td();
  if (options.test_qd) count += type_plan_dd_qd();
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) count += type_plan_edd();
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
