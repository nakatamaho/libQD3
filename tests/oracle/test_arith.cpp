#include "fn_registry.h"
#include "mpfr_oracle.h"
#include "qd_rng.h"
#include "tap.h"

#include <algorithm>
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

const int kSamples = 64;

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

void print_usage() {
  std::cout << "oracle_test_arith [-dd] [-td] [-qd] [-edd] [-all] [-v]"
            << " [--seed=N]\n";
}

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

template <class T>
std::vector<qd_oracle::Tap::Diagnostic>
base_diag(const char *program, const char *case_name, int iteration) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  std::ostringstream replay;
  replay << "tests/oracle/" << program << " -" << traits::name()
         << " --seed=" << qd_oracle::rng::active_seed();
  diag.push_back(qd_oracle::Tap::Diagnostic(
      "seed", std::to_string(qd_oracle::rng::active_seed())));
  diag.push_back(qd_oracle::Tap::Diagnostic("replay", replay.str()));
  diag.push_back(qd_oracle::Tap::Diagnostic("case", case_name));
  diag.push_back(qd_oracle::Tap::Diagnostic("iteration", int_text(iteration)));
  return diag;
}

template <class T>
void append_numeric_diag(std::vector<qd_oracle::Tap::Diagnostic> *diag,
                         const T &a, const T &b, const T &got,
                         mpfr_t ref, double relerr,
                         const qd_oracle::ErrorBound &bound) {
  diag->push_back(qd_oracle::Tap::Diagnostic("input_a_limbs",
                                             qd_oracle::limbs_hex(a)));
  diag->push_back(qd_oracle::Tap::Diagnostic("input_a_value",
                                             qd_oracle::value_to_mpfr_string(a)));
  diag->push_back(qd_oracle::Tap::Diagnostic("input_b_limbs",
                                             qd_oracle::limbs_hex(b)));
  diag->push_back(qd_oracle::Tap::Diagnostic("input_b_value",
                                             qd_oracle::value_to_mpfr_string(b)));
  diag->push_back(qd_oracle::Tap::Diagnostic("mpfr_reference",
                                             qd_oracle::mpfr_to_string(ref)));
  diag->push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                             qd_oracle::limbs_hex(got)));
  diag->push_back(qd_oracle::Tap::Diagnostic("got_value",
                                             qd_oracle::value_to_mpfr_string(got)));
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
T apply_arith(qd_oracle::ArithOp op, const T &a, const T &b) {
  switch (op) {
  case qd_oracle::arith_add:
    return a + b;
  case qd_oracle::arith_sub:
    return a - b;
  case qd_oracle::arith_mul:
    return a * b;
  case qd_oracle::arith_div:
    return a / b;
  case qd_oracle::arith_sqr:
    return sqr(a);
  }
  return T(0);
}

void apply_mpfr(qd_oracle::ArithOp op, mpfr_t ref, mpfr_t a, mpfr_t b) {
  switch (op) {
  case qd_oracle::arith_add:
    mpfr_add(ref, a, b, MPFR_RNDN);
    break;
  case qd_oracle::arith_sub:
    mpfr_sub(ref, a, b, MPFR_RNDN);
    break;
  case qd_oracle::arith_mul:
    mpfr_mul(ref, a, b, MPFR_RNDN);
    break;
  case qd_oracle::arith_div:
    mpfr_div(ref, a, b, MPFR_RNDN);
    break;
  case qd_oracle::arith_sqr:
    mpfr_sqr(ref, a, MPFR_RNDN);
    break;
  }
}

template <class T>
void make_inputs(qd_oracle::ArithOp op, T *a, T *b) {
  switch (op) {
  case qd_oracle::arith_add:
    *a = qd_oracle::rng::positive_type<T>(-220, 220);
    *b = qd_oracle::rng::positive_type<T>(-220, 220);
    if (qd_oracle::rng::coin()) {
      *a = -*a;
      *b = -*b;
    }
    break;
  case qd_oracle::arith_sub:
    *a = qd_oracle::rng::positive_type<T>(-220, 220);
    *b = qd_oracle::rng::positive_type<T>(-220, 220);
    if (*a < *b) {
      std::swap(*a, *b);
    }
    *b = *b / T(4);
    if (qd_oracle::rng::coin()) {
      *a = -*a;
      *b = -*b;
    }
    break;
  case qd_oracle::arith_mul:
    *a = qd_oracle::rng::uniform_type<T>(-140, 140);
    *b = qd_oracle::rng::uniform_type<T>(-140, 140);
    break;
  case qd_oracle::arith_div:
    *a = qd_oracle::rng::uniform_type<T>(-220, 220);
    *b = qd_oracle::rng::uniform_type<T>(-80, 80);
    break;
  case qd_oracle::arith_sqr:
    *a = qd_oracle::rng::uniform_type<T>(-140, 140);
    *b = T(0);
    break;
  }
}

template <class T>
bool run_exact_add_smoke(qd_oracle::Tap &tap, bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;

  const T a = T(1.25);
  const T b = T(2.5);
  const T got = a + b;
  const qd_oracle::ErrorBound bound = qd_oracle::exact_add_smoke_bound();

  mpfr_t a_mp;
  mpfr_t b_mp;
  mpfr_t ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, ref, (mpfr_ptr) 0);
  qd_oracle::to_mpfr(a_mp, a);
  qd_oracle::to_mpfr(b_mp, b);
  qd_oracle::require_exact(mpfr_add(ref, a_mp, b_mp, MPFR_RNDN),
                           "reference add");

  const double relerr = qd_oracle::relerr_in_eps(got, ref);
  const bool pass = relerr <= bound.eps_multiplier;

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = base_diag<T>("test_arith", "exact add smoke", 0);
    append_numeric_diag(&diag, a, b, got, ref, relerr, bound);
  }

  std::string name = std::string(traits::name()) + " exact add smoke";
  tap.ok(pass, name, diag, verbose);

  mpfr_clears(a_mp, b_mp, ref, (mpfr_ptr) 0);
  return pass;
}

template <class T>
bool run_arith_case(qd_oracle::Tap &tap,
                    const qd_oracle::ArithRegistryEntry &entry,
                    bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;

  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_a;
  T worst_b;
  T worst_got;
  mpfr_t worst_ref;
  mpfr_init2(worst_ref, qd_oracle::ref_prec<T>());

  mpfr_t a_mp;
  mpfr_t b_mp;
  mpfr_t ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, ref, (mpfr_ptr) 0);

  for (int i = 0; i < kSamples; ++i) {
    T a;
    T b;
    make_inputs(entry.op, &a, &b);
    T got = apply_arith(entry.op, a, b);

    qd_oracle::to_mpfr(a_mp, a);
    qd_oracle::to_mpfr(b_mp, b);
    apply_mpfr(entry.op, ref, a_mp, b_mp);

    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr;
      worst_iteration = i;
      worst_a = a;
      worst_b = b;
      worst_got = got;
      mpfr_set(worst_ref, ref, MPFR_RNDN);
    }
    if (!std::isfinite(relerr) || relerr > entry.bound.eps_multiplier) {
      pass = false;
    }
  }

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = base_diag<T>("test_arith", entry.name, worst_iteration);
    append_numeric_diag(&diag, worst_a, worst_b, worst_got, worst_ref, worst,
                        entry.bound);
  }

  std::string name = std::string(traits::name()) + " " + entry.name +
                     " random oracle";
  tap.ok(pass, name, diag, verbose);
  if (verbose || !pass) {
    std::cout << "# " << traits::name() << " " << entry.name
              << " worst_relerr_eps=" << worst << " samples=" << kSamples
              << "\n";
  }

  mpfr_clears(a_mp, b_mp, ref, (mpfr_ptr) 0);
  mpfr_clear(worst_ref);
  return pass;
}

template <class T>
bool run_compound_case(qd_oracle::Tap &tap) {
  typedef qd_oracle::TypeTraits<T> traits;
  bool pass = true;
  const char *failed_op = "";
  int failed_iteration = -1;
  T failed_a;
  T failed_b;
  T failed_got;
  T failed_expected;

  for (int i = 0; i < kSamples && pass; ++i) {
    T a;
    T b;
    make_inputs(qd_oracle::arith_add, &a, &b);
    T tmp = a;
    tmp += b;
    if (!qd_oracle::bit_equal(tmp, a + b)) {
      pass = false;
      failed_op = "+=";
      failed_iteration = i;
      failed_a = a;
      failed_b = b;
      failed_got = tmp;
      failed_expected = a + b;
      break;
    }

    make_inputs(qd_oracle::arith_sub, &a, &b);
    tmp = a;
    tmp -= b;
    if (!qd_oracle::bit_equal(tmp, a - b)) {
      pass = false;
      failed_op = "-=";
      failed_iteration = i;
      failed_a = a;
      failed_b = b;
      failed_got = tmp;
      failed_expected = a - b;
      break;
    }

    make_inputs(qd_oracle::arith_mul, &a, &b);
    tmp = a;
    tmp *= b;
    if (!qd_oracle::bit_equal(tmp, a * b)) {
      pass = false;
      failed_op = "*=";
      failed_iteration = i;
      failed_a = a;
      failed_b = b;
      failed_got = tmp;
      failed_expected = a * b;
      break;
    }

    make_inputs(qd_oracle::arith_div, &a, &b);
    tmp = a;
    tmp /= b;
    if (!qd_oracle::bit_equal(tmp, a / b)) {
      pass = false;
      failed_op = "/=";
      failed_iteration = i;
      failed_a = a;
      failed_b = b;
      failed_got = tmp;
      failed_expected = a / b;
      break;
    }
  }

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) {
    mpfr_t ref;
    mpfr_init2(ref, qd_oracle::ref_prec<T>());
    qd_oracle::to_mpfr(ref, failed_expected);
    diag = base_diag<T>("test_arith", "compound assignment", failed_iteration);
    diag.push_back(qd_oracle::Tap::Diagnostic("compound_op", failed_op));
    diag.push_back(qd_oracle::Tap::Diagnostic("input_a_limbs",
                                              qd_oracle::limbs_hex(failed_a)));
    diag.push_back(qd_oracle::Tap::Diagnostic("input_b_limbs",
                                              qd_oracle::limbs_hex(failed_b)));
    diag.push_back(qd_oracle::Tap::Diagnostic("mpfr_reference",
                                              qd_oracle::mpfr_to_string(ref)));
    diag.push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                              qd_oracle::limbs_hex(failed_got)));
    diag.push_back(qd_oracle::Tap::Diagnostic(
        "expected_limbs", qd_oracle::limbs_hex(failed_expected)));
    mpfr_clear(ref);
  }

  std::string name = std::string(traits::name()) + " compound assignment";
  tap.ok(pass, name, diag);
  return pass;
}

double random_double_operand() {
  const double mant = 1.0 +
      static_cast<double>(qd_oracle::rng::random_word52()) /
          static_cast<double>(1ULL << 52);
  double value = std::ldexp(mant, qd_oracle::rng::uniform_int(-40, 40));
  return qd_oracle::rng::coin() ? -value : value;
}

template <class T>
bool run_mixed_double_case(qd_oracle::Tap &tap) {
  typedef qd_oracle::TypeTraits<T> traits;
  const qd_oracle::ErrorBound add_bound = qd_oracle::add_sub_bound("mixed_add");
  const qd_oracle::ErrorBound mul_bound = qd_oracle::mul_bound("mixed_mul");

  bool pass = true;
  const char *failed_op = "";
  const qd_oracle::ErrorBound *failed_bound = &add_bound;
  int failed_iteration = -1;
  double failed_double = 0.0;
  T failed_a;
  T failed_got;
  mpfr_t failed_ref;
  mpfr_init2(failed_ref, qd_oracle::ref_prec<T>());
  double failed_relerr = 0.0;

  mpfr_t a_mp;
  mpfr_t d_mp;
  mpfr_t ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, d_mp, ref, (mpfr_ptr) 0);

  for (int i = 0; i < kSamples; ++i) {
    T a = qd_oracle::rng::uniform_type<T>(-120, 120);
    double d = random_double_operand();
    T got = a + d;
    qd_oracle::to_mpfr(a_mp, a);
    qd_oracle::require_exact(mpfr_set_d(d_mp, d, MPFR_RNDN),
                             "set mixed double");
    mpfr_add(ref, a_mp, d_mp, MPFR_RNDN);
    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (!std::isfinite(relerr) || relerr > add_bound.eps_multiplier) {
      pass = false;
      failed_op = "T+double";
      failed_bound = &add_bound;
      failed_iteration = i;
      failed_double = d;
      failed_a = a;
      failed_got = got;
      failed_relerr = relerr;
      mpfr_set(failed_ref, ref, MPFR_RNDN);
      break;
    }

    got = a * d;
    mpfr_mul(ref, a_mp, d_mp, MPFR_RNDN);
    relerr = qd_oracle::relerr_in_eps(got, ref);
    if (!std::isfinite(relerr) || relerr > mul_bound.eps_multiplier) {
      pass = false;
      failed_op = "T*double";
      failed_bound = &mul_bound;
      failed_iteration = i;
      failed_double = d;
      failed_a = a;
      failed_got = got;
      failed_relerr = relerr;
      mpfr_set(failed_ref, ref, MPFR_RNDN);
      break;
    }
  }

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) {
    T double_as_t(failed_double);
    diag = base_diag<T>("test_arith", "mixed double", failed_iteration);
    diag.push_back(qd_oracle::Tap::Diagnostic("mixed_op", failed_op));
    append_numeric_diag(&diag, failed_a, double_as_t, failed_got, failed_ref,
                        failed_relerr, *failed_bound);
  }

  std::string name = std::string(traits::name()) + " mixed double oracle";
  tap.ok(pass, name, diag);

  mpfr_clears(a_mp, d_mp, ref, (mpfr_ptr) 0);
  mpfr_clear(failed_ref);
  return pass;
}

template <class T>
bool run_comparison_case(qd_oracle::Tap &tap) {
  typedef qd_oracle::TypeTraits<T> traits;
  bool pass = true;
  int failed_iteration = -1;
  T failed_a;
  T failed_b;
  const char *failed_relation = "";

  mpfr_t a_mp;
  mpfr_t b_mp;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, (mpfr_ptr) 0);

  for (int i = 0; i < kSamples; ++i) {
    T a = qd_oracle::rng::uniform_type<T>(-160, 160);
    T b = (i % 8 == 0) ? a : qd_oracle::rng::uniform_type<T>(-160, 160);
    qd_oracle::to_mpfr(a_mp, a);
    qd_oracle::to_mpfr(b_mp, b);
    const int cmp = mpfr_cmp(a_mp, b_mp);

    if ((a == b) != (cmp == 0)) {
      pass = false;
      failed_relation = "==";
    } else if ((a != b) != (cmp != 0)) {
      pass = false;
      failed_relation = "!=";
    } else if ((a < b) != (cmp < 0)) {
      pass = false;
      failed_relation = "<";
    } else if ((a <= b) != (cmp <= 0)) {
      pass = false;
      failed_relation = "<=";
    } else if ((a > b) != (cmp > 0)) {
      pass = false;
      failed_relation = ">";
    } else if ((a >= b) != (cmp >= 0)) {
      pass = false;
      failed_relation = ">=";
    }

    if (!pass) {
      failed_iteration = i;
      failed_a = a;
      failed_b = b;
      break;
    }
  }

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) {
    diag = base_diag<T>("test_arith", "comparison", failed_iteration);
    diag.push_back(qd_oracle::Tap::Diagnostic("relation", failed_relation));
    diag.push_back(qd_oracle::Tap::Diagnostic("input_a_limbs",
                                              qd_oracle::limbs_hex(failed_a)));
    diag.push_back(qd_oracle::Tap::Diagnostic("input_b_limbs",
                                              qd_oracle::limbs_hex(failed_b)));
  }

  std::string name = std::string(traits::name()) + " comparison ordering";
  tap.ok(pass, name, diag);

  mpfr_clears(a_mp, b_mp, (mpfr_ptr) 0);
  return pass;
}

template <class T>
bool run_type(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  pass &= run_exact_add_smoke<T>(tap, verbose);

  std::size_t count = 0;
  const qd_oracle::ArithRegistryEntry *entries =
      qd_oracle::arithmetic_registry(&count);
  for (std::size_t i = 0; i < count; ++i) {
    pass &= run_arith_case<T>(tap, entries[i], verbose);
  }

  pass &= run_compound_case<T>(tap);
  pass &= run_mixed_double_case<T>(tap);
  pass &= run_comparison_case<T>(tap);
  return pass;
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

int planned_per_type() {
  std::size_t arith_count = 0;
  qd_oracle::arithmetic_registry(&arith_count);
  return 1 + static_cast<int>(arith_count) + 3;
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

  if (selected_count(*options) == 0) {
    select_all(options);
  }

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

  qd_oracle::Tap tap(selected_count(options) * planned_per_type());
  std::cout << "# seed: " << qd_oracle::rng::active_seed() << "\n";

  bool pass = true;
  unsigned int old_cw = 0;
  bool fpu_fixed = true;
  fpu_fix_start(&old_cw);

  if (options.test_dd) {
    pass &= run_type<dd_real>(tap, options.verbose);
  }
  if (options.test_td) {
    pass &= run_type<td_real>(tap, options.verbose);
  }
  if (options.test_qd) {
    pass &= run_type<qd_real>(tap, options.verbose);
  }

#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) {
    fpu_fix_end(&old_cw);
    fpu_fixed = false;
    pass &= run_type<edd_real>(tap, options.verbose);
  }
#endif

  if (fpu_fixed) {
    fpu_fix_end(&old_cw);
  }

  return pass ? tap.exit_status() : 1;
}
