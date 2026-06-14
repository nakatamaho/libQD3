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
const double kStableTrigThreshold = 0.25;

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

struct TrigContext {
  double condition;
  double sin_component_relerr;
  double cos_component_relerr;
  double tan_via_sin_cos_relerr;
  std::string sin_ref;
  std::string cos_ref;
  std::string tan_ref;

  TrigContext()
      : condition(0.0), sin_component_relerr(0.0), cos_component_relerr(0.0),
        tan_via_sin_cos_relerr(0.0) {}
};

void print_usage() {
  std::cout << "oracle_test_unary [-dd] [-td] [-qd] [-edd] [-all] [-v]"
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
base_diag(const char *case_name, int iteration) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  std::ostringstream replay;
  replay << "tests/oracle/test_unary -" << traits::name()
         << " --seed=" << qd_oracle::rng::active_seed();
  diag.push_back(qd_oracle::Tap::Diagnostic(
      "seed", std::to_string(qd_oracle::rng::active_seed())));
  diag.push_back(qd_oracle::Tap::Diagnostic("replay", replay.str()));
  diag.push_back(qd_oracle::Tap::Diagnostic("case", case_name));
  diag.push_back(qd_oracle::Tap::Diagnostic("iteration", int_text(iteration)));
  return diag;
}

bool is_trig_op(qd_oracle::UnaryOp op) {
  return op == qd_oracle::unary_sin || op == qd_oracle::unary_cos ||
         op == qd_oracle::unary_tan;
}

bool is_conditioned_trig_domain(qd_oracle::InputDomain domain) {
  return domain == qd_oracle::domain_trig_conditioned;
}

template <class T>
T apply_unary(qd_oracle::UnaryOp op, const T &a) {
  switch (op) {
  case qd_oracle::unary_sqrt:
    return sqrt(a);
  case qd_oracle::unary_sqr:
    return sqr(a);
  case qd_oracle::unary_exp:
    return exp(a);
  case qd_oracle::unary_log:
    return log(a);
  case qd_oracle::unary_log10:
    return log10(a);
  case qd_oracle::unary_sin:
    return sin(a);
  case qd_oracle::unary_cos:
    return cos(a);
  case qd_oracle::unary_tan:
    return tan(a);
  case qd_oracle::unary_asin:
    return asin(a);
  case qd_oracle::unary_acos:
    return acos(a);
  case qd_oracle::unary_atan:
    return atan(a);
  case qd_oracle::unary_sinh:
    return sinh(a);
  case qd_oracle::unary_cosh:
    return cosh(a);
  case qd_oracle::unary_tanh:
    return tanh(a);
  }
  return T(0);
}

void apply_mpfr(qd_oracle::UnaryOp op, mpfr_t ref, mpfr_t a) {
  switch (op) {
  case qd_oracle::unary_sqrt:
    mpfr_sqrt(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_sqr:
    mpfr_sqr(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_exp:
    mpfr_exp(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_log:
    mpfr_log(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_log10:
    mpfr_log10(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_sin:
    mpfr_sin(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_cos:
    mpfr_cos(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_tan:
    mpfr_tan(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_asin:
    mpfr_asin(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_acos:
    mpfr_acos(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_atan:
    mpfr_atan(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_sinh:
    mpfr_sinh(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_cosh:
    mpfr_cosh(ref, a, MPFR_RNDN);
    break;
  case qd_oracle::unary_tanh:
    mpfr_tanh(ref, a, MPFR_RNDN);
    break;
  }
}

template <class T>
void trig_refs(mpfr_t sin_ref, mpfr_t cos_ref, mpfr_t tan_ref, mpfr_t input) {
  mpfr_sin(sin_ref, input, MPFR_RNDN);
  mpfr_cos(cos_ref, input, MPFR_RNDN);
  mpfr_div(tan_ref, sin_ref, cos_ref, MPFR_RNDN);
}

bool stable_trig_refs(mpfr_t sin_ref, mpfr_t cos_ref) {
  const double sin_abs = std::fabs(mpfr_get_d(sin_ref, MPFR_RNDN));
  const double cos_abs = std::fabs(mpfr_get_d(cos_ref, MPFR_RNDN));
  return sin_abs >= kStableTrigThreshold && cos_abs >= kStableTrigThreshold;
}

double trig_condition_number(qd_oracle::UnaryOp op, mpfr_t input,
                             mpfr_t sin_ref, mpfr_t cos_ref) {
  mpfr_t tmp;
  mpfr_t denom;
  mpfr_inits2(mpfr_get_prec(input), tmp, denom, (mpfr_ptr) 0);
  mpfr_abs(tmp, input, MPFR_RNDN);

  switch (op) {
  case qd_oracle::unary_sin:
    mpfr_mul(tmp, tmp, cos_ref, MPFR_RNDN);
    mpfr_abs(tmp, tmp, MPFR_RNDN);
    mpfr_abs(denom, sin_ref, MPFR_RNDN);
    break;
  case qd_oracle::unary_cos:
    mpfr_mul(tmp, tmp, sin_ref, MPFR_RNDN);
    mpfr_abs(tmp, tmp, MPFR_RNDN);
    mpfr_abs(denom, cos_ref, MPFR_RNDN);
    break;
  case qd_oracle::unary_tan:
    mpfr_mul(denom, sin_ref, cos_ref, MPFR_RNDN);
    mpfr_abs(denom, denom, MPFR_RNDN);
    break;
  default:
    mpfr_set_ui(denom, 1, MPFR_RNDN);
    break;
  }

  double condition = 0.0;
  if (mpfr_zero_p(denom)) {
    condition = HUGE_VAL;
  } else {
    mpfr_div(tmp, tmp, denom, MPFR_RNDN);
    condition = mpfr_get_d(tmp, MPFR_RNDN);
  }

  mpfr_clears(tmp, denom, (mpfr_ptr) 0);
  return condition;
}

template <class T>
TrigContext make_trig_context(qd_oracle::UnaryOp op, const T &input,
                              mpfr_t input_mp) {
  TrigContext context;
  mpfr_t sin_ref;
  mpfr_t cos_ref;
  mpfr_t tan_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), sin_ref, cos_ref, tan_ref,
              (mpfr_ptr) 0);
  trig_refs<T>(sin_ref, cos_ref, tan_ref, input_mp);

  context.condition = trig_condition_number(op, input_mp, sin_ref, cos_ref);
  context.sin_ref = qd_oracle::mpfr_to_string(sin_ref);
  context.cos_ref = qd_oracle::mpfr_to_string(cos_ref);
  context.tan_ref = qd_oracle::mpfr_to_string(tan_ref);

  const T got_sin = sin(input);
  const T got_cos = cos(input);
  const T got_tan_via_sin_cos = got_sin / got_cos;
  context.sin_component_relerr = qd_oracle::relerr_in_eps(got_sin, sin_ref);
  context.cos_component_relerr = qd_oracle::relerr_in_eps(got_cos, cos_ref);
  context.tan_via_sin_cos_relerr =
      qd_oracle::relerr_in_eps(got_tan_via_sin_cos, tan_ref);

  mpfr_clears(sin_ref, cos_ref, tan_ref, (mpfr_ptr) 0);
  return context;
}

template <class T>
T make_trig_input(qd_oracle::InputDomain domain) {
  T fallback = T(0);
  mpfr_t input_mp;
  mpfr_t sin_ref;
  mpfr_t cos_ref;
  mpfr_t tan_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), input_mp, sin_ref, cos_ref, tan_ref,
              (mpfr_ptr) 0);

  for (int attempt = 0; attempt < 256; ++attempt) {
    T candidate = qd_oracle::rng::uniform_type<T>(-4, 4);
    fallback = candidate;
    qd_oracle::to_mpfr(input_mp, candidate);
    trig_refs<T>(sin_ref, cos_ref, tan_ref, input_mp);
    const bool stable = stable_trig_refs(sin_ref, cos_ref);
    if (domain == qd_oracle::domain_trig_stable && stable) {
      mpfr_clears(input_mp, sin_ref, cos_ref, tan_ref, (mpfr_ptr) 0);
      return candidate;
    }
    if (domain == qd_oracle::domain_trig_conditioned && !stable) {
      mpfr_clears(input_mp, sin_ref, cos_ref, tan_ref, (mpfr_ptr) 0);
      return candidate;
    }
  }

  mpfr_clears(input_mp, sin_ref, cos_ref, tan_ref, (mpfr_ptr) 0);
  return fallback;
}

template <class T>
T make_input(qd_oracle::InputDomain domain) {
  switch (domain) {
  case qd_oracle::domain_all_moderate:
    return qd_oracle::rng::uniform_type<T>(-24, 24);
  case qd_oracle::domain_nonnegative:
    return qd_oracle::rng::positive_type<T>(-160, 160);
  case qd_oracle::domain_positive:
    return qd_oracle::rng::positive_type<T>(-80, 80);
  case qd_oracle::domain_unit:
    return qd_oracle::rng::unit_type<T>();
  case qd_oracle::domain_exp_moderate:
    return qd_oracle::rng::uniform_type<T>(-3, 3);
  case qd_oracle::domain_trig_moderate:
    return qd_oracle::rng::uniform_type<T>(-4, 4);
  case qd_oracle::domain_trig_stable:
  case qd_oracle::domain_trig_conditioned:
    return make_trig_input<T>(domain);
  case qd_oracle::domain_hyperbolic_moderate:
    return qd_oracle::rng::uniform_type<T>(-2, 2);
  }
  return T(0);
}

void append_trig_diag(std::vector<qd_oracle::Tap::Diagnostic> *diag,
                      const TrigContext &context) {
  diag->push_back(qd_oracle::Tap::Diagnostic(
      "condition_number", double_text(context.condition)));
  diag->push_back(qd_oracle::Tap::Diagnostic("mpfr_sin", context.sin_ref));
  diag->push_back(qd_oracle::Tap::Diagnostic("mpfr_cos", context.cos_ref));
  diag->push_back(qd_oracle::Tap::Diagnostic("mpfr_tan", context.tan_ref));
  diag->push_back(qd_oracle::Tap::Diagnostic(
      "sin_component_relerr_eps", double_text(context.sin_component_relerr)));
  diag->push_back(qd_oracle::Tap::Diagnostic(
      "cos_component_relerr_eps", double_text(context.cos_component_relerr)));
  diag->push_back(qd_oracle::Tap::Diagnostic(
      "tan_via_sin_cos_relerr_eps",
      double_text(context.tan_via_sin_cos_relerr)));
}

std::string trig_comment_suffix(const TrigContext &context) {
  std::ostringstream os;
  os << " condition_number=" << context.condition
     << " sin_ref=" << context.sin_ref
     << " cos_ref=" << context.cos_ref
     << " tan_ref=" << context.tan_ref
     << " sin_component_relerr_eps=" << context.sin_component_relerr
     << " cos_component_relerr_eps=" << context.cos_component_relerr
     << " tan_via_sin_cos_relerr_eps=" << context.tan_via_sin_cos_relerr;
  return os.str();
}

bool worse_relerr(double relerr, double worst) {
  return relerr > worst || !std::isfinite(relerr);
}

template <class T>
bool run_unary_case(qd_oracle::Tap &tap,
                    const qd_oracle::UnaryRegistryEntry &entry,
                    bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;

  bool pass = true;
  double worst = -1.0;
  double worst_allowed = entry.bound.eps_multiplier;
  int worst_iteration = -1;
  T worst_input;
  T worst_got;
  TrigContext worst_trig_context;
  mpfr_t worst_ref;
  mpfr_init2(worst_ref, qd_oracle::ref_prec<T>());

  mpfr_t input_mp;
  mpfr_t ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), input_mp, ref, (mpfr_ptr) 0);

  for (int i = 0; i < kSamples; ++i) {
    T input = make_input<T>(entry.domain);
    T got = apply_unary(entry.op, input);
    qd_oracle::to_mpfr(input_mp, input);
    apply_mpfr(entry.op, ref, input_mp);

    TrigContext trig_context;
    double allowed = entry.bound.eps_multiplier;
    if (is_trig_op(entry.op)) {
      trig_context = make_trig_context<T>(entry.op, input, input_mp);
      if (is_conditioned_trig_domain(entry.domain)) {
        allowed *= std::max(1.0, trig_context.condition);
      }
    }

    double relerr = qd_oracle::relerr_in_eps(got, ref);
    if (worse_relerr(relerr, worst)) {
      worst = relerr;
      worst_allowed = allowed;
      worst_iteration = i;
      worst_input = input;
      worst_got = got;
      worst_trig_context = trig_context;
      mpfr_set(worst_ref, ref, MPFR_RNDN);
    }
    if (!std::isfinite(relerr) || relerr > allowed) {
      pass = false;
    }
  }

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) {
    diag = base_diag<T>(entry.name, worst_iteration);
    diag.push_back(qd_oracle::Tap::Diagnostic("input_limbs",
                                              qd_oracle::limbs_hex(worst_input)));
    diag.push_back(qd_oracle::Tap::Diagnostic("mpfr_reference",
                                              qd_oracle::mpfr_to_string(worst_ref)));
    diag.push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                              qd_oracle::limbs_hex(worst_got)));
    diag.push_back(qd_oracle::Tap::Diagnostic("relerr_eps",
                                              double_text(worst)));
    diag.push_back(qd_oracle::Tap::Diagnostic(
        "allowed_eps_multiplier", double_text(worst_allowed)));
    diag.push_back(qd_oracle::Tap::Diagnostic("bound_justification",
                                              entry.bound.justification));
    if (is_trig_op(entry.op)) {
      append_trig_diag(&diag, worst_trig_context);
    }
  }

  std::string name = std::string(traits::name()) + " " + entry.name +
                     " random oracle";
  tap.ok(pass, name, diag);
  if (verbose || pass) {
    std::cout << "# " << traits::name() << " " << entry.name
              << " worst_relerr_eps=" << worst
              << " allowed_eps=" << worst_allowed
              << " samples=" << kSamples
              << " worst_input_limbs=" << qd_oracle::limbs_hex(worst_input);
    if (is_trig_op(entry.op)) {
      std::cout << trig_comment_suffix(worst_trig_context);
    }
    std::cout << "\n";
  }

  mpfr_clears(input_mp, ref, (mpfr_ptr) 0);
  mpfr_clear(worst_ref);
  return pass;
}

template <class T>
bool run_type(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  std::size_t count = 0;
  const qd_oracle::UnaryRegistryEntry *entries = qd_oracle::unary_registry(&count);
  for (std::size_t i = 0; i < count; ++i) {
    pass &= run_unary_case<T>(tap, entries[i], verbose);
  }
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
  std::size_t count = 0;
  qd_oracle::unary_registry(&count);
  return static_cast<int>(count);
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
