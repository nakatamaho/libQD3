#include "fn_registry.h"
#include "mpfr_oracle.h"
#include "qd_rng.h"
#include "tap.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>
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
  std::string worst_report;

  Options()
      : test_dd(false), test_td(false), test_qd(false), test_edd(false),
        verbose(false), has_seed(false), seed(0), worst_report("") {}
};

struct TrigContext {
  double condition;
  double sin_component_relerr;
  double cos_component_relerr;
  double tan_via_sin_cos_relerr;
  int td_j_sector;
  int td_k_sector;
  std::string sin_ref;
  std::string cos_ref;
  std::string tan_ref;

  TrigContext()
      : condition(0.0), sin_component_relerr(0.0), cos_component_relerr(0.0),
        tan_via_sin_cos_relerr(0.0), td_j_sector(999), td_k_sector(999) {}
};

struct UnaryCaseRecord {
  std::string type;
  std::string function_name;
  std::string operation;
  std::string domain;
  std::string build_variant;
  std::string input_mpfr;
  std::string input_limbs;
  double relerr;
  double allowed;
  bool pass;
  double condition;
  std::string sin_ref;
  std::string cos_ref;
  std::string tan_ref;
  double sin_component_relerr;
  double cos_component_relerr;
  double tan_via_sin_cos_relerr;
  int td_j_sector;
  int td_k_sector;
};

const char *unary_op_id(qd_oracle::UnaryOp op) {
  switch (op) {
  case qd_oracle::unary_sqrt:
    return "sqrt";
  case qd_oracle::unary_sqr:
    return "sqr";
  case qd_oracle::unary_exp:
    return "exp";
  case qd_oracle::unary_log:
    return "log";
  case qd_oracle::unary_log10:
    return "log10";
  case qd_oracle::unary_sin:
    return "sin";
  case qd_oracle::unary_cos:
    return "cos";
  case qd_oracle::unary_tan:
    return "tan";
  case qd_oracle::unary_asin:
    return "asin";
  case qd_oracle::unary_acos:
    return "acos";
  case qd_oracle::unary_atan:
    return "atan";
  case qd_oracle::unary_sinh:
    return "sinh";
  case qd_oracle::unary_cosh:
    return "cosh";
  case qd_oracle::unary_tanh:
    return "tanh";
  }
  return "unknown";
}

const char *input_domain_name(qd_oracle::InputDomain domain) {
  switch (domain) {
  case qd_oracle::domain_all_moderate:
    return "all_moderate";
  case qd_oracle::domain_nonnegative:
    return "nonnegative";
  case qd_oracle::domain_positive:
    return "positive";
  case qd_oracle::domain_unit:
    return "unit";
  case qd_oracle::domain_exp_moderate:
    return "exp_moderate";
  case qd_oracle::domain_trig_moderate:
    return "trig_moderate";
  case qd_oracle::domain_trig_stable:
    return "trig_stable";
  case qd_oracle::domain_trig_conditioned:
    return "trig_conditioned";
  case qd_oracle::domain_hyperbolic_moderate:
    return "hyperbolic_moderate";
  }
  return "unknown";
}

std::string variant_signature() {
  std::ostringstream os;
#ifdef QD_IEEE_ADD
  os << "ieee_add=1";
#else
  os << "ieee_add=0";
#endif
#ifdef QD_SLOPPY_MUL
  os << "|sloppy_mul=1";
#else
  os << "|sloppy_mul=0";
#endif
#ifdef QD_SLOPPY_DIV
  os << "|sloppy_div=1";
#else
  os << "|sloppy_div=0";
#endif
#ifdef QD_FMA
  os << "|fma=1";
#else
  os << "|fma=0";
#endif
  return os.str();
}

std::string csv_safe(std::string value) {
  for (std::string::iterator it = value.begin(); it != value.end(); ++it) {
    if (*it == ' ' || *it == '/' || *it == '-' || *it == ',') {
      *it = '_';
    }
  }
  return value;
}

std::string extract_flag_from_variant(const std::string &variant,
                                     const char *key) {
  const std::size_t pos = variant.find(key);
  if (pos == std::string::npos) {
    return "0";
  }
  const std::size_t value_pos = pos + std::strlen(key);
  if (value_pos >= variant.size()) {
    return "0";
  }
  const char value = variant[value_pos];
  if (value == '0' || value == '1') {
    return std::string(1, value);
  }
  return "0";
}

void print_usage() {
  std::cout << "oracle_test_unary [-dd] [-td] [-qd] [-edd] [-all] [-v]"
            << " [--seed=N] [--worst-report=FILE]\n";
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
void set_td_trig_sector(TrigContext *context, const T &input);
void set_td_trig_sector(TrigContext *context, const td_real &input);

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
  set_td_trig_sector(&context, input);

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


void td_reduce_trig_arg_for_oracle(const td_real &a, int &j, int &k) {
  static const td_real td_pi16(1.963495408493620697e-01,
                               7.654042494670957545e-18,
                               -1.871731131073962291e-34);

  td_real z = nint(a / td_real::_2pi);
  td_real r = a - td_real::_2pi * z;

  td_real q = nint(r / td_real::_pi2);
  td_real t = r - td_real::_pi2 * q;
  j = static_cast<int>(q[0]);
  while (j > 2) j -= 4;
  while (j < -2) j += 4;

  q = nint(t / td_pi16);
  k = static_cast<int>(q[0]);
}

template <class T>
void set_td_trig_sector(TrigContext *context, const T &) {
  context->td_j_sector = 999;
  context->td_k_sector = 999;
}

void set_td_trig_sector(TrigContext *context, const td_real &input) {
  td_reduce_trig_arg_for_oracle(input, context->td_j_sector, context->td_k_sector);
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
  diag->push_back(qd_oracle::Tap::Diagnostic("td_j_sector",
                                            int_text(context.td_j_sector)));
  diag->push_back(qd_oracle::Tap::Diagnostic("td_k_sector",
                                            int_text(context.td_k_sector)));
}

std::string trig_comment_suffix(const TrigContext &context) {
  std::ostringstream os;
  os << " condition_number=" << context.condition
     << " sin_ref=" << context.sin_ref
     << " cos_ref=" << context.cos_ref
     << " tan_ref=" << context.tan_ref
     << " sin_component_relerr_eps=" << context.sin_component_relerr
     << " cos_component_relerr_eps=" << context.cos_component_relerr
     << " tan_via_sin_cos_relerr_eps=" << context.tan_via_sin_cos_relerr
     << " td_j=" << context.td_j_sector
     << " td_k=" << context.td_k_sector;
  return os.str();
}

template <class T>
void append_record(std::vector<UnaryCaseRecord> *records,
                  const qd_oracle::UnaryRegistryEntry &entry,
                  double relerr, double allowed, bool pass,
                  const T &worst_input, const std::string &worst_input_mpfr,
                  const TrigContext &trig_context) {
  if (!records) {
    return;
  }
  UnaryCaseRecord record;
  record.type = qd_oracle::TypeTraits<T>::name();
  record.function_name = csv_safe(entry.name);
  record.operation = unary_op_id(entry.op);
  record.domain = input_domain_name(entry.domain);
  record.build_variant = variant_signature();
  record.input_mpfr = worst_input_mpfr;
  record.input_limbs = qd_oracle::limbs_hex(worst_input);
  record.relerr = relerr;
  record.allowed = allowed;
  record.pass = pass;
  record.condition = trig_context.condition;
  record.sin_ref = trig_context.sin_ref;
  record.cos_ref = trig_context.cos_ref;
  record.tan_ref = trig_context.tan_ref;
  record.sin_component_relerr = trig_context.sin_component_relerr;
  record.cos_component_relerr = trig_context.cos_component_relerr;
  record.tan_via_sin_cos_relerr = trig_context.tan_via_sin_cos_relerr;
  record.td_j_sector = trig_context.td_j_sector;
  record.td_k_sector = trig_context.td_k_sector;
  records->push_back(record);
}

std::string csv_field(const std::string &value) {
  std::string out;
  out.reserve(value.size() + 2);
  out.push_back('"');
  for (std::size_t i = 0; i < value.size(); ++i) {
    const char c = value[i];
    if (c == '"') {
      out += "\"\"";
    } else {
      out += c;
    }
  }
  out.push_back('"');
  return out;
}

void emit_worst_report(const std::string &path,
                       const std::vector<UnaryCaseRecord> &records) {
  std::ofstream out(path.c_str());
  if (!out) {
    std::cerr << "Failed to open worst report: " << path << "\n";
    std::exit(1);
  }

  out << "build_variant,type,function,operation,domain,input_mpfr,input_limbs,relerr,allowed,pass,condition_number,mpfr_sin,mpfr_cos,mpfr_tan,sin_component_relerr_eps,cos_component_relerr_eps,tan_via_sin_cos_relerr_eps,td_j_sector,td_k_sector,seed,ieee_add,sloppy_mul,sloppy_div,fma\n";
  for (std::size_t i = 0; i < records.size(); ++i) {
    const UnaryCaseRecord &record = records[i];
    const std::string ieee_add = extract_flag_from_variant(record.build_variant, "ieee_add=");
    const std::string sloppy_mul = extract_flag_from_variant(record.build_variant, "sloppy_mul=");
    const std::string sloppy_div = extract_flag_from_variant(record.build_variant, "sloppy_div=");
    const std::string fma = extract_flag_from_variant(record.build_variant, "fma=");

    out << record.build_variant << "," << record.type << ","
        << record.function_name << "," << record.operation << ","
        << record.domain << "," << csv_field(record.input_mpfr) << ","
        << csv_field(record.input_limbs) << "," << std::setprecision(17)
        << record.relerr << "," << std::setprecision(17) << record.allowed
        << "," << (record.pass ? "pass" : "fail") << ","
        << record.condition << "," << record.sin_ref << ","
        << record.cos_ref << "," << record.tan_ref << ","
        << std::setprecision(17) << record.sin_component_relerr << ","
        << std::setprecision(17) << record.cos_component_relerr << ","
        << std::setprecision(17) << record.tan_via_sin_cos_relerr << ","
        << record.td_j_sector << "," << record.td_k_sector << ","
        << qd_oracle::rng::active_seed() << "," << ieee_add << ","
        << sloppy_mul << "," << sloppy_div << "," << fma << "\n";
  }
}

bool worse_relerr(double relerr, double worst) {
  return relerr > worst || !std::isfinite(relerr);
}

template <class T>
bool run_unary_case(qd_oracle::Tap &tap,
                    const qd_oracle::UnaryRegistryEntry &entry,
                    bool verbose,
                    std::vector<UnaryCaseRecord> *records) {
  typedef qd_oracle::TypeTraits<T> traits;

  bool pass = true;
  double worst = -1.0;
  double worst_allowed = entry.bound.eps_multiplier;
  int worst_iteration = -1;
  T worst_input;
  T worst_got;
  TrigContext worst_trig_context;
  mpfr_t worst_ref;
  std::string worst_input_mpfr;
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
      worst_input_mpfr = qd_oracle::mpfr_to_string(input_mp);
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
    diag.push_back(qd_oracle::Tap::Diagnostic("input_value",
                                              worst_input_mpfr));
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
              << " worst_input_limbs=" << qd_oracle::limbs_hex(worst_input)
              << " worst_input=" << worst_input_mpfr;
    if (is_trig_op(entry.op)) {
      std::cout << trig_comment_suffix(worst_trig_context);
    }
    std::cout << "\n";
  }

  append_record(records, entry, worst, worst_allowed, pass,
               worst_input, worst_input_mpfr, worst_trig_context);

  mpfr_clears(input_mp, ref, (mpfr_ptr) 0);
  mpfr_clear(worst_ref);
  return pass;
}

template <class T>
bool run_type(qd_oracle::Tap &tap, bool verbose,
             std::vector<UnaryCaseRecord> *records) {
  bool pass = true;
  std::size_t count = 0;
  const qd_oracle::UnaryRegistryEntry *entries = qd_oracle::unary_registry(&count);
  for (std::size_t i = 0; i < count; ++i) {
    pass &= run_unary_case<T>(tap, entries[i], verbose, records);
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
    } else if (std::strncmp(argv[i], "--worst-report=", 15) == 0) {
      options->worst_report = argv[i] + 15;
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

  if (options->worst_report.empty()) {
    const char *env_report = std::getenv("QD_ORACLE_WORST_REPORT");
    if (env_report != 0 && *env_report != '\0') {
      options->worst_report = env_report;
    }
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

  std::vector<UnaryCaseRecord> records;
  std::vector<UnaryCaseRecord> *record_sink =
      options.worst_report.empty() ? 0 : &records;

  qd_oracle::Tap tap(selected_count(options) * planned_per_type());
  std::cout << "# seed: " << qd_oracle::rng::active_seed() << "\n";

  bool pass = true;
  unsigned int old_cw = 0;
  bool fpu_fixed = true;
  fpu_fix_start(&old_cw);

  if (options.test_dd) {
    pass &= run_type<dd_real>(tap, options.verbose, record_sink);
  }
  if (options.test_td) {
    pass &= run_type<td_real>(tap, options.verbose, record_sink);
  }
  if (options.test_qd) {
    pass &= run_type<qd_real>(tap, options.verbose, record_sink);
  }

#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) {
    fpu_fix_end(&old_cw);
    fpu_fixed = false;
    pass &= run_type<edd_real>(tap, options.verbose, record_sink);
  }
#endif

  if (fpu_fixed) {
    fpu_fix_end(&old_cw);
  }

  if (!options.worst_report.empty()) {
    emit_worst_report(options.worst_report, records);
  }

  return pass ? tap.exit_status() : 1;
}
