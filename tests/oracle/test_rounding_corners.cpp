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
#include <vector>

#include <qd/fpu.h>

namespace {

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

struct RoundingCaseRecord {
  std::string type;
  std::string input_id;
  std::string operation;
  std::string variant;
  double relerr;
  double allowed;
  bool pass;
};

void print_usage() {
  std::cout << "oracle_test_rounding_corners [-dd|--dd] [-td|--td]"
            << " [-qd|--qd] [-edd|--edd] [-all|--all] [-v|--verbose]"
            << " [--seed=N] [--worst-report=FILE]" << std::endl;
}

const char *arith_op_id(qd_oracle::ArithOp op) {
  switch (op) {
  case qd_oracle::arith_add:
    return "add";
  case qd_oracle::arith_sub:
    return "sub";
  case qd_oracle::arith_mul:
    return "mul";
  case qd_oracle::arith_div:
    return "div";
  case qd_oracle::arith_sqr:
    return "sqr";
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
    if (*it == ' ' || *it == '/' || *it == '-') {
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

void append_record(std::vector<RoundingCaseRecord> *records,
                  const char *type_name, const char *case_name,
                  qd_oracle::ArithOp op, double relerr, double allowed,
                  bool pass) {
  if (!records) {
    return;
  }
  RoundingCaseRecord record;
  record.type = type_name;
  record.input_id = csv_safe(case_name);
  record.operation = arith_op_id(op);
  record.variant = variant_signature();
  record.relerr = relerr;
  record.allowed = allowed;
  record.pass = pass;
  records->push_back(record);
}

void emit_worst_report(const std::string &path,
                      const std::vector<RoundingCaseRecord> &records) {
  std::ofstream out(path.c_str());
  if (!out) {
    std::cerr << "Failed to open worst report: " << path << "\n";
    std::exit(1);
  }

  out << "build_variant,type,input_id,operation,relerr,allowed,pass,seed,ieee_add,sloppy_mul,sloppy_div,fma\n";
  for (std::size_t i = 0; i < records.size(); ++i) {
    const RoundingCaseRecord &r = records[i];
    const std::string ieee_add = extract_flag_from_variant(r.variant, "ieee_add=");
    const std::string sloppy_mul = extract_flag_from_variant(r.variant, "sloppy_mul=");
    const std::string sloppy_div = extract_flag_from_variant(r.variant, "sloppy_div=");
    const std::string fma = extract_flag_from_variant(r.variant, "fma=");

    out << r.variant << "," << r.type << "," << r.input_id << "," << r.operation
        << "," << std::setprecision(17) << r.relerr << "," << std::setprecision(17)
        << r.allowed << "," << (r.pass ? "pass" : "fail") << ","
        << qd_oracle::rng::active_seed() << "," << ieee_add << ","
        << sloppy_mul << "," << sloppy_div << "," << fma << "\n";
  }
}

std::string double_text(double value) {
  std::ostringstream os;
  os << value;
  return os.str();
}

template <class T>
std::vector<qd_oracle::Tap::Diagnostic> base_diag(const char *case_name) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  std::ostringstream replay;
  replay << "tests/oracle/test_rounding_corners -" << traits::name()
         << " --seed=" << qd_oracle::rng::active_seed();
  diag.push_back(qd_oracle::Tap::Diagnostic("seed",
                                             std::to_string(qd_oracle::rng::active_seed())));
  diag.push_back(qd_oracle::Tap::Diagnostic("replay", replay.str()));
  diag.push_back(qd_oracle::Tap::Diagnostic("case", case_name));
  return diag;
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

void apply_mpfr(qd_oracle::ArithOp op, mpfr_t out, mpfr_t a, mpfr_t b) {
  switch (op) {
  case qd_oracle::arith_add:
    mpfr_add(out, a, b, MPFR_RNDN);
    break;
  case qd_oracle::arith_sub:
    mpfr_sub(out, a, b, MPFR_RNDN);
    break;
  case qd_oracle::arith_mul:
    mpfr_mul(out, a, b, MPFR_RNDN);
    break;
  case qd_oracle::arith_div:
    mpfr_div(out, a, b, MPFR_RNDN);
    break;
  case qd_oracle::arith_sqr:
    mpfr_sqr(out, a, MPFR_RNDN);
    break;
  }
}

template <class T>
bool report_case(qd_oracle::Tap &tap, const char *case_name,
                 bool pass,
                 const std::vector<qd_oracle::Tap::Diagnostic> &diag,
                 bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;
  tap.ok(pass, std::string(traits::name()) + " " + case_name, diag,
         verbose);
  return pass;
}

template <class T>
bool exact_case(qd_oracle::Tap &tap, const char *case_name,
                qd_oracle::ArithOp op, const T &a, const T &b,
                std::vector<RoundingCaseRecord> *records, bool verbose) {
  bool pass = true;
  mpfr_t a_mp;
  mpfr_t b_mp;
  mpfr_t exact_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, exact_ref, (mpfr_ptr) 0);
  qd_oracle::to_mpfr(a_mp, a);
  qd_oracle::to_mpfr(b_mp, b);
  apply_mpfr(op, exact_ref, a_mp, b_mp);

  T got = apply_arith(op, a, b);
  T expected = qd_oracle::round_nearest_expansion<T>(exact_ref);

  const double relerr = qd_oracle::relerr_in_eps(got, exact_ref);
  pass = qd_oracle::bit_equal(got, expected);
  append_record(records, qd_oracle::TypeTraits<T>::name(), case_name, op, relerr,
                0.0, pass);

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = base_diag<T>(case_name);
    diag.push_back(qd_oracle::Tap::Diagnostic("input_a_limbs",
                                              qd_oracle::limbs_hex(a)));
    diag.push_back(qd_oracle::Tap::Diagnostic("input_a_value",
                                              qd_oracle::value_to_mpfr_string(a)));
    if (op == qd_oracle::arith_sqr) {
      diag.push_back(qd_oracle::Tap::Diagnostic("input_b_unused", "true"));
    } else {
      diag.push_back(qd_oracle::Tap::Diagnostic("input_b_limbs",
                                                qd_oracle::limbs_hex(b)));
      diag.push_back(qd_oracle::Tap::Diagnostic("input_b_value",
                                                qd_oracle::value_to_mpfr_string(b)));
    }
    diag.push_back(qd_oracle::Tap::Diagnostic("mpfr_reference",
                                              qd_oracle::mpfr_to_string(exact_ref)));
    diag.push_back(qd_oracle::Tap::Diagnostic("got_value",
                                              qd_oracle::value_to_mpfr_string(got)));
    diag.push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                              qd_oracle::limbs_hex(got)));
    diag.push_back(qd_oracle::Tap::Diagnostic("expected_limbs",
                                              qd_oracle::limbs_hex(expected)));
    diag.push_back(qd_oracle::Tap::Diagnostic("abs_error_mpfr",
                                              qd_oracle::abs_error_to_string(got, exact_ref)));
    diag.push_back(qd_oracle::Tap::Diagnostic("relerr_eps",
                                              double_text(relerr)));
    diag.push_back(qd_oracle::Tap::Diagnostic("ulp_error_estimate",
                                              double_text(qd_oracle::ulp_error_estimate(got, exact_ref))));
  }
  pass = pass && qd_oracle::bit_equal(expected, got);
  mpfr_clears(a_mp, b_mp, exact_ref, (mpfr_ptr) 0);
  return report_case<T>(tap, case_name, pass, diag, verbose);
}


template <class T>
bool bounded_case(qd_oracle::Tap &tap, const char *case_name,
                 qd_oracle::ArithOp op, const T &a, const T &b,
                 const qd_oracle::ErrorBound &bound,
                 std::vector<RoundingCaseRecord> *records, bool verbose) {
  bool pass = true;
  mpfr_t a_mp;
  mpfr_t b_mp;
  mpfr_t exact_ref;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, exact_ref, (mpfr_ptr) 0);
  qd_oracle::to_mpfr(a_mp, a);
  qd_oracle::to_mpfr(b_mp, b);
  apply_mpfr(op, exact_ref, a_mp, b_mp);

  T got = apply_arith(op, a, b);
  double relerr = qd_oracle::relerr_in_eps(got, exact_ref);
  if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) {
    pass = false;
  }
  append_record(records, qd_oracle::TypeTraits<T>::name(), case_name, op, relerr,
                bound.eps_multiplier, pass);

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = base_diag<T>(case_name);
    diag.push_back(qd_oracle::Tap::Diagnostic("input_a_limbs",
                                              qd_oracle::limbs_hex(a)));
    diag.push_back(qd_oracle::Tap::Diagnostic("input_a_value",
                                              qd_oracle::value_to_mpfr_string(a)));
    if (op == qd_oracle::arith_sqr) {
      diag.push_back(qd_oracle::Tap::Diagnostic("input_b_unused", "true"));
    } else {
      diag.push_back(qd_oracle::Tap::Diagnostic("input_b_limbs",
                                                qd_oracle::limbs_hex(b)));
      diag.push_back(qd_oracle::Tap::Diagnostic("input_b_value",
                                                qd_oracle::value_to_mpfr_string(b)));
    }
    diag.push_back(qd_oracle::Tap::Diagnostic("mpfr_reference",
                                              qd_oracle::mpfr_to_string(exact_ref)));
    diag.push_back(qd_oracle::Tap::Diagnostic("got_value",
                                              qd_oracle::value_to_mpfr_string(got)));
    diag.push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                              qd_oracle::limbs_hex(got)));
    diag.push_back(qd_oracle::Tap::Diagnostic("abs_error_mpfr",
                                              qd_oracle::abs_error_to_string(got, exact_ref)));
    diag.push_back(qd_oracle::Tap::Diagnostic("relerr_eps",
                                              double_text(relerr)));
    diag.push_back(qd_oracle::Tap::Diagnostic("ulp_error_estimate",
                                              double_text(qd_oracle::ulp_error_estimate(got, exact_ref))));
    diag.push_back(qd_oracle::Tap::Diagnostic("allowed_eps_multiplier",
                                              double_text(bound.eps_multiplier)));
    diag.push_back(qd_oracle::Tap::Diagnostic("bound_justification",
                                              bound.justification));
  }

  mpfr_clears(a_mp, b_mp, exact_ref, (mpfr_ptr) 0);
  return report_case<T>(tap, case_name, pass, diag, verbose);
}

// Class 1: unconditional exactness (non-overlapping arithmetic and zero identities).
template <class T>
bool run_exactness_cases(qd_oracle::Tap &tap,
                        std::vector<RoundingCaseRecord> *records,
                        bool verbose) {
  bool pass = true;
  pass &= exact_case<T>(
      tap, "sterbenz same-sign subtraction exact", qd_oracle::arith_sub,
      T(1.5), T(0.75), records, verbose);
  pass &= exact_case<T>(
      tap, "non-overlapping powers-of-two sum exact",
      qd_oracle::arith_add, ldexp(T(1), 48), ldexp(T(1), -120), records, verbose);
  pass &= exact_case<T>(
      tap, "exact scaling by power-of-two multiplication",
      qd_oracle::arith_mul, T(1.25), ldexp(T(1), 37), records, verbose);
  pass &= exact_case<T>(
      tap, "exact square for low-precision input",
      qd_oracle::arith_sqr, ldexp(T(1), -20), T(0), records, verbose);
  pass &= exact_case<T>(
      tap, "positive signed zero for a+(-a)", qd_oracle::arith_add, T(1), -T(1),
      records, verbose);
  pass &= exact_case<T>(
      tap, "positive signed zero for a-a", qd_oracle::arith_sub, T(1), T(1), records, verbose);
  return pass;
}

// Class 2: tie / faithful rounding smoke tests.
template <class T>
bool run_tie_cases(qd_oracle::Tap &tap,
                   std::vector<RoundingCaseRecord> *records,
                   bool verbose) {
  bool pass = true;
  const T half_eps = T(T::_eps) / T(2);

  pass &= bounded_case<T>(
      tap, "add tie-guarded midpoint case",
      qd_oracle::arith_add, T(1), half_eps,
      qd_oracle::stable_trig_bound("tie_add"), records, verbose);
  pass &= bounded_case<T>(
      tap, "mul tie-guarded midpoint case",
      qd_oracle::arith_mul, T(1), T(1) + half_eps,
      qd_oracle::stable_trig_bound("tie_mul"), records, verbose);
  return pass;
}

// Class 3: variant distinguishing cases across build variants.
template <class T>
bool run_add_variant_case(qd_oracle::Tap &tap,
                        std::vector<RoundingCaseRecord> *records,
                        bool verbose) {
  const T a = T(1.0) + ldexp(T(1.0), -40);
  const T b = -ldexp(T(1.0), -40) + ldexp(T(1.0), -100);

#if defined(QD_IEEE_ADD)
  return exact_case<T>(
      tap, "variant-distinguishing carry-ripple cancellation (IEEE)"
      " exact",
      qd_oracle::arith_add, a, b, records, verbose);
#else
  return bounded_case<T>(
      tap, "variant-distinguishing carry-ripple cancellation (Cray)",
      qd_oracle::arith_add, a, b,
      qd_oracle::add_sub_bound("carry_ripple"), records, verbose);
#endif
}

template <class T>
bool run_mul_variant_case(qd_oracle::Tap &tap,
                        std::vector<RoundingCaseRecord> *records,
                        bool verbose) {
  const T a = (T(1) + ldexp(T(1), -60)) + ldexp(T(1), -120);
  const T b = (T(1) + ldexp(T(1), -70)) + ldexp(T(1), -130);

#if defined(QD_SLOPPY_MUL)
  return bounded_case<T>(
      tap, "qd mul variant distinguishing (sloppy)",
      qd_oracle::arith_mul, a, b, qd_oracle::mul_bound("mul_variant"), records, verbose);
#else
  return exact_case<T>(
      tap, "qd mul variant distinguishing (accurate) exact",
      qd_oracle::arith_mul, a, b, records, verbose);
#endif
}

template <class T>
bool run_type(qd_oracle::Tap &tap,
              std::vector<RoundingCaseRecord> *records, bool verbose) {
  bool pass = true;
  pass &= run_exactness_cases<T>(tap, records, verbose);
  pass &= run_tie_cases<T>(tap, records, verbose);
  pass &= run_add_variant_case<T>(tap, records, verbose);
  return pass;
}

template <>
bool run_type<qd_real>(qd_oracle::Tap &tap,
                      std::vector<RoundingCaseRecord> *records, bool verbose) {
  bool pass = true;
  pass &= run_exactness_cases<qd_real>(tap, records, verbose);
  pass &= run_tie_cases<qd_real>(tap, records, verbose);
  pass &= run_add_variant_case<qd_real>(tap, records, verbose);
  pass &= run_mul_variant_case<qd_real>(tap, records, verbose);
  return pass;
}

int selected_count(const Options &options) {
  int count = 0;
  if (options.test_dd) count += 9;
  if (options.test_td) count += 9;
  if (options.test_qd) count += 10;
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) count += 9;
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
        std::strcmp(argv[i], "-help") == 0 ||
        std::strcmp(argv[i], "--help") == 0) {
      print_usage();
      std::exit(0);
    } else if (std::strcmp(argv[i], "-dd") == 0 ||
               std::strcmp(argv[i], "--dd") == 0) {
      options->test_dd = true;
    } else if (std::strcmp(argv[i], "-td") == 0 ||
               std::strcmp(argv[i], "--td") == 0) {
      options->test_td = true;
    } else if (std::strcmp(argv[i], "-qd") == 0 ||
               std::strcmp(argv[i], "--qd") == 0) {
      options->test_qd = true;
    } else if (std::strcmp(argv[i], "-edd") == 0 ||
               std::strcmp(argv[i], "--edd") == 0) {
#ifdef QD_HAVE_EDD_REAL
      options->test_edd = true;
#else
      std::cerr << "edd_real is not enabled in this build\n";
      return false;
#endif
    } else if (std::strcmp(argv[i], "-all") == 0 ||
               std::strcmp(argv[i], "--all") == 0) {
      select_all(options);
    } else if (std::strcmp(argv[i], "-v") == 0 ||
               std::strcmp(argv[i], "-verbose") == 0 ||
               std::strcmp(argv[i], "--verbose") == 0) {
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

  std::vector<RoundingCaseRecord> records;
  std::vector<RoundingCaseRecord> *record_sink =
      options.worst_report.empty() ? 0 : &records;

  qd_oracle::rng::configure(options.has_seed, options.seed);
  qd_oracle::Tap tap(selected_count(options));
  std::cout << "# seed: " << qd_oracle::rng::active_seed() << "\n";

  bool pass = true;
  unsigned int old_cw = 0;
  bool fpu_fixed = true;
  fpu_fix_start(&old_cw);

  if (options.test_dd) {
    pass &= run_type<dd_real>(tap, record_sink, options.verbose);
  }
  if (options.test_td) {
    pass &= run_type<td_real>(tap, record_sink, options.verbose);
  }
  if (options.test_qd) {
    pass &= run_type<qd_real>(tap, record_sink, options.verbose);
  }
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) {
    fpu_fix_end(&old_cw);
    fpu_fixed = false;
    pass &= run_type<edd_real>(tap, record_sink, options.verbose);
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
