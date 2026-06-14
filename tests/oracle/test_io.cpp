#include "mpfr_oracle.h"
#include "qd_rng.h"
#include "tap.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <ios>
#include <sstream>
#include <string>
#include <vector>

#include <qd/fpu.h>

namespace {

const int kSamples = 24;

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
  std::cout << "oracle_test_io [-dd] [-td] [-qd] [-edd] [-all] [-v]"
            << " [--seed=N]\n";
}

template <class T>
std::vector<qd_oracle::Tap::Diagnostic>
base_diag(const char *case_name, int iteration) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  std::ostringstream replay;
  replay << "tests/oracle/test_io -" << traits::name()
         << " --seed=" << qd_oracle::rng::active_seed();
  diag.push_back(qd_oracle::Tap::Diagnostic(
      "seed", std::to_string(qd_oracle::rng::active_seed())));
  diag.push_back(qd_oracle::Tap::Diagnostic("replay", replay.str()));
  diag.push_back(qd_oracle::Tap::Diagnostic("case", case_name));
  diag.push_back(qd_oracle::Tap::Diagnostic("iteration", int_text(iteration)));
  return diag;
}

template <class T>
void append_roundtrip_diag(std::vector<qd_oracle::Tap::Diagnostic> *diag,
                           const T &input, const std::string &text,
                           const T &parsed, double relerr,
                           double allowed_eps) {
  diag->push_back(qd_oracle::Tap::Diagnostic("input_limbs",
                                             qd_oracle::limbs_hex(input)));
  diag->push_back(qd_oracle::Tap::Diagnostic("text", text));
  diag->push_back(qd_oracle::Tap::Diagnostic("parsed_limbs",
                                             qd_oracle::limbs_hex(parsed)));
  diag->push_back(qd_oracle::Tap::Diagnostic("relerr_eps",
                                             double_text(relerr)));
  diag->push_back(qd_oracle::Tap::Diagnostic("allowed_eps_multiplier",
                                             double_text(allowed_eps)));
}

template <class T>
bool report_case(qd_oracle::Tap &tap, const char *case_name, bool pass,
                 double worst, int worst_iteration, const T &worst_input,
                 const std::string &worst_text, const T &worst_parsed,
                 double allowed_eps, bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) {
    diag = base_diag<T>(case_name, worst_iteration);
    append_roundtrip_diag(&diag, worst_input, worst_text, worst_parsed, worst,
                          allowed_eps);
  }
  std::string name = std::string(traits::name()) + " " + case_name;
  tap.ok(pass, name, diag);
  if (verbose || pass) {
    std::cout << "# " << traits::name() << " " << case_name
              << " worst_relerr_eps=" << worst
              << " allowed_eps=" << allowed_eps
              << " samples=" << kSamples << "\n";
  }
  return pass;
}

template <class T>
double relerr_between(const T &got, const T &ref_value) {
  mpfr_t ref;
  mpfr_init2(ref, qd_oracle::ref_prec<T>());
  qd_oracle::to_mpfr(ref, ref_value);
  double relerr = qd_oracle::relerr_in_eps(got, ref);
  mpfr_clear(ref);
  return relerr;
}

template <class T>
int full_precision() {
  return T::_ndigits + 8;
}

template <class T>
T make_io_input() {
  return qd_oracle::rng::uniform_type<T>(-8, 8);
}

template <class T>
bool run_scientific_roundtrip(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_input;
  T worst_parsed;
  std::string worst_text;
  const double allowed_eps = 2.0;

  for (int i = 0; i < kSamples; ++i) {
    T input = make_io_input<T>();
    std::string text = input.to_string(full_precision<T>(), 0,
                                       std::ios_base::scientific);
    T parsed(text.c_str());
    double relerr = relerr_between(parsed, input);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr;
      worst_iteration = i;
      worst_input = input;
      worst_parsed = parsed;
      worst_text = text;
    }
    if (!std::isfinite(relerr) || relerr > allowed_eps) {
      pass = false;
    }
  }

  return report_case(tap, "scientific full-precision bounded roundtrip", pass, worst,
                     worst_iteration, worst_input, worst_text, worst_parsed,
                     allowed_eps, verbose);
}

template <class T>
bool run_stream_roundtrip(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_input;
  T worst_parsed;
  std::string worst_text;
  const double allowed_eps = 2.0;

  for (int i = 0; i < kSamples; ++i) {
    T input = make_io_input<T>();
    std::ostringstream os;
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    os << std::setprecision(full_precision<T>()) << input;
    std::string text = os.str();
    T parsed;
    std::istringstream is(text);
    is >> parsed;
    double relerr = relerr_between(parsed, input);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr;
      worst_iteration = i;
      worst_input = input;
      worst_parsed = parsed;
      worst_text = text;
    }
    if (!is || !std::isfinite(relerr) || relerr > allowed_eps) {
      pass = false;
    }
  }

  return report_case(tap, "iostream full-precision bounded roundtrip", pass, worst,
                     worst_iteration, worst_input, worst_text, worst_parsed,
                     allowed_eps, verbose);
}

template <class T>
bool run_write_roundtrip(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_input;
  T worst_parsed;
  std::string worst_text;
  const double allowed_eps = 2.0;

  for (int i = 0; i < kSamples; ++i) {
    T input = make_io_input<T>();
    char buffer[256];
    input.write(buffer, static_cast<int>(sizeof(buffer)), full_precision<T>());
    std::string text(buffer);
    T parsed(text.c_str());
    double relerr = relerr_between(parsed, input);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr;
      worst_iteration = i;
      worst_input = input;
      worst_parsed = parsed;
      worst_text = text;
    }
    if (!std::isfinite(relerr) || relerr > allowed_eps) {
      pass = false;
    }
  }

  return report_case(tap, "write full-precision bounded roundtrip", pass, worst,
                     worst_iteration, worst_input, worst_text, worst_parsed,
                     allowed_eps, verbose);
}

template <class T>
bool run_reduced_precision(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  T worst_input;
  T worst_parsed;
  std::string worst_text;
  const int precision = std::max(6, T::_ndigits / 3);
  const double allowed_eps =
      std::pow(10.0, static_cast<double>(T::_ndigits - precision + 4));

  for (int i = 0; i < kSamples; ++i) {
    T input = make_io_input<T>();
    std::string text = input.to_string(precision, 0,
                                       std::ios_base::scientific);
    T parsed(text.c_str());
    double relerr = relerr_between(parsed, input);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr;
      worst_iteration = i;
      worst_input = input;
      worst_parsed = parsed;
      worst_text = text;
    }
    if (!std::isfinite(relerr) || relerr > allowed_eps) {
      pass = false;
    }
  }

  return report_case(tap, "reduced-precision bounded roundtrip", pass, worst,
                     worst_iteration, worst_input, worst_text, worst_parsed,
                     allowed_eps, verbose);
}

template <class T>
bool run_format_grid(qd_oracle::Tap &tap, bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  const T values[] = {T(0), T(1.25), T(-0.5), T(1024)};
  const int precisions[] = {0, 3, std::max(6, T::_ndigits / 2)};
  const int widths[] = {0, 24};
  const std::ios_base::fmtflags fmts[] = {
      static_cast<std::ios_base::fmtflags>(0), std::ios_base::scientific,
      std::ios_base::fixed};

  for (int vi = 0; vi < 4; ++vi) {
    for (int pi = 0; pi < 3; ++pi) {
      for (int wi = 0; wi < 2; ++wi) {
        for (int fi = 0; fi < 3; ++fi) {
          std::string text =
              values[vi].to_string(precisions[pi], widths[wi], fmts[fi], true,
                                   true, '_');
          if (text.empty() ||
              (widths[wi] > 0 && text.size() < static_cast<std::size_t>(widths[wi]))) {
            pass = false;
            if (diag.empty()) {
              diag = base_diag<T>("format grid invariants", vi);
              diag.push_back(qd_oracle::Tap::Diagnostic("value_limbs",
                                                        qd_oracle::limbs_hex(values[vi])));
              diag.push_back(qd_oracle::Tap::Diagnostic("text", text));
              diag.push_back(qd_oracle::Tap::Diagnostic("precision",
                                                        int_text(precisions[pi])));
              diag.push_back(qd_oracle::Tap::Diagnostic("width",
                                                        int_text(widths[wi])));
            }
          }
        }
      }
    }
  }

  std::string name = std::string(traits::name()) + " format grid invariants";
  tap.ok(pass, name, diag);
  if (verbose || pass) {
    std::cout << "# " << traits::name()
              << " format grid combinations=72\n";
  }
  return pass;
}

template <class T>
bool run_special_format(qd_oracle::Tap &tap) {
  typedef qd_oracle::TypeTraits<T> traits;
  bool pass = true;
  std::string zero = T(0).to_string(3, 0, std::ios_base::scientific);
  std::string inf = T::_inf.to_string(3, 0, std::ios_base::scientific);
  std::string nan = T::_nan.to_string(3, 0, std::ios_base::scientific);
  pass &= zero.find("0.000") != std::string::npos;
  pass &= inf.find("inf") != std::string::npos;
  pass &= nan.find("nan") != std::string::npos;

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) {
    diag = base_diag<T>("special format strings", 0);
    diag.push_back(qd_oracle::Tap::Diagnostic("zero", zero));
    diag.push_back(qd_oracle::Tap::Diagnostic("inf", inf));
    diag.push_back(qd_oracle::Tap::Diagnostic("nan", nan));
  }

  std::string name = std::string(traits::name()) + " special format strings";
  tap.ok(pass, name, diag);
  return pass;
}

template <class T>
bool run_special_parse(qd_oracle::Tap &) {
  return true;
}

template <>
bool run_special_parse<td_real>(qd_oracle::Tap &tap) {
  bool pass = true;
  pass &= isinf(td_real("inf"));
  pass &= isinf(td_real("-inf"));
  pass &= isnan(td_real("nan"));
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) diag = base_diag<td_real>("special parse strings", 0);
  tap.ok(pass, "td special parse strings", diag);
  return pass;
}

#ifdef QD_HAVE_EDD_REAL
template <>
bool run_special_parse<edd_real>(qd_oracle::Tap &tap) {
  bool pass = true;
  pass &= isinf(edd_real("inf"));
  pass &= isinf(edd_real("-inf"));
  pass &= isnan(edd_real("nan"));
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass) diag = base_diag<edd_real>("special parse strings", 0);
  tap.ok(pass, "edd special parse strings", diag);
  return pass;
}
#endif

template <class T>
bool run_type(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  pass &= run_scientific_roundtrip<T>(tap, verbose);
  pass &= run_stream_roundtrip<T>(tap, verbose);
  pass &= run_write_roundtrip<T>(tap, verbose);
  pass &= run_reduced_precision<T>(tap, verbose);
  pass &= run_format_grid<T>(tap, verbose);
  pass &= run_special_format<T>(tap);
  pass &= run_special_parse<T>(tap);
  return pass;
}

int type_plan_dd_qd() { return 6; }
int type_plan_td_edd() { return 7; }

int selected_plan(const Options &options) {
  int count = 0;
  if (options.test_dd) count += type_plan_dd_qd();
  if (options.test_td) count += type_plan_td_edd();
  if (options.test_qd) count += type_plan_dd_qd();
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) count += type_plan_td_edd();
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
