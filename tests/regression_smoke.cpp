/*
 * tests/regression_smoke.cpp
 *
 * Dependency-free regression coverage for public APIs which are easy to
 * exercise indirectly without noticing when they regress.
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include <qd/c_dd.h>
#include <qd/c_qd.h>
#include <qd/c_td.h>
#include <qd/dd_real.h>
#include <qd/fpu.h>
#include <qd/qd_real.h>
#include <qd/td_real.h>
#include <qd/ds_real.h>
#include <qd/ts_real.h>
#include <qd/qs_real.h>
#include <qd/qd_random.h>
#include <limits>

namespace {

struct TestContext {
  bool pass;

  TestContext() : pass(true) {}

  void check(const char *name, bool condition) {
    if (!condition) {
      std::cout << "FAIL regression_smoke: " << name << std::endl;
      pass = false;
    }
  }
};

template <class T>
qd_real td_reference(const T &a) {
  return qd_real(a[0]) + qd_real(a[1]) + qd_real(a[2]);
}

double qd_scale(const qd_real &a) {
  return std::max(1.0, std::abs(to_double(a)));
}

bool td_close(const td_real &got, const qd_real &reference, double factor) {
  const double error = to_double(abs(td_reference(got) - reference));
  return std::isfinite(error) &&
         error <= factor * td_real::_eps * qd_scale(reference);
}

template <int N>
bool equal_limbs(const double *got, const double *expected) {
  for (int i = 0; i < N; ++i) {
    if (got[i] != expected[i]) return false;
  }
  return true;
}

void check_c_api(TestContext &test) {
  double dd_2pi[2];
  double qd_2pi[4];
  double td_2pi[3];
  c_dd_2pi(dd_2pi);
  c_qd_2pi(qd_2pi);
  c_td_2pi(td_2pi);

  test.check("c_dd_2pi", equal_limbs<2>(dd_2pi, dd_real::_2pi.x));
  test.check("c_qd_2pi", equal_limbs<4>(qd_2pi, qd_real::_2pi.x));
  test.check("c_td_2pi", equal_limbs<3>(td_2pi, td_real::_2pi.x));
  test.check("c_dd_epsilon", c_dd_epsilon() == dd_real::_eps);
  test.check("c_qd_epsilon", c_qd_epsilon() == qd_real::_eps);
  test.check("c_td_epsilon", c_td_epsilon() == td_real::_eps);

  // This uses the rescaling path and verifies that the C return status is
  // initialized and propagated from fsqrt.
  const double input[4] = {std::ldexp(1.0, 1000), 0.0, 0.0, 0.0};
  double output[4] = {0.0, 0.0, 0.0, 0.0};
  const int status = c_qd_sqrt(input, output);
  const qd_real got(output);
  const qd_real expected(std::ldexp(1.0, 500));
  test.check("c_qd_sqrt status", status == 0);
  test.check("c_qd_sqrt rescaled result", got == expected);
}

class CerrCapture {
public:
  CerrCapture() : old_(std::cerr.rdbuf(buffer_.rdbuf())) {}

  ~CerrCapture() { std::cerr.rdbuf(old_); }

  void clear() {
    buffer_.str(std::string());
    buffer_.clear();
  }

  std::string str() const { return buffer_.str(); }

private:
  // buffer_ must be constructed before old_'s initializer uses it.
  std::ostringstream buffer_;
  std::streambuf *old_;
};

bool contains(const std::string &text, const char *needle) {
  return text.find(needle) != std::string::npos;
}

void check_error_suppression(TestContext &test) {
  const bool old_dd = dd_suppress_error_messages;
  const bool old_td = td_suppress_error_messages;
  const bool old_qd = qd_suppress_error_messages;
  CerrCapture capture;

  const auto invoke_all = []() {
    (void)::sqrt(dd_real(-1.0));
    (void)::sqrt(td_real(-1.0));
    (void)::sqrt(qd_real(-1.0));
  };

  dd_suppress_error_messages = false;
  td_suppress_error_messages = false;
  qd_suppress_error_messages = false;
  invoke_all();
  std::string output = capture.str();
  test.check("DD error message is emitted", contains(output, "(dd_real::sqrt)"));
  test.check("TD error message is emitted", contains(output, "(td_real::sqrt)"));
  test.check("QD error message is emitted", contains(output, "(qd_real::sqrt)"));

  capture.clear();
  dd_suppress_error_messages = true;
  td_suppress_error_messages = false;
  qd_suppress_error_messages = false;
  invoke_all();
  output = capture.str();
  test.check("DD suppression is independent", !contains(output, "(dd_real::sqrt)"));
  test.check("TD remains unsuppressed", contains(output, "(td_real::sqrt)"));
  test.check("QD remains unsuppressed", contains(output, "(qd_real::sqrt)"));

  capture.clear();
  dd_suppress_error_messages = false;
  td_suppress_error_messages = true;
  qd_suppress_error_messages = false;
  invoke_all();
  output = capture.str();
  test.check("TD suppression is independent", !contains(output, "(td_real::sqrt)"));
  test.check("DD remains unsuppressed", contains(output, "(dd_real::sqrt)"));
  test.check("QD remains unsuppressed", contains(output, "(qd_real::sqrt)"));

  capture.clear();
  dd_suppress_error_messages = false;
  td_suppress_error_messages = false;
  qd_suppress_error_messages = true;
  invoke_all();
  output = capture.str();
  test.check("QD suppression is independent", !contains(output, "(qd_real::sqrt)"));
  test.check("DD remains unsuppressed after QD toggle",
             contains(output, "(dd_real::sqrt)"));
  test.check("TD remains unsuppressed after QD toggle",
             contains(output, "(td_real::sqrt)"));

  capture.clear();
  dd_suppress_error_messages = true;
  td_suppress_error_messages = true;
  qd_suppress_error_messages = true;
  invoke_all();
  test.check("all error messages can be suppressed", capture.str().empty());

  dd_suppress_error_messages = old_dd;
  td_suppress_error_messages = old_td;
  qd_suppress_error_messages = old_qd;
}

template <class T>
void check_special_comparisons(TestContext &test, const char *name) {
  const T nan = T::_nan;
  const T inf = T::_inf;
  const T neg_inf = -T::_inf;
  const T one(1.0);

  const std::string prefix = std::string(name) + " ";
  const std::string nan_label = prefix + "NaN comparisons";
  test.check((nan_label + " is NaN").c_str(), nan.isnan());
  test.check((nan_label + " != itself").c_str(), nan != nan);
  test.check((nan_label + " is not equal").c_str(), !(nan == nan));
  test.check((nan_label + " is unordered below").c_str(), !(nan < one));
  test.check((nan_label + " is unordered above").c_str(), !(nan > one));
  test.check((nan_label + " is unordered at most").c_str(), !(nan <= one));
  test.check((nan_label + " is unordered at least").c_str(), !(nan >= one));

  test.check((prefix + "positive infinity equality").c_str(), inf == inf);
  test.check((prefix + "positive infinity <= itself").c_str(), inf <= inf);
  test.check((prefix + "positive infinity >= itself").c_str(), inf >= inf);
  test.check((prefix + "positive infinity is above one").c_str(), inf > one);
  test.check((prefix + "negative infinity is below one").c_str(), neg_inf < one);
  test.check((prefix + "negative infinity is not above one").c_str(),
             !(neg_inf > one));
}

void check_td_regressions(TestContext &test) {
  const td_real a(
      "1.0000000000000002220446049250313080847263336181640625");
  const td_real b(
      "1.00000000000000011102230246251565404236316680908203125");
  const qd_real add_reference = qd_real(a[0]) + qd_real(a[1]) + qd_real(a[2]) +
                                qd_real(b[0]) + qd_real(b[1]) + qd_real(b[2]);
  test.check("TD ieee_add", td_close(td_real::ieee_add(a, b), add_reference, 64.0));
  test.check("TD sloppy_add", td_close(td_real::sloppy_add(a, b), add_reference,
                                       4096.0));

  const char *sa = "1.2345678901234567890123456789012345678901e120";
  const char *sb = "9.8765432109876543210987654321098765432109e30";
  const td_real da(sa);
  const td_real db(sb);
  const qd_real div_reference = qd_real(sa) / qd_real(sb);
  test.check("TD accurate_div", td_close(td_real::accurate_div(da, db),
                                           div_reference, 512.0));
  test.check("TD sloppy_div", td_close(td_real::sloppy_div(da, db),
                                         div_reference, 4096.0));

  const td_real base("1.25");
  test.check("TD pow(int)", td_close(pow(base, 3), qd_real("1.953125"), 64.0));
  test.check("TD pow(real)",
             td_close(pow(base, td_real("2.5")),
                      pow(qd_real("1.25"), qd_real("2.5")), 512.0));

  const td_real value("1.2345678901234567890123456789012345678901");
  const td_real shifted = ldexp(value, 17);
  const td_real expected_shift(std::ldexp(value[0], 17),
                               std::ldexp(value[1], 17),
                               std::ldexp(value[2], 17));
  test.check("TD ldexp", shifted[0] == expected_shift[0] &&
                            shifted[1] == expected_shift[1] &&
                            shifted[2] == expected_shift[2]);
  test.check("TD to_int truncates", to_int(td_real("-7.75")) == -7);
  test.check("TD fabs", fabs(-value) == fabs(value));

  char buffer[128];
  value.write(buffer, sizeof(buffer));
  test.check("TD write round trip", td_close(td_real(buffer),
                                               td_reference(value), 64.0));
}

// Decimal literals that are exactly representable must parse exactly; the
// readers used to multiply by an inexact 10^-k, so td_real("3.0") was off.
template <class T>
void check_exact_decimal_parse(TestContext &test, const char *name) {
  static const char *const literals[] = {"0.5", "0.25", "1.5", "2.75", "3.0",
                                         "7.5", "0.0625", "-0.5", "123.375",
                                         "1e1", "5e-1", "1024"};
  for (const char *literal : literals) {
    const std::string label = std::string(name) + " parses " + literal + " exactly";
    test.check(label.c_str(), T(literal) == T(std::strtod(literal, 0)));
  }
}

// Long decimal literals must round to within one epsilon of the
// quad-double value for the binary32-based types.
template <class T>
void check_long_decimal_parse(TestContext &test, const char *name) {
  static const char *const literals[] = {
      "1.414213562373095048801688724209698078569671875376948073176679737990732e0",
      "1.259921049894873164767210607278228350570251464701507980081975112155299e0",
      "3.333333333333333333333333333333333333333333333333333333333333333333333e-1",
      "-2.718281828459045235360287471352662497757247093699959574966967627724076e-7"};
  const double eps = static_cast<double>(std::numeric_limits<T>::epsilon());
  for (const char *literal : literals) {
    const T got(literal);
    const qd_real reference(literal);
    qd_real sum(0.0);
    for (double limb : got.x) sum += limb;
    const double error = to_double(abs((sum - reference) / reference));
    const std::string label = std::string(name) + " long literal within 1 eps";
    test.check(label.c_str(), error <= eps);
  }
}

// Random functions must cover [0, 1) uniformly with every limb populated,
// independently of the platform's RAND_MAX, and be reproducible by seed.
template <class T, class Draw>
void check_random(TestContext &test, const char *name, Draw draw, int limbs) {
  const int samples = 20000;
  double lo = 1.0, hi = 0.0, sum = 0.0;
  int last_limb_set = 0;
  bool in_range = true;
  for (int i = 0; i < samples; ++i) {
    const T v = draw();
    const double d = to_double(v);
    in_range = in_range && d >= 0.0 && v < T(1.0);
    lo = std::min(lo, d);
    hi = std::max(hi, d);
    sum += d;
    if (static_cast<double>(v.x[limbs - 1]) != 0.0) ++last_limb_set;
  }
  const std::string prefix = std::string(name) + " random ";
  test.check((prefix + "in [0, 1)").c_str(), in_range);
  test.check((prefix + "covers the interval").c_str(), lo < 1e-3 && hi > 0.999);
  test.check((prefix + "mean").c_str(), std::fabs(sum / samples - 0.5) < 0.01);
  test.check((prefix + "fills the last limb").c_str(), last_limb_set > samples * 99 / 100);
  qd_srand(12345);
  const T first = draw();
  qd_srand(12345);
  test.check((prefix + "reproducible by seed").c_str(), draw() == first);
}

} // namespace

int main() {
  unsigned int old_cw = 0;
  fpu_fix_start(&old_cw);

  TestContext test;
  check_c_api(test);
  check_error_suppression(test);
  check_special_comparisons<dd_real>(test, "dd_real");
  check_special_comparisons<td_real>(test, "td_real");
  check_special_comparisons<qd_real>(test, "qd_real");
  check_td_regressions(test);
  check_exact_decimal_parse<dd_real>(test, "dd_real");
  check_exact_decimal_parse<td_real>(test, "td_real");
  check_exact_decimal_parse<qd_real>(test, "qd_real");
  check_exact_decimal_parse<ds_real>(test, "ds_real");
  check_exact_decimal_parse<ts_real>(test, "ts_real");
  check_exact_decimal_parse<qs_real>(test, "qs_real");
  check_long_decimal_parse<ds_real>(test, "ds_real");
  check_long_decimal_parse<ts_real>(test, "ts_real");
  check_long_decimal_parse<qs_real>(test, "qs_real");
  check_random<dd_real>(test, "dd_real", [] { return ddrand(); }, 2);
  check_random<td_real>(test, "td_real", [] { return tdrand(); }, 3);
  check_random<qd_real>(test, "qd_real", [] { return qdrand(); }, 4);
  check_random<ds_real>(test, "ds_real", [] { return ds_real::rand(); }, 2);
  check_random<ts_real>(test, "ts_real", [] { return ts_real::rand(); }, 3);
  check_random<qs_real>(test, "qs_real", [] { return qs_real::rand(); }, 4);

  fpu_fix_end(&old_cw);
  std::cout << (test.pass ? "PASS regression_smoke" : "FAIL regression_smoke")
            << std::endl;
  return test.pass ? 0 : 1;
}
