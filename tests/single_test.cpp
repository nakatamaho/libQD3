/* Detailed, backend-independent tests for the float expansion types. */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <qd/ds_real.h>
#include <qd/fpu.h>
#include <qd/qs_real.h>
#include <qd/ts_real.h>

namespace {

struct TestContext {
  bool pass;
  TestContext() : pass(true) {}

  void check(const std::string &name, bool condition) {
    if (!condition) {
      std::cout << "FAIL single_test: " << name << std::endl;
      pass = false;
    }
  }
};

template <class T>
bool close_to(const T &value, long double reference, long double factor = 64.0L) {
  const long double error = std::fabs(to_long_double(value) - reference);
  const long double scale = std::max(1.0L, std::fabs(reference));
  const long double epsilon = std::max(
      static_cast<long double>(T::_eps),
      4096.0L * std::numeric_limits<long double>::epsilon());
  return std::isfinite(error) && error <= factor * epsilon * scale;
}

template <class T>
void test_construction(TestContext &test, const char *name) {
  const std::string prefix(name);
  T limbs;
  for (int i = 0; i < static_cast<int>(sizeof(limbs.x) / sizeof(limbs.x[0]));
       ++i) {
    limbs[i] = static_cast<float>(i + 1);
  }
  for (int i = 0; i < static_cast<int>(sizeof(limbs.x) / sizeof(limbs.x[0]));
       ++i) {
    test.check(prefix + " limb constructor " + std::to_string(i),
               limbs[i] == static_cast<float>(i + 1));
  }

  const T decimal("1.2345678901234567890123456789");
  test.check(prefix + " decimal constructor",
             close_to(decimal, 1.2345678901234567890123456789L));
  test.check(prefix + " integer constructor", T(-17) == T("-17"));
  test.check(prefix + " copy constructor", T(decimal) == decimal);
  test.check(prefix + " cross expansion constructor",
             close_to(T(T("1.234567890123456789"), 0.125f),
                      to_long_double(T("1.234567890123456789")) + 0.125L));
  test.check(prefix + " finite state", decimal.isfinite());
  test.check(prefix + " numeric limits digits",
             std::numeric_limits<T>::digits > 0 &&
                 std::numeric_limits<T>::digits10 > 0);
  test.check(prefix + " numeric limits epsilon",
             std::numeric_limits<T>::epsilon() == T::_eps);
}

template <class T>
void test_arithmetic(TestContext &test, const char *name) {
  const std::string prefix(name);
  const T a("1.234567890123456789");
  const T b("-0.345678901234567890");
  const long double da = to_long_double(a);
  const long double db = to_long_double(b);

  test.check(prefix + "addition", close_to(a + b, da + db));
  test.check(prefix + "subtraction", close_to(a - b, da - db));
  test.check(prefix + "multiplication", close_to(a * b, da * db));
  test.check(prefix + "division", close_to(a / b, da / db));
  test.check(prefix + "negative", close_to(-a, -da));
  test.check(prefix + "self cancellation", (a - a).is_zero());
  const long double f01 = static_cast<long double>(0.1f);
  const long double f03 = static_cast<long double>(3.0f);
  test.check(prefix + "static add", close_to(T::add(1.0f, 0.1f),
                                             1.0L + f01));
  test.check(prefix + "static sub", close_to(T::sub(1.0f, 0.1f),
                                             1.0L - f01));
  test.check(prefix + "static mul", close_to(T::mul(1.0f, 0.1f),
                                             1.0L * f01));
  test.check(prefix + "static div", close_to(T::div(1.0f, 3.0f),
                                             1.0L / f03));
  test.check(prefix + "ieee add", close_to(T::ieee_add(a, b), da + db));
  test.check(prefix + "sloppy add", close_to(T::sloppy_add(a, b), da + db));
  test.check(prefix + "branch-free add", close_to(T::bf_add(a, b), da + db));
  test.check(prefix + "branch-free mul", close_to(T::bf_mul(a, b), da * db));
  test.check(prefix + "accurate div",
             close_to(T::accurate_div(a, b), da / db));
  test.check(prefix + "sloppy div", close_to(T::sloppy_div(a, b), da / db));
  test.check(prefix + "inverse", close_to(inv(a), 1.0L / da));
  test.check(prefix + "power of two scale", close_to(ldexp(a, 11),
                                                        std::ldexp(da, 11)));
}

template <class T>
void test_functions(TestContext &test, const char *name) {
  const std::string prefix(name);
  const T x("1.234567890123456789");
  const long double dx = to_long_double(x);
  const T positive("2.34567890123456789");
  const long double dp = to_long_double(positive);

  test.check(prefix + "sqr", close_to(sqr(x), dx * dx));
  test.check(prefix + "sqrt", close_to(sqrt(positive), std::sqrt(dp)));
  test.check(prefix + "sqrt squared",
             close_to(sqr(sqrt(positive)), dp, 256.0L));
  test.check(prefix + "npwr", close_to(npwr(x, 4), std::pow(dx, 4.0L),
                                          256.0L));
  test.check(prefix + "pow integer",
             close_to(pow(x, -3), std::pow(dx, -3.0L), 512.0L));
  test.check(prefix + "nroot", close_to(nroot(T("27"), 3), 3.0L,
                                          512.0L));
  test.check(prefix + "cbrt", close_to(cbrt(T("-8")), -2.0L, 512.0L));

  const T small("0.03125");
  test.check(prefix + "exp", close_to(exp(small), std::exp(0.03125L),
                                        512.0L));
  test.check(prefix + "log", close_to(log(positive), std::log(dp),
                                        1024.0L));
  test.check(prefix + "expm1", close_to(expm1(small), std::expm1(0.03125L),
                                          512.0L));
  test.check(prefix + "log1p", close_to(log1p(small), std::log1p(0.03125L),
                                          1024.0L));
  test.check(prefix + "exp/log identity", close_to(exp(log(positive)), dp,
                                                   4096.0L));
  test.check(prefix + "log10", close_to(log10(positive), std::log10(dp),
                                         2048.0L));
  test.check(prefix + "log2", close_to(log2(positive), std::log2(dp),
                                        2048.0L));
  test.check(prefix + "exp2", close_to(exp2(small), std::exp2(0.03125L),
                                        2048.0L));
}

template <class T>
void test_trigonometry(TestContext &test, const char *name) {
  const std::string prefix(name);
  const T x("0.375");
  const long double dx = to_long_double(x);
  T sine;
  T cosine;
  sincos(x, sine, cosine);
  test.check(prefix + "sincos sine", close_to(sine, std::sin(dx), 2048.0L));
  test.check(prefix + "sincos cosine",
             close_to(cosine, std::cos(dx), 2048.0L));
  test.check(prefix + "sincos consistency", sine == sin(x) &&
                                                  cosine == cos(x));
  test.check(prefix + "trig identity",
             close_to(sqr(sine) + sqr(cosine), 1.0L, 4096.0L));
  test.check(prefix + "tan", close_to(tan(x), std::tan(dx), 4096.0L));
  test.check(prefix + "asin", close_to(asin(sine), dx, 8192.0L));
  test.check(prefix + "acos", close_to(acos(cosine), dx, 8192.0L));
  test.check(prefix + "atan", close_to(atan(x), std::atan(dx), 4096.0L));
  test.check(prefix + "atan2", close_to(atan2(sine, cosine),
                                           std::atan2(std::sin(dx),
                                                      std::cos(dx)),
                                           8192.0L));

  const T h("0.125");
  const long double dh = to_long_double(h);
  T sinh_value;
  T cosh_value;
  sincosh(h, sinh_value, cosh_value);
  test.check(prefix + "sinh", close_to(sinh_value, std::sinh(dh), 4096.0L));
  test.check(prefix + "cosh", close_to(cosh_value, std::cosh(dh), 4096.0L));
  test.check(prefix + "tanh", close_to(tanh(h), std::tanh(dh), 4096.0L));
  test.check(prefix + "hyperbolic identity",
             close_to(sqr(cosh_value) - sqr(sinh_value), 1.0L, 8192.0L));
  test.check(prefix + "asinh", close_to(asinh(h), std::asinh(dh), 8192.0L));
  test.check(prefix + "acosh", close_to(acosh(T("1.25")), std::acosh(1.25L),
                                           8192.0L));
  test.check(prefix + "atanh", close_to(atanh(h), std::atanh(dh), 8192.0L));
}

template <class T>
void test_rounding_io_and_polynomial(TestContext &test, const char *name) {
  const std::string prefix(name);
  test.check(prefix + "nint positive", nint(T("2.5")) == T("3"));
  test.check(prefix + "nint negative", nint(T("-2.5")) == T("-2"));
  test.check(prefix + "floor", floor(T("-1.25")) == T("-2"));
  test.check(prefix + "ceil", ceil(T("-1.25")) == T("-1"));
  test.check(prefix + "aint", aint(T("-1.25")) == T("-1"));
  test.check(prefix + "trunc", trunc(T("-1.75")) == T("-1"));
  test.check(prefix + "round", round(T("-1.5")) == T("-2"));
  test.check(prefix + "fmod", close_to(fmod(T("7.5"), T("2")), 1.5L));
  test.check(prefix + "rem", close_to(rem(T("7.5"), T("2")), 1.5L));
  test.check(prefix + "drem", close_to(drem(T("7.5"), T("2")), -0.5L));
  T remainder;
  const T quotient = divrem(T("7.5"), T("2"), remainder);
  test.check(prefix + "divrem quotient", quotient == T("3"));
  test.check(prefix + "divrem remainder", remainder == T("1.5"));
  test.check(prefix + "hypot", close_to(hypot(T("3"), T("4")), 5.0L));

  T coefficients[3] = {T("-1"), T("0"), T("1")};
  const T root = polyroot(coefficients, 2, T("0.75"));
  test.check(prefix + "polyeval", close_to(polyeval(coefficients, 2, T("2")),
                                             3.0L));
  test.check(prefix + "polyroot", close_to(root, 1.0L, 8192.0L));

  const T value("-1.234567890123456789");
  char digits[128];
  int exponent = 0;
  value.to_digits(digits, exponent, T::_ndigits);
  test.check(prefix + "to_digits leading digit", digits[0] == '1');
  test.check(prefix + "to_digits exponent", exponent == 0);
  char written[256];
  value.write(written, sizeof(written));
  T read_value;
  T reader;
  test.check(prefix + "read/write", reader.read(written, read_value) == 0 &&
                                          close_to(read_value,
                                                   to_long_double(value), 128.0L));
  std::ostringstream output;
  output << std::setprecision(T::_ndigits) << value;
  std::istringstream input(output.str());
  T stream_value;
  input >> stream_value;
  test.check(prefix + "stream", !input.fail() &&
                                      close_to(stream_value,
                                               to_long_double(value), 128.0L));
  test.check(prefix + "rand finite", T::rand().isfinite());
  test.check(prefix + "debug rand finite", T::debug_rand().isfinite());
}

template <class T>
void test_special(TestContext &test, const char *name) {
  const std::string prefix(name);
  const T nan = T::_nan;
  const T inf = T::_inf;
  test.check(prefix + "nan classification", nan.isnan() && !nan.isfinite());
  test.check(prefix + "nan comparison", nan != nan && !(nan < T(0)));
  test.check(prefix + "inf classification", inf.isinf() && !inf.isfinite());
  test.check(prefix + "infinity arithmetic", (inf + inf).isinf());
  test.check(prefix + "infinity cancellation", (inf - inf).isnan());
  test.check(prefix + "zero division", (T(1) / T(0)).isinf());
  test.check(prefix + "zero sign", std::signbit(to_double(T(-0.0f))));
  test.check(prefix + "negative sqrt", sqrt(T(-1)).isnan());
  test.check(prefix + "log zero", log(T(0)).isinf());
  test.check(prefix + "sqrt four", sqrt(T(4)) == T(2));
}

template <class T>
void run_type(TestContext &test, const char *name) {
  test_construction<T>(test, name);
  test_arithmetic<T>(test, name);
  test_functions<T>(test, name);
  test_trigonometry<T>(test, name);
  test_rounding_io_and_polynomial<T>(test, name);
  test_special<T>(test, name);
}

} // namespace

int main() {
  unsigned int old_cw = 0;
  fpu_fix_start(&old_cw);
  const bool old_ds = ds_suppress_error_messages;
  const bool old_ts = ts_suppress_error_messages;
  const bool old_qs = qs_suppress_error_messages;
  ds_suppress_error_messages = true;
  ts_suppress_error_messages = true;
  qs_suppress_error_messages = true;

  TestContext test;
  run_type<ds_real>(test, "ds_real");
  run_type<ts_real>(test, "ts_real");
  run_type<qs_real>(test, "qs_real");

  ds_suppress_error_messages = old_ds;
  ts_suppress_error_messages = old_ts;
  qs_suppress_error_messages = old_qs;
  fpu_fix_end(&old_cw);
  std::cout << (test.pass ? "PASS single_test" : "FAIL single_test")
            << std::endl;
  return test.pass ? 0 : 1;
}
