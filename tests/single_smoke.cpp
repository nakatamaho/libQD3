/*
 * tests/single_smoke.cpp
 *
 * Smoke coverage for the float expansion public types.
 */

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

static_assert(sizeof(ds_real) == 2 * sizeof(float),
              "ds_real must contain two float limbs");
static_assert(sizeof(ts_real) == 3 * sizeof(float),
              "ts_real must contain three float limbs");
static_assert(sizeof(qs_real) == 4 * sizeof(float),
              "qs_real must contain four float limbs");

struct TestContext {
  bool pass;

  TestContext() : pass(true) {}

  void check(const char *name, bool condition) {
    if (!condition) {
      std::cout << "FAIL single_smoke: " << name << std::endl;
      pass = false;
    }
  }
};

template <class T>
bool close_to(const T &got, long double reference, float factor = 8.0f) {
  const long double error = std::fabs(to_long_double(got) - reference);
  const long double scale = std::max(1.0L, std::fabs(reference));
  const long double reference_epsilon = std::max(
      static_cast<long double>(T::_eps),
      4096.0L * std::numeric_limits<long double>::epsilon());
  return std::isfinite(error) &&
      error <= static_cast<long double>(factor) * reference_epsilon * scale;
}

template <class T>
void check_type(TestContext &test, const char *name) {
  const T a("1.234567890123456789");
  const T b("-0.345678901234567890");
  const long double da = to_long_double(a);
  const long double db = to_long_double(b);
  const std::string prefix = std::string(name) + " ";

  test.check((prefix + "finite construction").c_str(), a.isfinite());
  T limbs;
  for (int i = 0; i < static_cast<int>(sizeof(limbs.x) / sizeof(limbs.x[0]));
       ++i) {
    limbs[i] = static_cast<float>(i + 1);
  }
  for (int i = 0; i < static_cast<int>(sizeof(limbs.x) / sizeof(limbs.x[0]));
       ++i) {
    test.check((prefix + "float limb storage").c_str(),
               limbs[i] == static_cast<float>(i + 1));
  }
  test.check((prefix + "addition").c_str(), close_to(a + b, da + db));
  test.check((prefix + "subtraction").c_str(), close_to(a - b, da - db));
  test.check((prefix + "multiplication").c_str(), close_to(a * b, da * db));
  test.check((prefix + "division").c_str(), close_to(a / b, da / db));
  test.check((prefix + "sqr").c_str(), close_to(sqr(a), da * da));
  test.check((prefix + "sqrt").c_str(), close_to(sqrt(a), std::sqrt(da)));
  test.check((prefix + "pow integer").c_str(),
             close_to(pow(a, 3), std::pow(da, 3.0L)));
  test.check((prefix + "pow expansion").c_str(),
             close_to(pow(a, T("1.5")), std::pow(da, 1.5L)));
  test.check((prefix + "exp/log round trip").c_str(),
             close_to(exp(log(a)), da, 32.0f));
  test.check((prefix + "sin/cos identity").c_str(),
             close_to(sqr(sin(a)) + sqr(cos(a)), 1.0L, 64.0f));

  const T shifted = ldexp(a, 7);
  test.check((prefix + "ldexp").c_str(),
             close_to(shifted, std::ldexp(da, 7)));
  test.check((prefix + "absolute value").c_str(), fabs(-a) == fabs(a));
  test.check((prefix + "integer conversion").c_str(), to_int(T("-7.75")) == -7);
  test.check((prefix + "comparison").c_str(), a > b && b < a && a != b);

  const T nan = T::_nan;
  const T inf = T::_inf;
  test.check((prefix + "NaN").c_str(), nan.isnan() && nan != nan);
  test.check((prefix + "infinity").c_str(), inf.isinf() && inf > a);
  test.check((prefix + "zero division").c_str(), (T(1.0f) / T(0.0f)).isinf());
  test.check((prefix + "negative sqrt").c_str(), sqrt(T(-1.0f)).isnan());

  char buffer[128];
  a.write(buffer, sizeof(buffer));
  test.check((prefix + "write/read round trip").c_str(),
             close_to(T(buffer), da, 32.0f));
  test.check((prefix + "stream round trip").c_str(), [&]() {
    std::ostringstream output;
    output << std::setprecision(T::_ndigits) << a;
    std::istringstream input(output.str());
    T round_trip;
    input >> round_trip;
    return close_to(round_trip, da, 32.0f);
  }());

  test.check((prefix + "epsilon").c_str(),
             std::numeric_limits<T>::epsilon() == T::_eps);
  test.check((prefix + "digits").c_str(),
             std::numeric_limits<T>::digits > 0);
}

} // namespace

int main() {
  unsigned int old_cw = 0;
  fpu_fix_start(&old_cw);

  const bool old_ds_suppress = ds_suppress_error_messages;
  const bool old_ts_suppress = ts_suppress_error_messages;
  const bool old_qs_suppress = qs_suppress_error_messages;
  ds_suppress_error_messages = true;
  ts_suppress_error_messages = true;
  qs_suppress_error_messages = true;

  TestContext test;
  const ds_real pair(1.0f, 2.0f);
  const ts_real triple_limbs(1.0f, 2.0f, 3.0f);
  const qs_real quad_limbs(1.0f, 2.0f, 3.0f, 4.0f);
  test.check("ds_real has two float constructor limbs",
             pair[0] == 1.0f && pair[1] == 2.0f);
  test.check("ts_real has three float constructor limbs",
             triple_limbs[0] == 1.0f && triple_limbs[1] == 2.0f &&
                 triple_limbs[2] == 3.0f);
  test.check("qs_real has four float constructor limbs",
             quad_limbs[0] == 1.0f && quad_limbs[1] == 2.0f &&
                 quad_limbs[2] == 3.0f && quad_limbs[3] == 4.0f);
  check_type<ds_real>(test, "ds_real");
  check_type<ts_real>(test, "ts_real");
  check_type<qs_real>(test, "qs_real");

  const ts_real triple("1.234567890123456789");
  const qs_real quad(triple, 0.125f);
  test.check("qs_real accepts ts_real plus float precision",
             close_to(quad, to_long_double(triple) + 0.125L, 8.0f));
  ds_suppress_error_messages = old_ds_suppress;
  ts_suppress_error_messages = old_ts_suppress;
  qs_suppress_error_messages = old_qs_suppress;
  test.check("ds suppression flag restored",
             ds_suppress_error_messages == old_ds_suppress);
  test.check("ts suppression flag restored",
             ts_suppress_error_messages == old_ts_suppress);
  test.check("qs suppression flag restored",
             qs_suppress_error_messages == old_qs_suppress);

  fpu_fix_end(&old_cw);
  std::cout << (test.pass ? "PASS single_smoke" : "FAIL single_smoke")
            << std::endl;
  return test.pass ? 0 : 1;
}
