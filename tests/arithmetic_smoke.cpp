/*
 * tests/arithmetic_smoke.cpp
 *
 * Focused regression coverage for overflow-safe division and IEEE-style
 * special-value handling across the supported expansion types.
 */

#include <cmath>
#include <iostream>
#include <limits>
#include <qd/dd_real.h>
#include <qd/qd_real.h>
#include <qd/td_real.h>
#ifdef QD_HAVE_EDD_REAL
#include <qd/edd_real.h>
#endif
#include <qd/fpu.h>

using std::cout;
using std::endl;

inline bool signbit0(const dd_real &a) {
  return std::signbit(a._hi());
}

template <class T>
bool signbit0(const T &a) {
  return std::signbit(a[0]);
}

template <class T>
struct division_overflow_denominator {
  using type = double;

  static type value() { return std::numeric_limits<double>::denorm_min(); }
};

#ifdef QD_HAVE_EDD_REAL
template <>
struct division_overflow_denominator<edd_real> {
  using type = edd_word;

  static type value() { return std::numeric_limits<edd_word>::denorm_min(); }
};
#endif

template <class T>
bool check_type(const char *name) {
  const T zero(0.0);
  const T negative_zero(-0.0);
  const T one(1.0);
  const T negative_one(-1.0);
  const T positive_inf = T::_inf;
  const T negative_inf = -T::_inf;
  const T nan = T::_nan;
  bool pass = true;

  const auto check = [&](const char *label, bool condition) {
    if (!condition) {
      cout << "FAIL " << name << ": " << label << endl;
      pass = false;
    }
  };

  check("NaN / one is NaN", (nan / one).isnan());
  check("one / NaN is NaN", (one / nan).isnan());
  check("Inf / one is Inf", (positive_inf / one).isinf());
  check("-Inf / one is negative", (negative_inf / one).is_negative());
  check("one / Inf is positive zero", (one / positive_inf).is_zero() &&
        !signbit0(one / positive_inf));
  check("one / -Inf is negative zero", (one / negative_inf).is_zero() &&
        signbit0(one / negative_inf));
  check("Inf / Inf is NaN", (positive_inf / positive_inf).isnan());
  check("zero / zero is NaN", (zero / zero).isnan());
  check("zero / one is zero", (zero / one).is_zero());
  check("negative zero preserves sign", (negative_zero / one).is_zero() &&
        signbit0(negative_zero / one));
  check("one / zero is Inf", (one / zero).isinf());
  check("-one / zero is negative Inf", (negative_one / zero).isinf() &&
        (negative_one / zero).is_negative());

  check("sqrt(NaN) is NaN", sqrt(nan).isnan());
  check("sqrt(Inf) is Inf", sqrt(positive_inf).isinf());
  check("sqrt(-Inf) is NaN", sqrt(negative_inf).isnan());
  check("sqrt(-one) is NaN", sqrt(negative_one).isnan());
  check("sqrt(zero) is zero", sqrt(zero).is_zero());

  // These exact powers exercise the large-operand rescaling paths.
  const T p1000(std::ldexp(1.0, 1000));
  const T p995(std::ldexp(1.0, 995));
  const T p500(std::ldexp(1.0, 500));
  const T quotient = p1000 / T(32.0);
  const T root = sqrt(p1000);
  check("large division remains finite", quotient.isfinite());
  check("large division preserves exact power of two", quotient == p995);
  check("large sqrt remains finite", root.isfinite() && root.is_positive());
  const double root_error = to_double(abs(root - p500));
  check("large sqrt has the expected value",
        std::isfinite(root_error) &&
        root_error <= 128.0 * std::numeric_limits<double>::epsilon() *
                           to_double(p500));

  const T smallest = T(division_overflow_denominator<T>::value());
  const T direct_overflow = one / smallest;
  check("direct quotient overflow is canonical Inf",
        direct_overflow.isinf() && !direct_overflow.isnan());
  const T negative_direct_overflow = negative_one / smallest;
  check("negative direct quotient overflow is canonical Inf",
        negative_direct_overflow.isinf() &&
        !negative_direct_overflow.isnan() &&
        negative_direct_overflow.is_negative());

  using scalar_type = typename division_overflow_denominator<T>::type;
  const scalar_type scalar_smallest = division_overflow_denominator<T>::value();
  const T direct_scalar_overflow = one / scalar_smallest;
  check("scalar quotient overflow is canonical Inf",
        direct_scalar_overflow.isinf() && !direct_scalar_overflow.isnan());

  const T restored_overflow = T::_max / T(0.25);
  check("rescaled quotient overflow is canonical Inf",
        restored_overflow.isinf() && !restored_overflow.isnan());
  const T negative_restored_overflow = (-T::_max) / T(0.25);
  check("negative rescaled quotient overflow is canonical Inf",
        negative_restored_overflow.isinf() &&
        !negative_restored_overflow.isnan() &&
        negative_restored_overflow.is_negative());

  const T restored_scalar_overflow = T::_max / 0.25;
  check("rescaled scalar quotient overflow is canonical Inf",
        restored_scalar_overflow.isinf() && !restored_scalar_overflow.isnan());

  return pass;
}

int main() {
  bool pass = true;
  const bool old_dd_suppress = dd_suppress_error_messages;
  const bool old_td_suppress = td_suppress_error_messages;
  const bool old_qd_suppress = qd_suppress_error_messages;
  dd_suppress_error_messages = true;
  td_suppress_error_messages = true;
  qd_suppress_error_messages = true;
#ifdef QD_HAVE_EDD_REAL
  const bool old_edd_suppress = edd_suppress_error_messages;
  edd_suppress_error_messages = true;
#endif

  unsigned int old_cw = 0;
  fpu_fix_start(&old_cw);

  pass &= check_type<dd_real>("dd_real");
  pass &= check_type<td_real>("td_real");
  pass &= check_type<qd_real>("qd_real");

  {
    const dd_real dd_smallest(std::numeric_limits<double>::denorm_min());
    const qd_real mixed_direct_overflow = qd_real(1.0) / dd_smallest;
    if (!(mixed_direct_overflow.isinf() && !mixed_direct_overflow.isnan())) {
      cout << "FAIL qd_real: QD/DD quotient overflow is not canonical Inf"
           << endl;
      pass = false;
    }
    const qd_real mixed_restored_overflow = qd_real::_max / dd_real(0.25);
    if (!(mixed_restored_overflow.isinf() &&
          !mixed_restored_overflow.isnan())) {
      cout << "FAIL qd_real: rescaled QD/DD quotient overflow is not canonical Inf"
           << endl;
      pass = false;
    }
  }

#ifdef QD_HAVE_EDD_REAL
  fpu_fix_end(&old_cw);
  unsigned int old_edd_cw = 0;
  fpu_fix_start_80bit(&old_edd_cw);
  pass &= check_type<edd_real>("edd_real");
  fpu_fix_end(&old_edd_cw);
#else
  fpu_fix_end(&old_cw);
#endif

  dd_suppress_error_messages = old_dd_suppress;
  td_suppress_error_messages = old_td_suppress;
  qd_suppress_error_messages = old_qd_suppress;
#ifdef QD_HAVE_EDD_REAL
  edd_suppress_error_messages = old_edd_suppress;
#endif

  cout << (pass ? "PASS arithmetic_smoke" : "FAIL arithmetic_smoke")
       << endl;
  return pass ? 0 : 1;
}
