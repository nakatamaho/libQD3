/*
 * Float expansion arithmetic for ds_real, ts_real, and qs_real.
 *
 * This file is intentionally independent of the legacy double-based real
 * types.  All arithmetic is performed with binary32 limbs and error-free
 * float primitives.
 */
#ifndef _QD_SINGLE_REAL_DIRECT_H
#define _QD_SINGLE_REAL_DIRECT_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>

#include <qd/qd_config.h>
#include <qd/inline.h>

#ifdef isnan
#undef isnan
#endif
#ifdef isfinite
#undef isfinite
#endif
#ifdef isinf
#undef isinf
#endif
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

QD_API extern bool ds_suppress_error_messages;
QD_API extern bool ts_suppress_error_messages;
QD_API extern bool qs_suppress_error_messages;

template <int N>
struct single_real;

template <int N> single_real<N> abs(const single_real<N> &);
template <int N> single_real<N> sqr(const single_real<N> &);
template <int N> single_real<N> sqrt(const single_real<N> &);
template <int N> single_real<N> npwr(const single_real<N> &, int);
template <int N> single_real<N> nroot(const single_real<N> &, int);
template <int N> single_real<N> pow(const single_real<N> &, int);
template <int N>
single_real<N> pow(const single_real<N> &, const single_real<N> &);
template <int N> single_real<N> exp(const single_real<N> &);
template <int N> single_real<N> log(const single_real<N> &);
template <int N> single_real<N> nint(const single_real<N> &);
template <int N>
single_real<N> atan2(const single_real<N> &, const single_real<N> &);
template <int N> int to_int(const single_real<N> &);

namespace qd_single_detail {

template <int N> struct traits;

template <> struct traits<2> {
  static const char *name() { return "ds_real"; }
  enum { digits = 46, digits10 = 13, ndigits = 14 };
};

template <> struct traits<3> {
  static const char *name() { return "ts_real"; }
  enum { digits = 70, digits10 = 21, ndigits = 21 };
};

template <> struct traits<4> {
  static const char *name() { return "qs_real"; }
  enum { digits = 94, digits10 = 28, ndigits = 28 };
};

template <int N> struct error_state;
template <> struct error_state<2> {
  static bool &suppress() { return ds_suppress_error_messages; }
};
template <> struct error_state<3> {
  static bool &suppress() { return ts_suppress_error_messages; }
};
template <> struct error_state<4> {
  static bool &suppress() { return qs_suppress_error_messages; }
};

template <int N>
inline void report_error(const char *message) {
  if (!error_state<N>::suppress() && message != 0) {
    std::string text(message);
    const std::string generic("single_real");
    const std::string::size_type position = text.find(generic);
    if (position != std::string::npos) {
      text.replace(position, generic.size(), traits<N>::name());
    }
    std::cerr << "ERROR " << text << std::endl;
  }
}

inline float quick_two_sum(float a, float b, float &error) {
  const float sum = a + b;
  error = b - (sum - a);
  return sum;
}

inline float two_sum(float a, float b, float &error) {
  const float sum = a + b;
  const float bb = sum - a;
  error = (a - (sum - bb)) + (b - bb);
  return sum;
}

inline float two_diff(float a, float b, float &error) {
  const float difference = a - b;
  const float bb = difference - a;
  error = (a - (difference - bb)) - (b + bb);
  return difference;
}

inline float two_prod(float a, float b, float &error) {
  const float product = a * b;
  if (!std::isfinite(static_cast<double>(product))) {
    error = 0.0f;
    return product;
  }
  error = std::fmaf(a, b, -product);
  return product;
}

inline float nint_float(float value) {
  return value == std::floor(value) ? value : std::floor(value + 0.5f);
}

inline float signed_zero(float value) {
  return std::copysign(0.0f, value);
}

static const int kExpansionCapacity = 64;

template <int Capacity>
inline int grow_expansion(const float *expansion, int length, float value,
                          float *out) {
  float q = value;
  int out_length = 0;
  for (int i = 0; i < length; ++i) {
    float sum;
    float error;
    sum = two_sum(q, expansion[i], error);
    if (error != 0.0f) out[out_length++] = error;
    q = sum;
  }
  if (q != 0.0f || out_length == 0) out[out_length++] = q;
  return out_length;
}

template <int N>
inline void normalize_terms(const float *terms, int count, float out[N]) {
  float sorted[kExpansionCapacity];
  float expansion[kExpansionCapacity];
  float next[kExpansionCapacity];
  const int length = std::min(count, kExpansionCapacity);
  for (int i = 0; i < length; ++i) sorted[i] = terms[i];
  std::sort(sorted, sorted + length,
            [](float a, float b) { return std::fabs(a) < std::fabs(b); });

  int expansion_length = 0;
  for (int i = 0; i < length; ++i) {
    if (sorted[i] == 0.0f) continue;
    const int next_length = grow_expansion<kExpansionCapacity>(
        expansion, expansion_length, sorted[i], next);
    for (int j = 0; j < next_length; ++j) expansion[j] = next[j];
    expansion_length = next_length;
  }

  float current[kExpansionCapacity];
  float folded[kExpansionCapacity];
  for (int i = 0; i < expansion_length; ++i) {
    current[i] = expansion[expansion_length - 1 - i];
  }
  for (int pass = 0; pass < 8 && expansion_length > 1; ++pass) {
    float remainder[kExpansionCapacity];
    int remainder_length = 0;
    float sum = current[expansion_length - 1];
    for (int i = expansion_length - 2; i >= 0; --i) {
      float error;
      sum = two_sum(current[i], sum, error);
      if (error != 0.0f) remainder[remainder_length++] = error;
    }
    folded[0] = sum;
    for (int i = remainder_length - 1; i >= 0; --i) {
      folded[remainder_length - i] = remainder[i];
    }
    const int folded_length = remainder_length + 1;
    for (int i = 0; i < folded_length; ++i) current[i] = folded[i];
    expansion_length = folded_length;
  }

  for (int i = 0; i < N; ++i) out[i] = 0.0f;
  // Accumulate the distilled terms from low to high into the requested
  // number of high-to-low limbs.  Any remainder beyond N limbs is the
  // correctly rounded truncation of the expansion.
  for (int i = expansion_length - 1; i >= 0; --i) {
    float remainder = current[i];
    for (int j = 0; j < N && remainder != 0.0f; ++j) {
      float error;
      out[j] = two_sum(out[j], remainder, error);
      remainder = error;
    }
  }
}

template <int N>
inline void component_terms(const float *a, const float *b,
                            float *terms, int &count) {
  count = 0;
  for (int i = N - 1; i >= 0; --i) terms[count++] = a[i];
  for (int i = N - 1; i >= 0; --i) terms[count++] = b[i];
}

template <int N>
inline void product_terms(const float *a, const float *b,
                          float *terms, int &count) {
  count = 0;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      float error;
      const float product = two_prod(a[i], b[j], error);
      terms[count++] = error;
      terms[count++] = product;
    }
  }
}

} // namespace qd_single_detail

template <int N>
struct QD_API single_real {
  static_assert(N >= 2 && N <= 4,
                "single_real supports two to four float limbs");

  float x[N];

  single_real() {
    for (int i = 0; i < N; ++i) x[i] = 0.0f;
  }
  single_real(float value) {
    x[0] = value;
    for (int i = 1; i < N; ++i) x[i] = 0.0f;
  }
  single_real(double value);
  single_real(long double value);

  template <class I,
            typename std::enable_if<qd::is_supported_integral<I>::value,
                                    int>::type = 0>
  single_real(I value);

  template <class... Rest>
  single_real(float first, Rest... rest) {
    static_assert(sizeof...(Rest) + 1 <= N,
                  "too many limbs for this single expansion");
    const float values[] = {first, static_cast<float>(rest)...};
    const int count = static_cast<int>(sizeof...(Rest)) + 1;
    for (int i = 0; i < count; ++i) x[i] = values[i];
    for (int i = count; i < N; ++i) x[i] = 0.0f;
  }

  explicit single_real(const float *values) {
    for (int i = 0; i < N; ++i) x[i] = values[i];
  }
  single_real(const char *value);

  template <int M> single_real(const single_real<M> &value);
  template <int M> single_real(const single_real<M> &value, float tail);

  float operator[](int i) const { return x[i]; }
  float &operator[](int i) { return x[i]; }
  float _hi() const { return x[0]; }
  float _lo() const { return x[1]; }

  static void error(const char *message) {
    qd_single_detail::report_error<N>(message);
  }

  static const single_real _nan;
  static const single_real _inf;
  static const single_real _2pi;
  static const single_real _pi;
  static const single_real _3pi4;
  static const single_real _pi2;
  static const single_real _pi4;
  static const single_real _e;
  static const single_real _log2;
  static const single_real _log10;
  static const single_real _max;
  static const single_real _safe_max;
  static const float _eps;
  static const float _min_normalized;
  static const int _ndigits;

  bool isnan() const {
    for (int i = 0; i < N; ++i) {
      if (std::isnan(static_cast<double>(x[i]))) return true;
    }
    return false;
  }
  bool isfinite() const { return std::isfinite(static_cast<double>(x[0])); }
  bool isinf() const { return std::isinf(static_cast<double>(x[0])); }

  static single_real add(float, float);
  static single_real ieee_add(const single_real &, const single_real &);
  static single_real sloppy_add(const single_real &, const single_real &);
  static single_real bf_add(const single_real &, const single_real &);
  single_real &operator+=(float);
  single_real &operator+=(const single_real &);
  static single_real sub(float, float);
  single_real &operator-=(float);
  single_real &operator-=(const single_real &);
  single_real operator-() const;
  static single_real mul(float, float);
  static single_real bf_mul(const single_real &, const single_real &);
  single_real &operator*=(float);
  single_real &operator*=(const single_real &);
  static single_real div(float, float);
  static single_real sloppy_div(const single_real &, const single_real &);
  static single_real accurate_div(const single_real &, const single_real &);
  single_real &operator/=(float);
  single_real &operator/=(const single_real &);

  single_real &operator=(float);
  single_real &operator=(double);
  single_real &operator=(long double);
  template <class I,
            typename std::enable_if<qd::is_supported_integral<I>::value,
                                    int>::type = 0>
  single_real &operator=(I);
  single_real &operator=(const char *);
  template <int M> single_real &operator=(const single_real<M> &);

  template <typename Integer,
            typename std::enable_if<qd::is_supported_integral<Integer>::value,
                                    int>::type = 0>
  operator Integer() const {
    return static_cast<Integer>(x[0]);
  }

  single_real operator^(int) const;
  static single_real sqr(float);
  static single_real sqrt(float);
  bool is_zero() const {
    for (int i = 0; i < N; ++i) if (x[i] != 0.0f) return false;
    return true;
  }
  bool is_one() const {
    if (x[0] != 1.0f) return false;
    for (int i = 1; i < N; ++i) if (x[i] != 0.0f) return false;
    return true;
  }
  bool is_positive() const { return x[0] > 0.0f; }
  bool is_negative() const { return x[0] < 0.0f; }
  static single_real rand();
  static single_real debug_rand();

  void to_digits(char *, int &, int precision = _ndigits) const;
  void write(char *, int, int precision = _ndigits,
             bool showpos = false, bool uppercase = false) const;
  std::string to_string(
      int precision = _ndigits, int width = 0,
      std::ios_base::fmtflags fmt =
          static_cast<std::ios_base::fmtflags>(0),
      bool showpos = false, bool uppercase = false,
      char fill = ' ') const;
  int read(const char *, single_real &);
  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
  void dump_bits(const std::string &name = "",
                 std::ostream &os = std::cerr) const;
};

// The template implementation follows the complete class declaration.

namespace qd_single_detail {

template <int N>
inline void component_copy(const single_real<N> &value, float *out) {
  for (int i = 0; i < N; ++i) out[i] = value.x[i];
}

template <int N>
inline single_real<N> from_terms(const float *terms, int count) {
  for (int i = 0; i < count; ++i) {
    if (std::isnan(static_cast<double>(terms[i]))) {
      single_real<N> result;
      const float nan = std::numeric_limits<float>::quiet_NaN();
      for (int j = 0; j < N; ++j) result.x[j] = nan;
      return result;
    }
  }
  for (int i = 0; i < count; ++i) {
    if (std::isinf(static_cast<double>(terms[i]))) {
      single_real<N> result;
      const float inf = std::signbit(terms[i])
          ? -std::numeric_limits<float>::infinity()
          : std::numeric_limits<float>::infinity();
      for (int j = 0; j < N; ++j) result.x[j] = inf;
      return result;
    }
  }
  single_real<N> result;
  qd_single_detail::normalize_terms<N>(terms, count, result.x);
  return result;
}

template <int N>
inline single_real<N> from_long_double(long double value) {
  if (std::isnan(value)) {
    return single_real<N>(std::numeric_limits<float>::quiet_NaN());
  }
  if (std::isinf(value)) {
    return single_real<N>(value < 0.0L
                              ? -std::numeric_limits<float>::infinity()
                              : std::numeric_limits<float>::infinity());
  }
  float terms[N];
  long double remainder = value;
  for (int i = 0; i < N; ++i) {
    terms[i] = static_cast<float>(remainder);
    remainder -= static_cast<long double>(terms[i]);
  }
  return from_terms<N>(terms, N);
}

template <int N>
inline long double to_long_double(const single_real<N> &value) {
  long double result = 0.0L;
  for (int i = N - 1; i >= 0; --i) {
    result += static_cast<long double>(value.x[i]);
  }
  if (result == 0.0L) {
    return std::copysign(0.0L, value.x[0]);
  }
  return result;
}

template <int N>
inline single_real<N> canonical(const single_real<N> &value) {
  float terms[N];
  component_copy(value, terms);
  return from_terms<N>(terms, N);
}

template <int N>
inline single_real<N> scale_components(const single_real<N> &value,
                                       int exponent) {
  single_real<N> result;
  for (int i = 0; i < N; ++i) {
    result.x[i] = std::scalbn(value.x[i], exponent);
  }
  return result;
}

template <int N>
inline single_real<N> parse_decimal(const char *text) {
  if (text == 0) return single_real<N>::_nan;
  std::string input(text);
  const std::string::size_type first_non_space =
      input.find_first_not_of(" \t\n\r\f\v");
  if (first_non_space == std::string::npos) return single_real<N>::_nan;
  const std::string::size_type last_non_space =
      input.find_last_not_of(" \t\n\r\f\v");
  input = input.substr(first_non_space, last_non_space - first_non_space + 1);

  std::string lower = input;
  std::transform(lower.begin(), lower.end(), lower.begin(), [](char c) {
    return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  if (lower == "nan" || lower == "+nan" || lower == "-nan") {
    return single_real<N>(std::numeric_limits<float>::quiet_NaN());
  }
  if (lower == "inf" || lower == "+inf" || lower == "infinity" ||
      lower == "+infinity") {
    return single_real<N>(std::numeric_limits<float>::infinity());
  }
  if (lower == "-inf" || lower == "-infinity") {
    return single_real<N>(-std::numeric_limits<float>::infinity());
  }

  std::string::size_type pos = 0;
  bool negative = false;
  if (input[pos] == '+' || input[pos] == '-') {
    negative = input[pos] == '-';
    ++pos;
  }
  std::string digits;
  int integer_digits = 0;
  bool after_decimal = false;
  bool any_digit = false;
  for (; pos < input.size(); ++pos) {
    const char c = input[pos];
    if (c >= '0' && c <= '9') {
      digits.push_back(c);
      if (!after_decimal) ++integer_digits;
      any_digit = true;
    } else if (c == '.' && !after_decimal) {
      after_decimal = true;
    } else {
      break;
    }
  }
  if (!any_digit) return single_real<N>::_nan;

  int exponent = 0;
  if (pos < input.size() && (input[pos] == 'e' || input[pos] == 'E')) {
    ++pos;
    bool exponent_negative = false;
    if (pos < input.size() && (input[pos] == '+' || input[pos] == '-')) {
      exponent_negative = input[pos] == '-';
      ++pos;
    }
    bool exponent_digit = false;
    while (pos < input.size() && input[pos] >= '0' && input[pos] <= '9') {
      exponent_digit = true;
      if (exponent < 10000) exponent = exponent * 10 + input[pos] - '0';
      ++pos;
    }
    if (!exponent_digit) return single_real<N>::_nan;
    if (exponent_negative) exponent = -exponent;
  }
  if (input.find_first_not_of(" \t\n\r\f\v", pos) != std::string::npos) {
    return single_real<N>::_nan;
  }

  const std::string::size_type first_digit = digits.find_first_not_of('0');
  if (first_digit == std::string::npos) {
    return single_real<N>(negative ? -0.0f : 0.0f);
  }
  const int total_significant =
      static_cast<int>(digits.size() - first_digit);
  const int significant = std::min(total_significant,
                                   static_cast<int>(traits<N>::ndigits));
  // Build an integer prefix using exact binary powers of ten, then scale it
  // by the decimal position.  Repeatedly multiplying a fractional expansion
  // by the float literal 0.1f would permanently retain the binary32
  // approximation of one tenth and lose the decimal digits being parsed.
  single_real<N> mantissa(0.0f);
  for (int i = 0; i < significant; ++i) {
    mantissa *= 10.0f;
    mantissa += static_cast<float>(digits[first_digit + i] - '0');
  }
  const int fractional_digits = after_decimal
      ? static_cast<int>(digits.size()) - integer_digits : 0;
  int decimal_scale = fractional_digits - (total_significant - significant);
  decimal_scale -= exponent;
  if (decimal_scale > 10000) {
    return negative ? -single_real<N>::_inf : single_real<N>::_inf;
  }
  if (decimal_scale < -10000) {
    return single_real<N>(negative ? -0.0f : 0.0f);
  }
  if (decimal_scale > 0) {
    mantissa /= npwr(single_real<N>(10.0f), decimal_scale);
  } else {
    mantissa *= npwr(single_real<N>(10.0f), -decimal_scale);
  }
  return negative ? -mantissa : mantissa;
}

template <int N>
inline single_real<N> binary_add(const single_real<N> &a,
                                 const single_real<N> &b) {
  if (a.isnan() || b.isnan()) return single_real<N>::_nan;
  if (a.isinf() || b.isinf()) {
    if (a.isinf() && b.isinf() && a.is_negative() != b.is_negative()) {
      return single_real<N>::_nan;
    }
    return a.isinf() ? a : b;
  }
  float terms[2 * N];
  int count;
  component_terms<N>(a.x, b.x, terms, count);
  return from_terms<N>(terms, count);
}

template <int N>
inline single_real<N> binary_mul(const single_real<N> &a,
                                 const single_real<N> &b) {
  if (a.isnan() || b.isnan()) return single_real<N>::_nan;
  if ((a.isinf() && b.is_zero()) || (b.isinf() && a.is_zero())) {
    return single_real<N>::_nan;
  }
  if (a.isinf() || b.isinf()) {
    return a.is_negative() != b.is_negative()
        ? -single_real<N>::_inf : single_real<N>::_inf;
  }
  int scale = 0;
  if (!a.is_zero() && !b.is_zero()) {
    const int ea = std::ilogb(std::fabs(a.x[0]));
    const int eb = std::ilogb(std::fabs(b.x[0]));
    if (ea + eb > 100) scale = -32;
    if (ea + eb < -100) scale = 32;
  }
  const single_real<N> aa = scale_components(a, scale);
  const single_real<N> bb = scale_components(b, scale);
  float terms[2 * N * N];
  int count;
  product_terms<N>(aa.x, bb.x, terms, count);
  const single_real<N> result = from_terms<N>(terms, count);
  return scale == 0 ? result : scale_components(result, -2 * scale);
}

template <int N>
inline single_real<N> division_special(const single_real<N> &a,
                                       const single_real<N> &b,
                                       bool &handled) {
  const bool negative = std::signbit(a.x[0]) != std::signbit(b.x[0]);
  if (a.isnan() || b.isnan() || (a.isinf() && b.isinf()) ||
      (a.is_zero() && b.is_zero())) {
    handled = true;
    return single_real<N>::_nan;
  }
  if (b.is_zero() || a.isinf()) {
    handled = true;
    return negative ? -single_real<N>::_inf : single_real<N>::_inf;
  }
  if (a.is_zero() || b.isinf()) {
    handled = true;
    return single_real<N>(std::copysign(0.0f, negative ? -1.0f : 1.0f));
  }
  handled = false;
  return single_real<N>();
}

template <int N>
inline single_real<N> division_impl(const single_real<N> &a,
                                     const single_real<N> &b, int digits) {
  bool handled = false;
  const single_real<N> special = division_special(a, b, handled);
  if (handled) return special;
  int scale = 0;
  const int ea = std::ilogb(std::fabs(a.x[0]));
  const int eb = std::ilogb(std::fabs(b.x[0]));
  if (std::max(ea, eb) > 100) scale = -32;
  if (std::min(ea, eb) < -100) scale = 32;
  const single_real<N> aa = scale_components(a, scale);
  const single_real<N> bb = scale_components(b, scale);
  float quotient[5] = {0, 0, 0, 0, 0};
  single_real<N> remainder = aa;
  const int count = std::min(digits, 5);
  for (int i = 0; i < count; ++i) {
    quotient[i] = remainder.x[0] / bb.x[0];
    if (!std::isfinite(static_cast<double>(quotient[i]))) {
      return quotient[i] < 0.0f ? -single_real<N>::_inf
                                : single_real<N>::_inf;
    }
    remainder -= bb * quotient[i];
  }
  return from_terms<N>(quotient, count);
}

template <int N>
inline single_real<N> signed_zero_result(const single_real<N> &value) {
  return single_real<N>(std::copysign(0.0f, value.x[0]));
}

} // namespace qd_single_detail

template <int N>
inline single_real<N>::single_real(double value)
    : single_real(static_cast<long double>(value)) {}

template <int N>
inline single_real<N>::single_real(long double value) {
  *this = qd_single_detail::from_long_double<N>(value);
}

template <int N>
template <class I,
          typename std::enable_if<qd::is_supported_integral<I>::value,
                                  int>::type>
inline single_real<N>::single_real(I value) {
  *this = qd_single_detail::from_long_double<N>(static_cast<long double>(value));
}

template <int N>
inline single_real<N>::single_real(const char *value) {
  *this = qd_single_detail::parse_decimal<N>(value);
  if (isnan()) error("(single_real::single_real): INPUT ERROR.");
}

template <int N>
template <int M>
inline single_real<N>::single_real(const single_real<M> &value) {
  float terms[M];
  qd_single_detail::component_copy(value, terms);
  *this = qd_single_detail::from_terms<N>(terms, M);
}

template <int N>
template <int M>
inline single_real<N>::single_real(const single_real<M> &value, float tail) {
  float terms[M + 1];
  qd_single_detail::component_copy(value, terms);
  terms[M] = tail;
  *this = qd_single_detail::from_terms<N>(terms, M + 1);
}

template <int N>
inline single_real<N> operator+(const single_real<N> &a,
                                const single_real<N> &b) {
  return qd_single_detail::binary_add(a, b);
}
template <int N>
inline single_real<N> operator+(const single_real<N> &a, float b) {
  return a + single_real<N>(b);
}
template <int N>
inline single_real<N> operator+(float a, const single_real<N> &b) {
  return single_real<N>(a) + b;
}
template <int N>
inline single_real<N> operator+(const single_real<N> &a, double b) {
  return a + single_real<N>(b);
}
template <int N>
inline single_real<N> operator+(double a, const single_real<N> &b) {
  return single_real<N>(a) + b;
}
template <int N>
inline single_real<N> operator-(const single_real<N> &a,
                                const single_real<N> &b) {
  return a + (-b);
}
template <int N>
inline single_real<N> operator-(const single_real<N> &a, float b) {
  return a - single_real<N>(b);
}
template <int N>
inline single_real<N> operator-(float a, const single_real<N> &b) {
  return single_real<N>(a) - b;
}
template <int N>
inline single_real<N> operator-(const single_real<N> &a, double b) {
  return a - single_real<N>(b);
}
template <int N>
inline single_real<N> operator-(double a, const single_real<N> &b) {
  return single_real<N>(a) - b;
}
template <int N>
inline single_real<N> operator*(const single_real<N> &a,
                                const single_real<N> &b) {
  return qd_single_detail::binary_mul(a, b);
}
template <int N>
inline single_real<N> operator*(const single_real<N> &a, float b) {
  return a * single_real<N>(b);
}
template <int N>
inline single_real<N> operator*(float a, const single_real<N> &b) {
  return single_real<N>(a) * b;
}
template <int N>
inline single_real<N> operator*(const single_real<N> &a, double b) {
  return a * single_real<N>(b);
}
template <int N>
inline single_real<N> operator*(double a, const single_real<N> &b) {
  return single_real<N>(a) * b;
}
template <int N>
inline single_real<N> operator/(const single_real<N> &a,
                                const single_real<N> &b) {
#ifdef QD_SLOPPY_DIV
  return single_real<N>::sloppy_div(a, b);
#else
  return single_real<N>::accurate_div(a, b);
#endif
}
template <int N>
inline single_real<N> operator/(const single_real<N> &a, float b) {
  return a / single_real<N>(b);
}
template <int N>
inline single_real<N> operator/(float a, const single_real<N> &b) {
  return single_real<N>(a) / b;
}
template <int N>
inline single_real<N> operator/(const single_real<N> &a, double b) {
  return a / single_real<N>(b);
}
template <int N>
inline single_real<N> operator/(double a, const single_real<N> &b) {
  return single_real<N>(a) / b;
}

template <int N>
inline single_real<N> &single_real<N>::operator+=(float value) {
  *this = *this + value;
  return *this;
}
template <int N>
inline single_real<N> &single_real<N>::operator+=(const single_real &value) {
  *this = *this + value;
  return *this;
}
template <int N>
inline single_real<N> &single_real<N>::operator-=(float value) {
  *this = *this - value;
  return *this;
}
template <int N>
inline single_real<N> &single_real<N>::operator-=(const single_real &value) {
  *this = *this - value;
  return *this;
}
template <int N>
inline single_real<N> single_real<N>::operator-() const {
  single_real result;
  for (int i = 0; i < N; ++i) result.x[i] = -x[i];
  return result;
}
template <int N>
inline single_real<N> &single_real<N>::operator*=(float value) {
  *this = *this * value;
  return *this;
}
template <int N>
inline single_real<N> &single_real<N>::operator*=(const single_real &value) {
  *this = *this * value;
  return *this;
}
template <int N>
inline single_real<N> &single_real<N>::operator/=(float value) {
  *this = *this / value;
  return *this;
}
template <int N>
inline single_real<N> &single_real<N>::operator/=(const single_real &value) {
  *this = *this / value;
  return *this;
}
template <int N>
inline single_real<N> &single_real<N>::operator=(float value) {
  x[0] = value;
  for (int i = 1; i < N; ++i) x[i] = 0.0f;
  return *this;
}
template <int N>
inline single_real<N> &single_real<N>::operator=(double value) {
  return *this = static_cast<long double>(value);
}
template <int N>
inline single_real<N> &single_real<N>::operator=(long double value) {
  *this = qd_single_detail::from_long_double<N>(value);
  return *this;
}
template <int N>
template <class I,
          typename std::enable_if<qd::is_supported_integral<I>::value,
                                  int>::type>
inline single_real<N> &single_real<N>::operator=(I value) {
  *this = qd_single_detail::from_long_double<N>(static_cast<long double>(value));
  return *this;
}
template <int N>
inline single_real<N> &single_real<N>::operator=(const char *value) {
  *this = qd_single_detail::parse_decimal<N>(value);
  if (isnan()) error("(single_real::operator=): INPUT ERROR.");
  return *this;
}
template <int N>
template <int M>
inline single_real<N> &single_real<N>::operator=(const single_real<M> &value) {
  float terms[M];
  qd_single_detail::component_copy(value, terms);
  *this = qd_single_detail::from_terms<N>(terms, M);
  return *this;
}

template <int N>
inline single_real<N> single_real<N>::add(float a, float b) {
  float error;
  const float sum = qd_single_detail::two_sum(a, b, error);
  const float terms[2] = {sum, error};
  return qd_single_detail::from_terms<N>(terms, 2);
}
template <int N>
inline single_real<N> single_real<N>::sub(float a, float b) {
  float error;
  const float difference = qd_single_detail::two_diff(a, b, error);
  const float terms[2] = {difference, error};
  return qd_single_detail::from_terms<N>(terms, 2);
}
template <int N>
inline single_real<N> single_real<N>::mul(float a, float b) {
  float error;
  const float product = qd_single_detail::two_prod(a, b, error);
  const float terms[2] = {product, error};
  return qd_single_detail::from_terms<N>(terms, 2);
}
template <int N>
inline single_real<N> single_real<N>::div(float a, float b) {
  return single_real<N>(a) / single_real<N>(b);
}
template <int N>
inline single_real<N> single_real<N>::ieee_add(const single_real &a,
                                               const single_real &b) {
  return a + b;
}
template <int N>
inline single_real<N> single_real<N>::sloppy_add(const single_real &a,
                                                 const single_real &b) {
  return a + b;
}
template <int N>
inline single_real<N> single_real<N>::bf_add(const single_real &a,
                                              const single_real &b) {
  return a + b;
}
template <int N>
inline single_real<N> single_real<N>::bf_mul(const single_real &a,
                                              const single_real &b) {
  return a * b;
}
template <int N>
inline single_real<N> single_real<N>::sloppy_div(const single_real &a,
                                                 const single_real &b) {
  return qd_single_detail::division_impl(a, b, N);
}
template <int N>
inline single_real<N> single_real<N>::accurate_div(const single_real &a,
                                                   const single_real &b) {
  return qd_single_detail::division_impl(a, b, N + 1);
}

template <int N>
inline single_real<N> sqr(const single_real<N> &a) { return a * a; }

template <int N>
inline single_real<N> single_real<N>::sqr(float value) {
  return ::sqr(single_real<N>(value));
}

template <int N>
inline single_real<N> abs(const single_real<N> &a) {
  return a.is_negative() ? -a : a;
}
template <int N>
inline single_real<N> fabs(const single_real<N> &a) { return abs(a); }

template <int N>
inline single_real<N> ldexp(const single_real<N> &a, int exponent) {
  return qd_single_detail::scale_components(a, exponent);
}

template <int N>
inline single_real<N> mul_pwr2(const single_real<N> &a, float value) {
  single_real<N> result;
  for (int i = 0; i < N; ++i) result.x[i] = a.x[i] * value;
  return result;
}

template <int N>
inline single_real<N> sqrt(const single_real<N> &a) {
  if (a.isnan()) return single_real<N>::_nan;
  if (a.isinf()) return a.is_negative() ? single_real<N>::_nan : a;
  if (a.is_zero()) return a;
  if (a.is_negative()) {
    single_real<N>::error("(single_real::sqrt): Negative argument.");
    return single_real<N>::_nan;
  }

  int exponent = 0;
  std::frexp(a.x[0], &exponent);
  const int scale = (exponent & 1) ? exponent - 1 : exponent;
  const single_real<N> scaled =
      qd_single_detail::scale_components(a, -scale);
  single_real<N> result(static_cast<float>(std::sqrt(scaled.x[0])));
  for (int i = 0; i < N + 2; ++i) {
    result += (scaled - sqr(result)) / (2.0f * result);
  }
  return qd_single_detail::scale_components(result, scale / 2);
}

template <int N>
inline single_real<N> single_real<N>::sqrt(float value) {
  return ::sqrt(single_real<N>(value));
}

template <int N>
inline single_real<N> npwr(const single_real<N> &a, int n) {
  if (n == 0) {
    if (a.is_zero()) {
      single_real<N>::error("(single_real::npwr): Invalid argument.");
      return single_real<N>::_nan;
    }
    return single_real<N>(1.0f);
  }
  single_real<N> base = a;
  single_real<N> result(1.0f);
  unsigned int exponent = static_cast<unsigned int>(n < 0 ? -n : n);
  while (exponent != 0) {
    if (exponent & 1u) result *= base;
    exponent >>= 1u;
    if (exponent != 0) base = sqr(base);
  }
  return n < 0 ? single_real<N>(1.0f) / result : result;
}

template <int N>
inline single_real<N> pow(const single_real<N> &a, int n) {
  return npwr(a, n);
}

template <int N>
inline single_real<N> nroot(const single_real<N> &a, int n) {
  if (n <= 0) {
    single_real<N>::error("(single_real::nroot): N must be positive.");
    return single_real<N>::_nan;
  }
  if (n % 2 == 0 && a.is_negative()) {
    single_real<N>::error("(single_real::nroot): Negative argument.");
    return single_real<N>::_nan;
  }
  if (n == 1) return a;
  if (n == 2) return sqrt(a);
  if (a.is_zero()) return a;
  const bool negative = a.is_negative();
  const single_real<N> magnitude = negative ? -a : a;
  single_real<N> result(qd_single_detail::from_long_double<N>(
      std::exp(std::log(qd_single_detail::to_long_double(magnitude)) / n)));
  const single_real<N> nn(static_cast<float>(n));
  for (int i = 0; i < N + 2; ++i) {
    result = ((nn - 1.0f) * result + magnitude / npwr(result, n - 1)) / nn;
  }
  return negative ? -result : result;
}

template <int N>
inline single_real<N> expm1(const single_real<N> &a) {
  if (a.isnan()) return single_real<N>::_nan;
  if (a.is_zero()) return a;
  if (abs(a) < 0.75f) {
    single_real<N> term = a;
    single_real<N> sum = a;
    single_real<N> threshold = abs(a) * single_real<N>::_eps;
    if (threshold.is_zero()) threshold = single_real<N>::_eps;
    float n = 1.0f;
    for (int i = 0; i < 1024; ++i) {
      n += 1.0f;
      term *= a;
      term /= n;
      sum += term;
      if (abs(term) <= threshold) break;
    }
    return sum;
  }
  return exp(a) - 1.0f;
}

template <int N>
inline single_real<N> exp(const single_real<N> &a) {
  if (a.isnan()) return single_real<N>::_nan;
  if (a.is_zero()) return single_real<N>(1.0f);
  if (a.isinf()) return a.is_positive() ? single_real<N>::_inf
                                        : single_real<N>(0.0f);
  if (a.x[0] <= -89.0f) return single_real<N>(0.0f);
  if (a.x[0] >= 89.0f) return single_real<N>::_inf;
  if (a.is_one()) return single_real<N>::_e;
  const long double approximate = qd_single_detail::to_long_double(a) /
      qd_single_detail::to_long_double(single_real<N>::_log2);
  const int m = static_cast<int>(std::floor(approximate + 0.5L));
  const single_real<N> reduced =
      mul_pwr2(a - single_real<N>::_log2 * static_cast<float>(m),
               0x1p-12f);
  single_real<N> s = expm1(reduced);
  for (int i = 0; i < 12; ++i) s = 2.0f * s + sqr(s);
  return ldexp(s + 1.0f, m);
}

template <int N>
inline single_real<N> log(const single_real<N> &a) {
  if (a.isnan()) return single_real<N>::_nan;
  if (a.is_zero()) return -single_real<N>::_inf;
  if (a.is_negative()) {
    single_real<N>::error("(single_real::log): Non-positive argument.");
    return single_real<N>::_nan;
  }
  if (a.isinf()) return single_real<N>::_inf;
  if (a.is_one()) return single_real<N>(0.0f);
  int exponent = 0;
  const long double mantissa = std::frexp(
      static_cast<long double>(a.x[0]), &exponent);
  single_real<N> result =
      qd_single_detail::from_long_double<N>(std::log(mantissa)) +
      single_real<N>::_log2 * static_cast<float>(exponent);
  for (int i = 0; i < N + 1; ++i) result += a * exp(-result) - 1.0f;
  return result;
}

template <int N>
inline single_real<N> log10(const single_real<N> &a) {
  return log(a) / single_real<N>::_log10;
}
template <int N>
inline single_real<N> log2(const single_real<N> &a) {
  return log(a) / single_real<N>::_log2;
}
template <int N>
inline single_real<N> exp2(const single_real<N> &a) {
  return exp(single_real<N>::_log2 * a);
}
template <int N>
inline single_real<N> pow(const single_real<N> &a,
                          const single_real<N> &b) {
  return exp(b * log(a));
}

template <int N>
inline single_real<N> log1p(const single_real<N> &a) {
  if (a.isnan()) return single_real<N>::_nan;
  if (a == -1.0) return -single_real<N>::_inf;
  if (a < -1.0) {
    single_real<N>::error("(single_real::log1p): Argument out of domain.");
    return single_real<N>::_nan;
  }
  if (a.is_zero()) return a;
  if (abs(a) < 0.75f) {
    single_real<N> term = a;
    single_real<N> sum = a;
    single_real<N> threshold = abs(a) * single_real<N>::_eps;
    if (threshold.is_zero()) threshold = single_real<N>::_eps;
    float n = 1.0f;
    for (int i = 0; i < 1024; ++i) {
      n += 1.0f;
      term *= -a;
      const single_real<N> add = term / n;
      sum += add;
      if (abs(add) <= threshold) break;
    }
    return sum;
  }
  return log(1.0f + a);
}

template <int N>
inline single_real<N> nint(const single_real<N> &a) {
  if (a.isnan() || a.isinf()) return a;
  float terms[N] = {};
  bool exact = true;
  for (int i = 0; i < N && exact; ++i) {
    const float rounded = qd_single_detail::nint_float(a.x[i]);
    terms[i] = rounded;
    exact = rounded == a.x[i];
    if (!exact && i + 1 < N &&
        std::fabs(rounded - a.x[i]) == 0.5f && a.x[i + 1] < 0.0f) {
      terms[i] = rounded - 1.0f;
    }
  }
  return qd_single_detail::from_terms<N>(terms, N);
}

template <int N>
inline single_real<N> floor(const single_real<N> &a) {
  if (a.isnan() || a.isinf()) return a;
  float terms[N] = {};
  bool exact = true;
  for (int i = 0; i < N && exact; ++i) {
    terms[i] = std::floor(a.x[i]);
    exact = terms[i] == a.x[i];
  }
  return qd_single_detail::from_terms<N>(terms, N);
}
template <int N>
inline single_real<N> ceil(const single_real<N> &a) {
  if (a.isnan() || a.isinf()) return a;
  float terms[N] = {};
  bool exact = true;
  for (int i = 0; i < N && exact; ++i) {
    terms[i] = std::ceil(a.x[i]);
    exact = terms[i] == a.x[i];
  }
  return qd_single_detail::from_terms<N>(terms, N);
}
template <int N>
inline single_real<N> aint(const single_real<N> &a) {
  return a.is_negative() ? ceil(a) : floor(a);
}
template <int N>
inline single_real<N> fmod(const single_real<N> &a,
                           const single_real<N> &b) {
  return a - b * aint(a / b);
}
template <int N>
inline single_real<N> hypot(const single_real<N> &a,
                            const single_real<N> &b) {
  if (a.isnan() || b.isnan()) return single_real<N>::_nan;
  if (a.isinf() || b.isinf()) return single_real<N>::_inf;
  single_real<N> x = abs(a);
  single_real<N> y = abs(b);
  if (x < y) std::swap(x, y);
  if (x.is_zero()) return x;
  return x * sqrt(1.0f + sqr(y / x));
}

namespace qd_single_detail {
template <int N>
inline single_real<N> pi_over_sixteen() {
  static const single_real<N> value(
      "0.196349540849362077403915211454968930262323");
  return value;
}
template <int N>
inline single_real<N> sin_pi_over_sixteen() {
  static const single_real<N> value(
      "0.195090322016128267848284868477022576602");
  return value;
}
template <int N>
inline single_real<N> cos_pi_over_sixteen() {
  static const single_real<N> value(
      "0.980785280403230449126182236134239036974");
  return value;
}
template <int N>
inline void reduce_trig_argument(const single_real<N> &a,
                                 single_real<N> &t, int &j, int &k) {
  const single_real<N> z = nint(a / single_real<N>::_2pi);
  single_real<N> remainder = a - single_real<N>::_2pi * z;
  const single_real<N> q = nint(remainder / single_real<N>::_pi2);
  j = to_int(q);
  while (j > 2) j -= 4;
  while (j < -2) j += 4;
  remainder -= single_real<N>::_pi2 * static_cast<float>(to_int(q));
  // The remaining argument is within pi/4, which is small enough for the
  // direct Taylor kernels.  Keeping the full remainder avoids a table
  // rotation whose single-step implementation would lose the sector count.
  k = 0;
  t = remainder;
}
template <int N>
inline single_real<N> sin_taylor(const single_real<N> &a) {
  if (a.is_zero()) return single_real<N>(signed_zero(a.x[0]));
  const single_real<N> threshold = abs(a) * single_real<N>::_eps;
  const single_real<N> square = -sqr(a);
  single_real<N> term = a;
  single_real<N> sum = a;
  float n = 1.0f;
  for (int i = 0; i < 128; ++i) {
    term *= square;
    term /= (n + 1.0f) * (n + 2.0f);
    sum += term;
    n += 2.0f;
    if (abs(term) <= threshold) break;
  }
  return sum;
}
template <int N>
inline single_real<N> cos_taylor(const single_real<N> &a) {
  if (a.is_zero()) return single_real<N>(1.0f);
  const single_real<N> threshold = single_real<N>::_eps;
  const single_real<N> square = -sqr(a);
  single_real<N> term(1.0f);
  single_real<N> sum(1.0f);
  float n = 0.0f;
  for (int i = 0; i < 128; ++i) {
    term *= square;
    term /= (n + 1.0f) * (n + 2.0f);
    sum += term;
    n += 2.0f;
    if (abs(term) <= threshold) break;
  }
  return sum;
}
} // namespace qd_single_detail

template <int N>
inline void sincos(const single_real<N> &a, single_real<N> &s,
                   single_real<N> &c) {
  if (a.isnan()) {
    s = c = single_real<N>::_nan;
    return;
  }
  if (a.isinf()) {
    single_real<N>::error("(single_real::sincos): Infinite argument.");
    s = c = single_real<N>::_nan;
    return;
  }
  if (a.is_zero()) {
    s = single_real<N>(qd_single_detail::signed_zero(a.x[0]));
    c = single_real<N>(1.0f);
    return;
  }
  single_real<N> t;
  int j = 0;
  int k = 0;
  qd_single_detail::reduce_trig_argument(a, t, j, k);
  if (std::abs(k) > 4) {
    single_real<N>::error("(single_real::sincos): Cannot reduce argument.");
    s = c = single_real<N>::_nan;
    return;
  }
  const single_real<N> sin_t = qd_single_detail::sin_taylor(t);
  const single_real<N> cos_t = qd_single_detail::cos_taylor(t);
  single_real<N> sin_sector = sin_t;
  single_real<N> cos_sector = cos_t;
  if (k != 0) {
    const single_real<N> u = qd_single_detail::cos_pi_over_sixteen<N>();
    const single_real<N> v = qd_single_detail::sin_pi_over_sixteen<N>();
    if (k > 0) {
      sin_sector = u * sin_t + v * cos_t;
      cos_sector = u * cos_t - v * sin_t;
    } else {
      sin_sector = u * sin_t - v * cos_t;
      cos_sector = u * cos_t + v * sin_t;
    }
  }
  if (std::abs(j) == 0) {
    s = sin_sector;
    c = cos_sector;
  } else if (j == 1) {
    s = cos_sector;
    c = -sin_sector;
  } else if (j == -1) {
    s = -cos_sector;
    c = sin_sector;
  } else {
    s = -sin_sector;
    c = -cos_sector;
  }
}
template <int N>
inline single_real<N> sin(const single_real<N> &a) {
  single_real<N> s, c;
  sincos(a, s, c);
  return s;
}
template <int N>
inline single_real<N> cos(const single_real<N> &a) {
  single_real<N> s, c;
  sincos(a, s, c);
  return c;
}
template <int N>
inline single_real<N> tan(const single_real<N> &a) {
  single_real<N> s, c;
  sincos(a, s, c);
  return s / c;
}

template <int N>
inline single_real<N> atan2(const single_real<N> &y,
                            const single_real<N> &x) {
  if (x.isnan() || y.isnan()) return single_real<N>::_nan;
  if (x.is_zero()) {
    if (y.is_zero()) {
      single_real<N>::error("(single_real::atan2): Both arguments zero.");
      return single_real<N>::_nan;
    }
    return y.is_negative() ? -single_real<N>::_pi2 : single_real<N>::_pi2;
  }
  if (y.is_zero()) {
    if (x.is_positive()) return single_real<N>(
        qd_single_detail::signed_zero(y.x[0]));
    return y.is_negative() ? -single_real<N>::_pi : single_real<N>::_pi;
  }
  if (x == y) return y.is_positive() ? single_real<N>::_pi4
                                     : -single_real<N>::_3pi4;
  if (x == -y) return y.is_positive() ? single_real<N>::_3pi4
                                      : -single_real<N>::_pi4;
  const single_real<N> radius = sqrt(sqr(x) + sqr(y));
  const single_real<N> xx = x / radius;
  const single_real<N> yy = y / radius;
  single_real<N> result(qd_single_detail::from_long_double<N>(
      std::atan2(qd_single_detail::to_long_double(y),
                 qd_single_detail::to_long_double(x))));
  single_real<N> sin_result;
  single_real<N> cos_result;
  if (std::fabs(xx.x[0]) > std::fabs(yy.x[0])) {
    for (int i = 0; i < N + 1; ++i) {
      sincos(result, sin_result, cos_result);
      result += (yy - sin_result) / cos_result;
    }
  } else {
    for (int i = 0; i < N + 1; ++i) {
      sincos(result, sin_result, cos_result);
      result -= (xx - cos_result) / sin_result;
    }
  }
  return result;
}
template <int N>
inline single_real<N> atan(const single_real<N> &a) {
  return atan2(a, single_real<N>(1.0f));
}
template <int N>
inline single_real<N> asin(const single_real<N> &a) {
  const single_real<N> magnitude = abs(a);
  if (a.isnan() || magnitude > 1.0) return single_real<N>::_nan;
  if (magnitude.is_one()) {
    return a.is_positive() ? single_real<N>::_pi2
                           : -single_real<N>::_pi2;
  }
  return atan2(a, sqrt(1.0f - sqr(a)));
}
template <int N>
inline single_real<N> acos(const single_real<N> &a) {
  const single_real<N> magnitude = abs(a);
  if (a.isnan() || magnitude > 1.0) return single_real<N>::_nan;
  if (magnitude.is_one()) {
    return a.is_positive() ? single_real<N>(0.0f) : single_real<N>::_pi;
  }
  return atan2(sqrt(1.0f - sqr(a)), a);
}

template <int N>
inline single_real<N> sinh(const single_real<N> &a) {
  if (a.isnan()) return single_real<N>::_nan;
  if (a.is_zero()) return single_real<N>(
      qd_single_detail::signed_zero(a.x[0]));
  if (a.isinf()) return a;
  if (abs(a) > 0.05) {
    const single_real<N> ea = exp(a);
    return mul_pwr2(ea - 1.0f / ea, 0.5f);
  }
  single_real<N> term = a;
  single_real<N> sum = a;
  const single_real<N> square = sqr(a);
  const single_real<N> threshold = abs(a) * single_real<N>::_eps;
  float n = 1.0f;
  for (int i = 0; i < 256; ++i) {
    n += 2.0f;
    term *= square;
    term /= (n - 1.0f) * n;
    sum += term;
    if (abs(term) <= threshold) break;
  }
  return sum;
}
template <int N>
inline single_real<N> cosh(const single_real<N> &a) {
  if (a.isnan()) return single_real<N>::_nan;
  if (a.is_zero()) return single_real<N>(1.0f);
  if (a.isinf()) return single_real<N>::_inf;
  const single_real<N> ea = exp(a);
  return mul_pwr2(ea + 1.0f / ea, 0.5f);
}
template <int N>
inline single_real<N> tanh(const single_real<N> &a) {
  if (a.isnan()) return single_real<N>::_nan;
  if (a.is_zero()) return single_real<N>(
      qd_single_detail::signed_zero(a.x[0]));
  if (a.isinf()) return a.is_positive() ? single_real<N>(1.0f)
                                        : single_real<N>(-1.0f);
  if (abs(a) > 0.05) {
    const single_real<N> ea = exp(a);
    const single_real<N> inv_ea = 1.0f / ea;
    return (ea - inv_ea) / (ea + inv_ea);
  }
  const single_real<N> s = sinh(a);
  return s / sqrt(1.0f + sqr(s));
}
template <int N>
inline void sincosh(const single_real<N> &a, single_real<N> &s,
                    single_real<N> &c) {
  s = sinh(a);
  c = cosh(a);
}
template <int N>
inline single_real<N> asinh(const single_real<N> &a) {
  return log(a + sqrt(sqr(a) + 1.0f));
}
template <int N>
inline single_real<N> acosh(const single_real<N> &a) {
  return a < 1.0 ? single_real<N>::_nan
                 : log(a + sqrt(sqr(a) - 1.0f));
}
template <int N>
inline single_real<N> atanh(const single_real<N> &a) {
  if (abs(a) > 1.0) return single_real<N>::_nan;
  if (a == 1.0) return single_real<N>::_inf;
  if (a == -1.0) return -single_real<N>::_inf;
  return 0.5f * log((1.0f + a) / (1.0f - a));
}

template <int N>
inline single_real<N> inv(const single_real<N> &a) {
  return single_real<N>(1.0f) / a;
}

template <int N>
inline single_real<N> rem(const single_real<N> &a,
                         const single_real<N> &b) {
  return a - b * aint(a / b);
}

template <int N>
inline single_real<N> drem(const single_real<N> &a,
                           const single_real<N> &b) {
  return a - b * nint(a / b);
}

template <int N>
inline single_real<N> divrem(const single_real<N> &a,
                             const single_real<N> &b,
                             single_real<N> &remainder) {
  const single_real<N> quotient = aint(a / b);
  remainder = a - b * quotient;
  return quotient;
}

template <int N>
inline single_real<N> cbrt(const single_real<N> &a) {
  if (a.isnan()) return single_real<N>::_nan;
  if (a.isinf() || a.is_zero()) return a;
  return nroot(a, 3);
}

template <int N>
inline single_real<N> trunc(const single_real<N> &a) {
  return a.is_negative() ? ceil(a) : floor(a);
}

template <int N>
inline single_real<N> round(const single_real<N> &a) {
  if (a.isnan() || a.isinf() || a.is_zero()) return a;
  return a.is_negative() ? ceil(a - 0.5f) : floor(a + 0.5f);
}

template <int N>
inline single_real<N> single_real<N>::operator^(int exponent) const {
  return npwr(*this, exponent);
}

template <int N>
inline int qd_single_compare(const single_real<N> &a,
                             const single_real<N> &b) {
  if (a.isnan() || b.isnan()) return 0;
  if (a.is_zero() && b.is_zero()) return 0;
  if (a.is_negative() != b.is_negative()) return a.is_negative() ? -1 : 1;
  for (int i = 0; i < N; ++i) {
    if (a.x[i] < b.x[i]) return -1;
    if (a.x[i] > b.x[i]) return 1;
  }
  return 0;
}

template <int N>
inline bool operator==(const single_real<N> &a,
                       const single_real<N> &b) {
  return !a.isnan() && !b.isnan() && qd_single_compare(a, b) == 0;
}
template <int N>
inline bool operator!=(const single_real<N> &a,
                       const single_real<N> &b) { return !(a == b); }
template <int N>
inline bool operator<(const single_real<N> &a,
                      const single_real<N> &b) {
  return qd_single_compare(a, b) < 0;
}
template <int N>
inline bool operator>(const single_real<N> &a,
                      const single_real<N> &b) { return b < a; }
template <int N>
inline bool operator<=(const single_real<N> &a,
                       const single_real<N> &b) { return !(a > b); }
template <int N>
inline bool operator>=(const single_real<N> &a,
                       const single_real<N> &b) { return !(a < b); }

#define QD_SINGLE_SCALAR_COMPARE(op)                                      \
  template <int N>                                                        \
  inline bool operator op(const single_real<N> &a, double b) {             \
    return a op single_real<N>(b);                                        \
  }                                                                       \
  template <int N>                                                        \
  inline bool operator op(double a, const single_real<N> &b) {             \
    return single_real<N>(a) op b;                                        \
  }

QD_SINGLE_SCALAR_COMPARE(==)
QD_SINGLE_SCALAR_COMPARE(!=)
QD_SINGLE_SCALAR_COMPARE(<)
QD_SINGLE_SCALAR_COMPARE(>)
QD_SINGLE_SCALAR_COMPARE(<=)
QD_SINGLE_SCALAR_COMPARE(>=)

#undef QD_SINGLE_SCALAR_COMPARE

template <int N>
inline bool isnan(const single_real<N> &a) { return a.isnan(); }
template <int N>
inline bool isfinite(const single_real<N> &a) { return a.isfinite(); }
template <int N>
inline bool isinf(const single_real<N> &a) { return a.isinf(); }

template <int N>
inline single_real<N> polyeval(const single_real<N> *coefficients, int degree,
                               const single_real<N> &value) {
  if (coefficients == 0 || degree < 0) return single_real<N>::_nan;
  single_real<N> result = coefficients[degree];
  for (int i = degree - 1; i >= 0; --i) {
    result *= value;
    result += coefficients[i];
  }
  return result;
}

template <int N>
inline single_real<N> polyroot(const single_real<N> *coefficients, int degree,
                               const single_real<N> &initial, int max_iter = 32,
                               double threshold = 0.0) {
  if (coefficients == 0 || degree <= 0 || max_iter <= 0) {
    single_real<N>::error("(single_real::polyroot): Invalid argument.");
    return single_real<N>::_nan;
  }
  single_real<N> derivative[64];
  if (degree > 64) {
    single_real<N>::error("(single_real::polyroot): Degree too large.");
    return single_real<N>::_nan;
  }
  long double scale = std::fabs(qd_single_detail::to_long_double(coefficients[0]));
  for (int i = 1; i <= degree; ++i) {
    derivative[i - 1] = coefficients[i] * static_cast<float>(i);
    scale = std::max(scale,
                     std::fabs(qd_single_detail::to_long_double(coefficients[i])));
  }
  const long double limit = (threshold == 0.0 ? single_real<N>::_eps
                                               : threshold) * scale;
  single_real<N> result = initial;
  for (int i = 0; i < max_iter; ++i) {
    const single_real<N> value = polyeval(coefficients, degree, result);
    if (std::fabs(qd_single_detail::to_long_double(value)) < limit) {
      return result;
    }
    const single_real<N> slope = polyeval(derivative, degree - 1, result);
    if (slope.is_zero() || slope.isnan()) break;
    result -= value / slope;
  }
  single_real<N>::error("(single_real::polyroot): Failed to converge.");
  return single_real<N>::_nan;
}

template <int N>
inline single_real<N> single_real<N>::rand() {
  single_real result(0.0f);
  float scale = 0x1p-24f;
  for (int i = 0; i < N; ++i) {
    result += static_cast<float>(std::rand()) /
              static_cast<float>(RAND_MAX) * scale;
    scale *= 0x1p-24f;
  }
  return result;
}

template <int N>
inline single_real<N> single_real<N>::debug_rand() {
  if (std::rand() % 2 == 0) return rand();
  single_real result(0.0f);
  int exponent = 0;
  for (int i = 0; i < N; ++i) {
    const float term = std::ldexp(
        static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX),
        -exponent);
    result += term;
    exponent += 24 + std::rand() % 80;
  }
  return result;
}

namespace qd_single_detail {

template <int N>
inline std::string decimal_digits(const single_real<N> &value,
                                  int precision, int &exponent) {
  const int count = std::max(1, precision);
  if (value.is_zero()) {
    exponent = 0;
    return std::string(static_cast<std::size_t>(count), '0');
  }
  single_real<N> remainder = abs(value);
  const long double estimate = std::fabs(qd_single_detail::to_long_double(remainder));
  exponent = static_cast<int>(std::floor(std::log10(estimate)));
  if (exponent > 0) {
    for (int i = 0; i < exponent; ++i) remainder /= 10.0f;
  } else if (exponent < 0) {
    for (int i = 0; i > exponent; --i) remainder *= 10.0f;
  }
  while (remainder >= 10.0f) {
    remainder /= 10.0f;
    ++exponent;
  }
  while (remainder < 1.0f) {
    remainder *= 10.0f;
    --exponent;
  }

  std::string digits;
  digits.reserve(static_cast<std::size_t>(count + 1));
  for (int i = 0; i < count + 1; ++i) {
    int digit = static_cast<int>(remainder.x[0]);
    if (digit < 0) digit = 0;
    if (digit > 9) digit = 9;
    digits.push_back(static_cast<char>('0' + digit));
    remainder = (remainder - static_cast<float>(digit)) * 10.0f;
  }
  if (digits[count] >= '5') {
    int i = count - 1;
    while (i >= 0 && digits[i] == '9') {
      digits[i] = '0';
      --i;
    }
    if (i < 0) {
      digits.insert(digits.begin(), '1');
      ++exponent;
    } else {
      ++digits[i];
    }
  }
  digits.resize(static_cast<std::size_t>(count));
  return digits;
}

inline void append_exponent(std::string &output, int exponent,
                            bool uppercase) {
  output += uppercase ? 'E' : 'e';
  output += exponent >= 0 ? '+' : '-';
  int magnitude = exponent >= 0 ? exponent : -exponent;
  std::ostringstream stream;
  stream << std::setw(2) << std::setfill('0') << magnitude;
  output += stream.str();
}

template <int N>
inline std::string format_value(const single_real<N> &value, int precision,
                                std::ios_base::fmtflags format, bool showpos,
                                bool uppercase) {
  const bool fixed = (format & std::ios_base::fixed) != 0;
  std::string output;
  if (value.isnan()) return uppercase ? "NAN" : "nan";
  const bool negative = value.is_negative() ||
      (value.is_zero() && std::signbit(value.x[0]));
  if (negative) output += '-';
  else if (showpos) output += '+';
  if (value.isinf()) return output + (uppercase ? "INF" : "inf");
  if (value.is_zero()) {
    output += '0';
    if (precision > 0) {
      output += '.';
      output.append(static_cast<std::size_t>(precision), '0');
    }
    return output;
  }

  int exponent = 0;
  const int requested = fixed
      ? std::max(1, precision + (std::max(0, exponent) + 1))
      : std::max(1, precision + 1);
  // Estimate the exponent before choosing the fixed-point digit count.
  const long double estimate = std::fabs(qd_single_detail::to_long_double(value));
  const int estimated_exponent =
      static_cast<int>(std::floor(std::log10(estimate)));
  const int count = fixed
      ? std::max(1, precision + (estimated_exponent >= 0
                                     ? estimated_exponent + 1
                                     : -estimated_exponent))
      : requested;
  const std::string digits = decimal_digits(value, count, exponent);
  if (fixed) {
    const int point = exponent + 1;
    if (point <= 0) {
      output += "0.";
      output.append(static_cast<std::size_t>(-point), '0');
      output += digits;
    } else if (point >= static_cast<int>(digits.size())) {
      output += digits;
      output.append(static_cast<std::size_t>(point - digits.size()), '0');
      if (precision > 0) {
        output += '.';
        output.append(static_cast<std::size_t>(precision), '0');
      }
    } else {
      output.append(digits, 0, static_cast<std::size_t>(point));
      if (precision > 0) {
        output += '.';
        output.append(digits, static_cast<std::size_t>(point),
                      std::string::npos);
      }
    }
    return output;
  }

  output += digits[0];
  if (precision > 0) {
    output += '.';
    output.append(digits, 1, std::string::npos);
  }
  append_exponent(output, exponent, uppercase);
  return output;
}

} // namespace qd_single_detail

template <int N>
inline void single_real<N>::to_digits(char *output, int &exponent,
                                      int precision) const {
  if (output == 0 || precision <= 0) {
    exponent = 0;
    return;
  }
  const std::string digits = qd_single_detail::decimal_digits(
      *this, precision, exponent);
  for (int i = 0; i < precision; ++i) output[i] = digits[i];
  output[precision] = 0;
}

template <int N>
inline void single_real<N>::write(char *output, int length, int precision,
                                  bool showpos, bool uppercase) const {
  if (output == 0 || length <= 0) return;
  const std::string text = qd_single_detail::format_value(
      *this, precision, std::ios_base::scientific, showpos, uppercase);
  const std::size_t count = std::min<std::size_t>(
      text.size(), static_cast<std::size_t>(length - 1));
  for (std::size_t i = 0; i < count; ++i) output[i] = text[i];
  output[count] = 0;
}

template <int N>
inline std::string single_real<N>::to_string(
    int precision, int width, std::ios_base::fmtflags format, bool showpos,
    bool uppercase, char fill) const {
  std::string text = qd_single_detail::format_value(
      *this, precision, format, showpos, uppercase);
  if (width > static_cast<int>(text.size())) {
    const std::size_t count = static_cast<std::size_t>(width - text.size());
    if (format & std::ios_base::left) text.append(count, fill);
    else if (format & std::ios_base::internal &&
             (text[0] == '+' || text[0] == '-'))
      text.insert(1, count, fill);
    else text.insert(0, count, fill);
  }
  return text;
}

template <int N>
inline int single_real<N>::read(const char *input, single_real<N> &value) {
  if (input == 0) return -1;
  const single_real<N> parsed = qd_single_detail::parse_decimal<N>(input);
  std::string lower(input);
  std::transform(lower.begin(), lower.end(), lower.begin(), [](char c) {
    return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  const bool special = lower == "nan" || lower == "+nan" || lower == "-nan" ||
      lower == "inf" || lower == "+inf" || lower == "-inf" ||
      lower == "infinity" || lower == "+infinity" || lower == "-infinity";
  if (parsed.isnan() && !special) return -1;
  value = parsed;
  return 0;
}

template <int N>
inline void single_real<N>::dump(const std::string &name,
                                 std::ostream &stream) const {
  if (!name.empty()) stream << name << " = ";
  stream << "[ ";
  for (int i = 0; i < N; ++i) {
    if (i != 0) stream << ", ";
    stream << std::scientific << std::setprecision(9) << std::setw(16)
           << x[i];
  }
  stream << " ]" << std::endl;
}

template <int N>
inline void single_real<N>::dump_bits(const std::string &name,
                                      std::ostream &stream) const {
  if (!name.empty()) stream << name << " = ";
  stream << "[ ";
  for (int i = 0; i < N; ++i) {
    if (i != 0) stream << ", ";
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(x[i]), "float must be binary32");
    std::memcpy(&bits, &x[i], sizeof(bits));
    stream << "0x" << std::hex << bits << std::dec;
  }
  stream << " ]" << std::endl;
}

template <int N>
inline std::ostream &operator<<(std::ostream &stream,
                                const single_real<N> &value) {
  const bool showpos = (stream.flags() & std::ios_base::showpos) != 0;
  const bool uppercase = (stream.flags() & std::ios_base::uppercase) != 0;
  return stream << value.to_string(static_cast<int>(stream.precision()),
                                   static_cast<int>(stream.width()),
                                   stream.flags(), showpos, uppercase,
                                   stream.fill());
}

template <int N>
inline std::istream &operator>>(std::istream &stream, single_real<N> &value) {
  std::string text;
  stream >> text;
  if (!stream.fail() && value.read(text.c_str(), value) != 0) {
    stream.setstate(std::ios_base::failbit);
  }
  return stream;
}

template <int N>
inline double to_double(const single_real<N> &value) {
  return static_cast<double>(qd_single_detail::to_long_double(value));
}
template <int N>
inline long double to_long_double(const single_real<N> &value) {
  return qd_single_detail::to_long_double(value);
}
template <int N>
inline int to_int(const single_real<N> &value) {
  return static_cast<int>(qd_single_detail::to_long_double(value));
}
template <int N>
inline long to_long(const single_real<N> &value) {
  return static_cast<long>(qd_single_detail::to_long_double(value));
}
template <int N>
inline unsigned long to_unsigned_long(const single_real<N> &value) {
  return static_cast<unsigned long>(qd_single_detail::to_long_double(value));
}
template <int N>
inline long long to_long_long(const single_real<N> &value) {
  return static_cast<long long>(qd_single_detail::to_long_double(value));
}
template <int N>
inline unsigned long long to_unsigned_long_long(const single_real<N> &value) {
  return static_cast<unsigned long long>(
      qd_single_detail::to_long_double(value));
}
template <int N>
inline std::int64_t to_int64_t(const single_real<N> &value) {
  return static_cast<std::int64_t>(qd_single_detail::to_long_double(value));
}
template <int N>
inline std::uint64_t to_uint64_t(const single_real<N> &value) {
  return static_cast<std::uint64_t>(qd_single_detail::to_long_double(value));
}

namespace std {
template <int N>
class numeric_limits<single_real<N> > : public numeric_limits<float> {
 public:
  static inline float epsilon() { return single_real<N>::_eps; }
  static inline single_real<N> max() { return single_real<N>::_max; }
  static inline single_real<N> safe_max() { return single_real<N>::_safe_max; }
  static inline float min() { return single_real<N>::_min_normalized; }
  static const int digits = qd_single_detail::traits<N>::digits;
  static const int digits10 = qd_single_detail::traits<N>::digits10;
};
} // namespace std

template <int N>
inline single_real<N> mul_pwr2(const single_real<N> &value, double scale) {
  return ldexp(value, static_cast<int>(std::log2(scale)));
}

#endif /* _QD_SINGLE_REAL_DIRECT_H */
