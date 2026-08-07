/*
 * include/dd_real.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Double-double precision (>= 106-bit significand) floating point
 * arithmetic package based on David Bailey's Fortran-90 double-double
 * package, with some changes. See  
 *
 *   http://www.nersc.gov/~dhbailey/mpdist/mpdist.html
 *   
 * for the original Fortran-90 version.
 *
 * Overall structure is similar to that of Keith Brigg's C++ double-double
 * package.  See  
 *
 *   http://www-epidem.plansci.cam.ac.uk/~kbriggs/doubledouble.html
 *
 * for more details.  In particular, the fix for x86 computers is borrowed
 * from his code.
 *
 * Yozo Hida
 */

#ifndef _QD_DD_REAL_H
#define _QD_DD_REAL_H

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <limits>
#include <type_traits>
#include <qd/qd_config.h>
#include <qd/fpu.h>
#include <qd/inline.h>

#if defined(QD_BF)
#  if !defined(QD_BF_ADD)
#    define QD_BF_ADD
#  endif
#  if !defined(QD_BF_MUL)
#    define QD_BF_MUL
#  endif
#endif

// Some compilers define isnan, isfinite, and isinf as macros, even for
// C++ codes, which cause havoc when overloading these functions.  We undef
// them here.
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

struct td_real;

struct QD_API dd_real {
  double x[2];

  dd_real(double hi, double lo) { x[0] = hi; x[1] = lo; }
  dd_real() {x[0] = 0.0; x[1] = 0.0; }
  dd_real(double h) { x[0] = h; x[1] = 0.0; }
  dd_real(int h) {
    qd::int64_to_double_expansion(static_cast<std::int64_t>(h), x, 2);
  }
  dd_real(std::uint64_t h);
  dd_real(std::int64_t h);
  template <class I,
            typename std::enable_if<qd::is_supported_integral<I>::value,
                                    int>::type = 0>
  dd_real(I h) {
    qd::integer_to_double_expansion(h, x, 2);
  }

  dd_real (const char *s);
  explicit dd_real (const double *d) {
    x[0] = d[0]; x[1] = d[1];
  }

  static void error(const char *msg);

  double _hi() const { return x[0]; }
  double _lo() const { return x[1]; }

  static const dd_real _2pi;
  static const dd_real _pi;
  static const dd_real _3pi4;
  static const dd_real _pi2;
  static const dd_real _pi4;
  static const dd_real _e;
  static const dd_real _log2;
  static const dd_real _log10;
  static const dd_real _nan;
  static const dd_real _inf;

  static const double _eps;
  static const double _min_normalized;
  static const dd_real _max;
  static const dd_real _safe_max;
  static const int _ndigits;

  bool isnan() const { return QD_ISNAN(x[0]) || QD_ISNAN(x[1]); }
  bool isfinite() const { return QD_ISFINITE(x[0]); }
  bool isinf() const { return QD_ISINF(x[0]); }

  static dd_real add(double a, double b);
  static dd_real ieee_add(const dd_real &a, const dd_real &b);
  static dd_real sloppy_add(const dd_real &a, const dd_real &b);
  static dd_real bf_add(const dd_real &a, const dd_real &b);

  dd_real &operator+=(double a);
  dd_real &operator+=(const dd_real &a);
  dd_real &operator+=(const td_real &a);

  static dd_real sub(double a, double b);

  dd_real &operator-=(double a);
  dd_real &operator-=(const dd_real &a);
  dd_real &operator-=(const td_real &a);

  dd_real operator-() const;

  static dd_real mul(double a, double b);
  static dd_real bf_mul(const dd_real &a, const dd_real &b);

  dd_real &operator*=(double a);
  dd_real &operator*=(const dd_real &a);
  dd_real &operator*=(const td_real &a);

  static dd_real div(double a, double b);
  static dd_real sloppy_div(const dd_real &a, const dd_real &b);
  static dd_real accurate_div(const dd_real &a, const dd_real &b);
  
  dd_real &operator/=(double a);
  dd_real &operator/=(const dd_real &a);
  dd_real &operator/=(const td_real &a);

  dd_real &operator=(double a);
  dd_real &operator=(std::uint64_t a);
  dd_real &operator=(std::int64_t a);
  template <class I>
  typename std::enable_if<qd::is_supported_integral<I>::value,
                          dd_real &>::type
  operator=(I a) {
    qd::integer_to_double_expansion(a, x, 2);
    return *this;
  }
  dd_real &operator=(const char *s);

  template <typename Integer,
            typename std::enable_if<
                std::is_same<Integer, int>::value ||
                std::is_same<Integer, long>::value ||
                std::is_same<Integer, unsigned long>::value ||
                std::is_same<Integer, long long>::value ||
                std::is_same<Integer, unsigned long long>::value ||
                std::is_same<Integer, std::int64_t>::value ||
                std::is_same<Integer, std::uint64_t>::value,
                int>::type = 0>
  operator Integer() const {
    return static_cast<Integer>(x[0]);
  }

  dd_real operator^(int n);
  static dd_real sqr(double d);

  static dd_real sqrt(double a);
  
  bool is_zero() const;
  bool is_one() const;
  bool is_positive() const;
  bool is_negative() const;

  static dd_real rand(void);

  void to_digits(char *s, int &expn, int precision = _ndigits) const;
  void write(char *s, int len, int precision = _ndigits, 
      bool showpos = false, bool uppercase = false) const;
  std::string to_string(int precision = _ndigits, int width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;
  int read(const char *s, dd_real &a);

  /* Debugging Methods */
  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
  void dump_bits(const std::string &name = "", 
                 std::ostream &os = std::cerr) const;

  static dd_real debug_rand();
};


namespace std {
  template <>
  class numeric_limits<dd_real> : public numeric_limits<double> {
  public:
    inline static double epsilon() { return dd_real::_eps; }
    inline static dd_real max() { return dd_real::_max; }
    inline static dd_real safe_max() { return dd_real::_safe_max; }
    inline static double min() { return dd_real::_min_normalized; }
    static const int digits = 104;
    static const int digits10 = 31;
  };
}

QD_API dd_real ddrand(void);
QD_API dd_real sqrt(const dd_real &a);

QD_API dd_real polyeval(const dd_real *c, int n, const dd_real &x);
QD_API dd_real polyroot(const dd_real *c, int n, 
    const dd_real &x0, int max_iter = 32, double thresh = 0.0);

QD_API inline bool isnan(const dd_real &a) { return a.isnan(); }
QD_API inline bool isfinite(const dd_real &a) { return a.isfinite(); }
QD_API inline bool isinf(const dd_real &a) { return a.isinf(); }

/* Computes  dd * d  where d is known to be a power of 2. */
QD_API dd_real mul_pwr2(const dd_real &dd, double d);

QD_API dd_real operator+(const dd_real &a, double b);
QD_API dd_real operator+(double a, const dd_real &b);
QD_API dd_real operator+(const dd_real &a, const dd_real &b);
QD_API td_real operator+(const dd_real &a, const td_real &b);
QD_API td_real operator+(const td_real &a, const dd_real &b);

QD_API dd_real operator-(const dd_real &a, double b);
QD_API dd_real operator-(double a, const dd_real &b);
QD_API dd_real operator-(const dd_real &a, const dd_real &b);
QD_API td_real operator-(const dd_real &a, const td_real &b);
QD_API td_real operator-(const td_real &a, const dd_real &b);

QD_API dd_real operator*(const dd_real &a, double b);
QD_API dd_real operator*(double a, const dd_real &b);
QD_API dd_real operator*(const dd_real &a, const dd_real &b);
QD_API td_real operator*(const dd_real &a, const td_real &b);
QD_API td_real operator*(const td_real &a, const dd_real &b);

QD_API dd_real operator/(const dd_real &a, double b);
QD_API dd_real operator/(double a, const dd_real &b);
QD_API dd_real operator/(const dd_real &a, const dd_real &b);
QD_API td_real operator/(const dd_real &a, const td_real &b);
QD_API td_real operator/(const td_real &a, const dd_real &b);

QD_API dd_real inv(const dd_real &a);

QD_API dd_real rem(const dd_real &a, const dd_real &b);
QD_API dd_real drem(const dd_real &a, const dd_real &b);
QD_API dd_real divrem(const dd_real &a, const dd_real &b, dd_real &r);

QD_API dd_real pow(const dd_real &a, int n);
QD_API dd_real pow(const dd_real &a, const dd_real &b);
QD_API dd_real npwr(const dd_real &a, int n);
QD_API dd_real sqr(const dd_real &a);

QD_API dd_real sqrt(const dd_real &a);
QD_API dd_real nroot(const dd_real &a, int n);

QD_API bool operator==(const dd_real &a, double b);
QD_API bool operator==(double a, const dd_real &b);
QD_API bool operator==(const dd_real &a, const dd_real &b);
QD_API bool operator==(const dd_real &a, const td_real &b);
QD_API bool operator==(const td_real &a, const dd_real &b);

QD_API bool operator<=(const dd_real &a, double b);
QD_API bool operator<=(double a, const dd_real &b);
QD_API bool operator<=(const dd_real &a, const dd_real &b);
QD_API bool operator<=(const dd_real &a, const td_real &b);
QD_API bool operator<=(const td_real &a, const dd_real &b);

QD_API bool operator>=(const dd_real &a, double b);
QD_API bool operator>=(double a, const dd_real &b);
QD_API bool operator>=(const dd_real &a, const dd_real &b);
QD_API bool operator>=(const dd_real &a, const td_real &b);
QD_API bool operator>=(const td_real &a, const dd_real &b);

QD_API bool operator<(const dd_real &a, double b);
QD_API bool operator<(double a, const dd_real &b);
QD_API bool operator<(const dd_real &a, const dd_real &b);
QD_API bool operator<(const dd_real &a, const td_real &b);
QD_API bool operator<(const td_real &a, const dd_real &b);

QD_API bool operator>(const dd_real &a, double b);
QD_API bool operator>(double a, const dd_real &b);
QD_API bool operator>(const dd_real &a, const dd_real &b);
QD_API bool operator>(const dd_real &a, const td_real &b);
QD_API bool operator>(const td_real &a, const dd_real &b);

QD_API bool operator!=(const dd_real &a, double b);
QD_API bool operator!=(double a, const dd_real &b);
QD_API bool operator!=(const dd_real &a, const dd_real &b);
QD_API bool operator!=(const dd_real &a, const td_real &b);
QD_API bool operator!=(const td_real &a, const dd_real &b);

template <typename Integer>
struct dd_real_integer_operand
    : std::integral_constant<bool,
          std::is_same<Integer, int>::value ||
          std::is_same<Integer, unsigned int>::value ||
          std::is_same<Integer, long>::value ||
          std::is_same<Integer, unsigned long>::value ||
          std::is_same<Integer, long long>::value ||
          std::is_same<Integer, unsigned long long>::value ||
          std::is_same<Integer, std::int64_t>::value ||
          std::is_same<Integer, std::uint64_t>::value> {};

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator==(const dd_real &a, Integer b) {
  return a == dd_real(b);
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator==(Integer a, const dd_real &b) {
  return dd_real(a) == b;
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator!=(const dd_real &a, Integer b) {
  return a != dd_real(b);
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator!=(Integer a, const dd_real &b) {
  return dd_real(a) != b;
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator<(const dd_real &a, Integer b) {
  return a < dd_real(b);
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator<(Integer a, const dd_real &b) {
  return dd_real(a) < b;
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator>(const dd_real &a, Integer b) {
  return a > dd_real(b);
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator>(Integer a, const dd_real &b) {
  return dd_real(a) > b;
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator<=(const dd_real &a, Integer b) {
  return a <= dd_real(b);
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator<=(Integer a, const dd_real &b) {
  return dd_real(a) <= b;
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator>=(const dd_real &a, Integer b) {
  return a >= dd_real(b);
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline bool operator>=(Integer a, const dd_real &b) {
  return dd_real(a) >= b;
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline dd_real operator+(const dd_real &a, Integer b) {
  return a + dd_real(b);
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline dd_real operator+(Integer a, const dd_real &b) {
  return dd_real(a) + b;
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline dd_real operator-(const dd_real &a, Integer b) {
  return a - dd_real(b);
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline dd_real operator-(Integer a, const dd_real &b) {
  return dd_real(a) - b;
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline dd_real operator*(const dd_real &a, Integer b) {
  return a * dd_real(b);
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline dd_real operator*(Integer a, const dd_real &b) {
  return dd_real(a) * b;
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline dd_real operator/(const dd_real &a, Integer b) {
  return a / dd_real(b);
}

template <typename Integer,
          typename std::enable_if<dd_real_integer_operand<Integer>::value, int>::type = 0>
inline dd_real operator/(Integer a, const dd_real &b) {
  return dd_real(a) / b;
}

QD_API dd_real nint(const dd_real &a);
QD_API dd_real floor(const dd_real &a);
QD_API dd_real ceil(const dd_real &a);
QD_API dd_real aint(const dd_real &a);

QD_API dd_real ddrand(void);

double to_double(const dd_real &a);
int    to_int(const dd_real &a);
long   to_long(const dd_real &a);
unsigned long to_unsigned_long(const dd_real &a);
long long to_long_long(const dd_real &a);
unsigned long long to_unsigned_long_long(const dd_real &a);
std::int64_t to_int64_t(const dd_real &a);
std::uint64_t to_uint64_t(const dd_real &a);

QD_API dd_real exp(const dd_real &a);
QD_API dd_real ldexp(const dd_real &a, int exp);
QD_API dd_real log(const dd_real &a);
QD_API dd_real log10(const dd_real &a);
QD_API dd_real log2(const dd_real &a);
QD_API dd_real exp2(const dd_real &a);
QD_API dd_real expm1(const dd_real &a);
QD_API dd_real log1p(const dd_real &a);

QD_API dd_real sin(const dd_real &a);
QD_API dd_real cos(const dd_real &a);
QD_API dd_real tan(const dd_real &a);
QD_API void sincos(const dd_real &a, dd_real &sin_a, dd_real &cos_a);

QD_API dd_real asin(const dd_real &a);
QD_API dd_real acos(const dd_real &a);
QD_API dd_real atan(const dd_real &a);
QD_API dd_real atan2(const dd_real &y, const dd_real &x);

QD_API dd_real sinh(const dd_real &a);
QD_API dd_real cosh(const dd_real &a);
QD_API dd_real tanh(const dd_real &a);
QD_API void sincosh(const dd_real &a, 
                      dd_real &sinh_a, dd_real &cosh_a);

QD_API dd_real asinh(const dd_real &a);
QD_API dd_real acosh(const dd_real &a);
QD_API dd_real atanh(const dd_real &a);

QD_API dd_real fabs(const dd_real &a);
QD_API dd_real abs(const dd_real &a);   /* same as fabs */

QD_API dd_real fmod(const dd_real &a, const dd_real &b);
QD_API dd_real hypot(const dd_real &a, const dd_real &b);
QD_API dd_real cbrt(const dd_real &a);
QD_API dd_real trunc(const dd_real &a);
QD_API dd_real round(const dd_real &a);

QD_API std::ostream& operator<<(std::ostream &s, const dd_real &a);
QD_API std::istream& operator>>(std::istream &s, dd_real &a);
#ifdef QD_INLINE
#include <qd/dd_inline.h>
#endif



QD_API extern bool dd_suppress_error_messages; /** Set to true to suppress error messages. Initialized in dd_real.cpp */

#endif /* _QD_DD_REAL_H */
