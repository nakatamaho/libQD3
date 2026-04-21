/*
 * include/edd_real.h
 *
 * Extended-double-double precision arithmetic package based on a native
 * two-limb _Float64x expansion for GNU C++ binary80 targets.
 */
#ifndef _QD_EDD_REAL_H
#define _QD_EDD_REAL_H

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

#include <qd/fpu.h>
#include <qd/qd_config.h>
#include <qd/dd_real.h>
#include <qd/qd_real.h>

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

#ifndef QD_HAVE_EDD_REAL
#error "edd_real requires a configure result with GNU C++ _Float64x binary80 support enabled."
#endif

#if QD_EDD_FLT64X_MANT_DIG != 64
#error "edd_real requires _Float64x with binary80 mantissa precision (64 bits)."
#endif

#if QD_EDD_FLT64X_MAX_EXP != 16384
#error "edd_real requires _Float64x with binary80 exponent range (max exp 16384)."
#endif

using edd_word = _Float64x;

struct QD_API edd_real {
  edd_word x[2];

  edd_real() {
    x[0] = (edd_word) 0.0;
    x[1] = (edd_word) 0.0;
  }

  edd_real(edd_word hi, edd_word lo) {
    x[0] = hi;
    x[1] = lo;
  }

  explicit edd_real(const edd_word *xx) {
    x[0] = xx[0];
    x[1] = xx[1];
  }

  edd_real(edd_word d) {
    x[0] = d;
    x[1] = (edd_word) 0.0;
  }

  edd_real(double d) {
    x[0] = (_Float64x) d;
    x[1] = (_Float64x) 0.0;
  }

  edd_real(int i) {
    x[0] = (edd_word) i;
    x[1] = (edd_word) 0.0;
  }

  edd_real(const char *s);
  edd_real(const dd_real &dd);
  explicit edd_real(const qd_real &qd);

  edd_word operator[](int i) const { return x[i]; }
  edd_word &operator[](int i) { return x[i]; }

  edd_word _hi() const { return x[0]; }
  edd_word _lo() const { return x[1]; }

  static void error(const char *msg);

  static const edd_real _2pi;
  static const edd_real _pi;
  static const edd_real _3pi4;
  static const edd_real _pi2;
  static const edd_real _pi4;
  static const edd_real _e;
  static const edd_real _log2;
  static const edd_real _log10;
  static const edd_real _nan;
  static const edd_real _inf;
  static const edd_real _max;
  static const edd_real _safe_max;

  static const edd_word _eps;
  static const edd_word _min_normalized;
  static const int _ndigits;

  bool isnan() const;
  bool isfinite() const;
  bool isinf() const;

  static edd_real add(edd_word a, edd_word b);
  static edd_real sub(edd_word a, edd_word b);
  static edd_real mul(edd_word a, edd_word b);
  static edd_real div(edd_word a, edd_word b);
  static edd_real sqr(edd_word a);
  static edd_real sqrt(edd_word a);

  edd_real &operator+=(edd_word a);
  edd_real &operator+=(const edd_real &a);
  edd_real &operator-=(edd_word a);
  edd_real &operator-=(const edd_real &a);
  edd_real &operator*=(edd_word a);
  edd_real &operator*=(const edd_real &a);
  edd_real &operator/=(edd_word a);
  edd_real &operator/=(const edd_real &a);

  edd_real operator-() const;

  edd_real &operator=(edd_word a);
  edd_real &operator=(double a);
  edd_real &operator=(int a);
  edd_real &operator=(const char *s);
  edd_real &operator=(const dd_real &a);
  edd_real &operator=(const qd_real &a);

  edd_real operator^(int n) const;

  bool is_zero() const;
  bool is_one() const;
  bool is_positive() const;
  bool is_negative() const;

  void to_digits(char *s, int &expn, int precision = _ndigits) const;
  void write(char *s, int len, int precision = _ndigits,
      bool showpos = false, bool uppercase = false) const;
  std::string to_string(int precision = _ndigits, int width = 0,
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0),
      bool showpos = false, bool uppercase = false, char fill = ' ') const;
  static int read(const char *s, edd_real &a);
};

namespace std {
  template <>
  class numeric_limits<edd_real> : public numeric_limits<long double> {
  public:
    inline static _Float64x epsilon() { return edd_real::_eps; }
    inline static _Float64x min() { return edd_real::_min_normalized; }
    inline static edd_real max() { return edd_real::_max; }
    inline static edd_real safe_max() { return edd_real::_safe_max; }
    static const int digits = 126;
    static const int digits10 = 37;
  };
}

QD_API inline bool isnan(const edd_real &a) { return a.isnan(); }
QD_API inline bool isfinite(const edd_real &a) { return a.isfinite(); }
QD_API inline bool isinf(const edd_real &a) { return a.isinf(); }

QD_API edd_real operator+(const edd_real &a, edd_word b);
QD_API edd_real operator+(edd_word a, const edd_real &b);
QD_API edd_real operator+(const edd_real &a, const edd_real &b);

QD_API edd_real operator-(const edd_real &a, edd_word b);
QD_API edd_real operator-(edd_word a, const edd_real &b);
QD_API edd_real operator-(const edd_real &a, const edd_real &b);

QD_API edd_real operator*(const edd_real &a, edd_word b);
QD_API edd_real operator*(edd_word a, const edd_real &b);
QD_API edd_real operator*(const edd_real &a, const edd_real &b);

QD_API edd_real operator/(const edd_real &a, edd_word b);
QD_API edd_real operator/(edd_word a, const edd_real &b);
QD_API edd_real operator/(const edd_real &a, const edd_real &b);

QD_API edd_real sqr(const edd_real &a);
QD_API edd_real sqrt(const edd_real &a);
QD_API edd_real npwr(const edd_real &a, int n);
QD_API edd_real pow(const edd_real &a, int n);

QD_API edd_real abs(const edd_real &a);
QD_API edd_real fabs(const edd_real &a);

QD_API edd_real ldexp(const edd_real &a, int n);
QD_API edd_real mul_pwr2(const edd_real &a, edd_word b);

QD_API bool operator==(const edd_real &a, edd_word b);
QD_API bool operator==(edd_word a, const edd_real &b);
QD_API bool operator==(const edd_real &a, const edd_real &b);

QD_API bool operator!=(const edd_real &a, edd_word b);
QD_API bool operator!=(edd_word a, const edd_real &b);
QD_API bool operator!=(const edd_real &a, const edd_real &b);

QD_API bool operator<(const edd_real &a, edd_word b);
QD_API bool operator<(edd_word a, const edd_real &b);
QD_API bool operator<(const edd_real &a, const edd_real &b);

QD_API bool operator>(const edd_real &a, edd_word b);
QD_API bool operator>(edd_word a, const edd_real &b);
QD_API bool operator>(const edd_real &a, const edd_real &b);

QD_API bool operator<=(const edd_real &a, edd_word b);
QD_API bool operator<=(edd_word a, const edd_real &b);
QD_API bool operator<=(const edd_real &a, const edd_real &b);

QD_API bool operator>=(const edd_real &a, edd_word b);
QD_API bool operator>=(edd_word a, const edd_real &b);
QD_API bool operator>=(const edd_real &a, const edd_real &b);

QD_API edd_word to_float64x(const edd_real &a);
QD_API dd_real to_dd_real(const edd_real &a);
QD_API edd_real to_edd_real(const dd_real &a);
QD_API edd_real to_edd_real(const qd_real &a);
QD_API qd_real to_qd_real(const edd_real &a);
QD_API double to_double(const edd_real &a);
QD_API int to_int(const edd_real &a);

QD_API std::ostream &operator<<(std::ostream &s, const edd_real &a);
QD_API std::istream &operator>>(std::istream &s, edd_real &a);

#ifdef QD_INLINE
#include <qd/edd_inline.h>
#endif

QD_API extern bool edd_suppress_error_messages;

#endif /* _QD_EDD_REAL_H */
