/*
 * include/inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This file contains the basic functions used both by double-double
 * and quad-double package.  These are declared as inline functions as
 * they are the smallest building blocks of the double-double and 
 * quad-double arithmetic.
 */
#ifndef _QD_INLINE_H
#define _QD_INLINE_H

#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

#ifdef QD_VACPP_BUILTINS_H
/* For VisualAge C++ __fmadd */
#include <builtins.h>
#endif

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace qd {

static const double _d_nan = std::numeric_limits<double>::quiet_NaN();
static const double _d_inf = std::numeric_limits<double>::infinity();

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
inline double quick_two_sum(double a, double b, double &err) {
  double s = a + b;
  err = b - (s - a);
  return s;
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
inline double quick_two_diff(double a, double b, double &err) {
  double s = a - b;
  err = (a - s) - b;
  return s;
}

/* Computes fl(a+b) and err(a+b).  */
inline double two_sum(double a, double b, double &err) {
  double s = a + b;
  double bb = s - a;
  err = (a - (s - bb)) + (b - bb);
  return s;
}

/* Computes fl(a-b) and err(a-b).  */
inline double two_diff(double a, double b, double &err) {
  double s = a - b;
  double bb = s - a;
  err = (a - (s - bb)) - (b + bb);
  return s;
}

inline void uint64_to_double_expansion(std::uint64_t value, double *out,
                                       int limbs) {
  const double hi =
      std::ldexp(static_cast<double>(value >> 32), 32);
  const double lo = static_cast<double>(
      value & static_cast<std::uint64_t>(0xffffffffu));
  double err = 0.0;
  const double sum = qd::two_sum(hi, lo, err);

  if (limbs < 2 && err != 0.0) {
    throw std::range_error("uint64_t constructor would round");
  }

  if (limbs > 0) {
    out[0] = sum;
  }
  if (limbs > 1) {
    out[1] = err;
  }
  for (int i = 2; i < limbs; i++) {
    out[i] = 0.0;
  }
}

inline void int64_to_double_expansion(std::int64_t value, double *out,
                                      int limbs) {
  const bool negative = value < 0;
  const std::uint64_t magnitude = negative
      ? static_cast<std::uint64_t>(-(value + 1)) + 1u
      : static_cast<std::uint64_t>(value);

  qd::uint64_to_double_expansion(magnitude, out, limbs);
  if (negative) {
    for (int i = 0; i < limbs; i++) {
      out[i] = -out[i];
    }
  }
}

template <class I>
struct integral_conversion_type {
  typedef typename std::remove_cv<
      typename std::remove_reference<I>::type>::type type;
};

template <class I,
          bool IsIntegral =
              std::is_integral<typename integral_conversion_type<I>::type>::value>
struct is_supported_integral : std::false_type {};

template <class I>
struct is_supported_integral<I, true>
    : std::integral_constant<
          bool,
          !std::is_same<typename integral_conversion_type<I>::type, bool>::value &&
              (sizeof(typename integral_conversion_type<I>::type) <=
               sizeof(std::uint64_t))> {};

template <class I, bool IsSupported = is_supported_integral<I>::value>
struct is_supported_signed_integral : std::false_type {};

template <class I>
struct is_supported_signed_integral<I, true>
    : std::integral_constant<
          bool,
          std::is_signed<typename integral_conversion_type<I>::type>::value> {};

template <class I, bool IsSupported = is_supported_integral<I>::value>
struct is_supported_unsigned_integral : std::false_type {};

template <class I>
struct is_supported_unsigned_integral<I, true>
    : std::integral_constant<
          bool,
          !std::is_signed<typename integral_conversion_type<I>::type>::value> {};

template <class I>
inline typename std::enable_if<is_supported_signed_integral<I>::value,
                               void>::type
integer_to_double_expansion(I value, double *out, int limbs) {
  qd::int64_to_double_expansion(static_cast<std::int64_t>(value), out, limbs);
}

template <class I>
inline typename std::enable_if<is_supported_unsigned_integral<I>::value,
                               void>::type
integer_to_double_expansion(I value, double *out, int limbs) {
  qd::uint64_to_double_expansion(static_cast<std::uint64_t>(value), out, limbs);
}

#ifndef QD_FMS
/* Computes high word and lo word of a */
inline void split(double a, double &hi, double &lo) {
  double temp;
  if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
    a *= 3.7252902984619140625e-09;  // 2^-28
    temp = _QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
    hi *= 268435456.0;          // 2^28
    lo *= 268435456.0;          // 2^28
  } else {
    temp = _QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
  }
}
#endif

/* Computes fl(a*b) and err(a*b). */
inline double two_prod(double a, double b, double &err) {
#ifdef QD_FMS
  double p = a * b;
  err = QD_FMS(a, b, p);
  return p;
#else
  double a_hi, a_lo, b_hi, b_lo;
  double p = a * b;
  split(a, a_hi, a_lo);
  split(b, b_hi, b_lo);
  err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
  return p;
#endif
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
inline double two_sqr(double a, double &err) {
#ifdef QD_FMS
  double p = a * a;
  err = QD_FMS(a, a, p);
  return p;
#else
  double hi, lo;
  double q = a * a;
  split(a, hi, lo);
  err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;
  return q;
#endif
}

/* Computes the nearest integer to d. */
inline double nint(double d) {
  if (d == std::floor(d))
    return d;
  return std::floor(d + 0.5);
}

/* Computes the truncated integer. */
inline double aint(double d) {
  return (d >= 0.0) ? std::floor(d) : std::ceil(d);
}

/* These are provided to give consistent 
   interface for double with double-double and quad-double. */
inline void sincosh(double t, double &sinh_t, double &cosh_t) {
  sinh_t = std::sinh(t);
  cosh_t = std::cosh(t);
}

inline double sqr(double t) {
  return t * t;
}

inline double to_double(double a) { return a; }
inline int    to_int(double a) { return static_cast<int>(a); }

}

#endif /* _QD_INLINE_H */
