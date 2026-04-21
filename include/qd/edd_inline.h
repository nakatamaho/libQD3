/*
 * include/edd_inline.h
 *
 * Inline building blocks for edd_real based on native _Float64x limbs.
 */
#ifndef _QD_EDD_INLINE_H
#define _QD_EDD_INLINE_H

#include <qd/edd_real.h>

#ifndef QD_INLINE
#define inline
#endif

namespace edd {

inline edd_word fabsx(edd_word x) { return __builtin_fabsf64x(x); }
inline edd_word sqrtx(edd_word x) { return __builtin_sqrtf64x(x); }
inline edd_word floorx(edd_word x) { return __builtin_floorf64x(x); }
inline edd_word ldexpx(edd_word x, int e) { return __builtin_ldexpf64x(x, e); }
inline edd_word copysignx(edd_word x, edd_word y) {
  return __builtin_copysignf64x(x, y);
}
inline edd_word log10x(edd_word x) { return __builtin_log10f64x(x); }

inline bool isnanx(edd_word x) { return __builtin_isnan(x); }
inline bool isinfx(edd_word x) { return __builtin_isinf_sign(x) != 0; }
inline bool isfinitex(edd_word x) { return __builtin_isfinite(x); }

inline edd_word d_nan() { return __builtin_nanf64x(""); }
inline edd_word d_inf() { return __builtin_huge_valf64x(); }

inline edd_word splitter() {
  return ldexpx((edd_word) 1.0, QD_EDD_SPLIT_BITS) + (edd_word) 1.0;
}

inline edd_word split_thresh() {
  return ldexpx((edd_word) 1.0, QD_EDD_FLT64X_MAX_EXP - QD_EDD_SPLIT_SCALE_BITS);
}

inline edd_word div_rescale_thresh() {
  return ldexpx((edd_word) 1.0,
      QD_EDD_FLT64X_MAX_EXP - QD_EDD_FLT64X_MANT_DIG);
}

inline edd_word sqrt_rescale_thresh() {
  return ldexpx((edd_word) 1.0, QD_EDD_FLT64X_MAX_EXP - 3);
}

inline edd_word quick_two_sum(edd_word a, edd_word b, edd_word &err) {
  edd_word s = a + b;
  err = b - (s - a);
  return s;
}

inline edd_word quick_two_diff(edd_word a, edd_word b, edd_word &err) {
  edd_word s = a - b;
  err = (a - s) - b;
  return s;
}

inline edd_word two_sum(edd_word a, edd_word b, edd_word &err) {
  edd_word s = a + b;
  edd_word bb = s - a;
  err = (a - (s - bb)) + (b - bb);
  return s;
}

inline edd_word two_diff(edd_word a, edd_word b, edd_word &err) {
  edd_word s = a - b;
  edd_word bb = s - a;
  err = (a - (s - bb)) - (b + bb);
  return s;
}

inline void split(edd_word a, edd_word &hi, edd_word &lo) {
  edd_word temp;
  if (fabsx(a) > split_thresh()) {
    a = ldexpx(a, -QD_EDD_SPLIT_SCALE_BITS);
    temp = splitter() * a;
    hi = temp - (temp - a);
    lo = a - hi;
    hi = ldexpx(hi, QD_EDD_SPLIT_SCALE_BITS);
    lo = ldexpx(lo, QD_EDD_SPLIT_SCALE_BITS);
  } else {
    temp = splitter() * a;
    hi = temp - (temp - a);
    lo = a - hi;
  }
}

inline edd_word two_prod(edd_word a, edd_word b, edd_word &err) {
  edd_word a_hi, a_lo, b_hi, b_lo;
  edd_word p = a * b;
  split(a, a_hi, a_lo);
  split(b, b_hi, b_lo);
  err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
  return p;
}

inline edd_word two_sqr(edd_word a, edd_word &err) {
  edd_word hi, lo;
  edd_word q = a * a;
  split(a, hi, lo);
  err = ((hi * hi - q) + ((edd_word) 2.0) * hi * lo) + lo * lo;
  return q;
}

inline edd_word nint(edd_word x) {
  if (x == floorx(x))
    return x;
  return floorx(x + (edd_word) 0.5);
}

inline void renorm2(edd_word &c0, edd_word &c1) {
  if (isinfx(c0))
    return;
  c0 = quick_two_sum(c0, c1, c1);
}

inline void renorm4(edd_word &c0, edd_word &c1, edd_word &c2, edd_word &c3) {
  edd_word s0, s1, s2 = (edd_word) 0.0;

  if (isinfx(c0))
    return;

  s0 = quick_two_sum(c2, c3, c3);
  s0 = quick_two_sum(c1, s0, c2);
  c0 = quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;
  if (s1 != (edd_word) 0.0) {
    s1 = quick_two_sum(s1, c2, s2);
    if (s2 != (edd_word) 0.0)
      s1 = quick_two_sum(s1, s2 + c3, s2);
    else
      s1 = quick_two_sum(s1, c3, s2);
  } else {
    s0 = quick_two_sum(s0, c2, s1);
    if (s1 != (edd_word) 0.0)
      s1 = quick_two_sum(s1, c3, s2);
    else
      s0 = quick_two_sum(s0, c3, s1);
  }

  c0 = s0;
  c1 = s1;
}

inline edd_real from_qd_truncate(const qd_real &a) {
  edd_word c0 = (edd_word) a[0];
  edd_word c1 = (edd_word) a[1];
  edd_word c2 = (edd_word) a[2];
  edd_word c3 = (edd_word) a[3];
  renorm4(c0, c1, c2, c3);
  return edd_real(c0, c1);
}

inline void word_to_doubles(edd_word x, double &hi, double &lo) {
  hi = static_cast<double>(x);
  lo = static_cast<double>(x - (edd_word) hi);
}

inline qd_real to_qd_conversion(const edd_real &a) {
  double x0, x1, x2, x3;
  word_to_doubles(a[0], x0, x1);
  word_to_doubles(a[1], x2, x3);
  qd_real q(x0, x1, x2, x3);
  q.renorm();
  return q;
}

template <class A, class B>
inline bool comparison_eq(const A &a, const B &b) {
  edd_real ea(a);
  edd_real eb(b);
  if (ea.isnan() || eb.isnan())
    return false;
  if (ea.isinf() || eb.isinf())
    return ea[0] == eb[0];
  edd_real d = ea - eb;
  return d.is_zero();
}

template <class A, class B>
inline bool comparison_lt(const A &a, const B &b) {
  edd_real ea(a);
  edd_real eb(b);
  if (ea.isnan() || eb.isnan())
    return false;
  if (ea.isinf() || eb.isinf())
    return ea[0] < eb[0];
  edd_real d = ea - eb;
  return d.is_negative();
}

template <class A, class B>
inline bool comparison_gt(const A &a, const B &b) {
  edd_real ea(a);
  edd_real eb(b);
  if (ea.isnan() || eb.isnan())
    return false;
  if (ea.isinf() || eb.isinf())
    return ea[0] > eb[0];
  edd_real d = ea - eb;
  return d.is_positive();
}

template <class A, class B>
inline bool comparison_le(const A &a, const B &b) {
  edd_real ea(a);
  edd_real eb(b);
  if (ea.isnan() || eb.isnan())
    return false;
  if (ea.isinf() || eb.isinf())
    return ea[0] <= eb[0];
  edd_real d = ea - eb;
  return d.is_zero() || d.is_negative();
}

template <class A, class B>
inline bool comparison_ge(const A &a, const B &b) {
  edd_real ea(a);
  edd_real eb(b);
  if (ea.isnan() || eb.isnan())
    return false;
  if (ea.isinf() || eb.isinf())
    return ea[0] >= eb[0];
  edd_real d = ea - eb;
  return d.is_zero() || d.is_positive();
}

} // namespace edd

inline bool edd_real::isnan() const {
  return edd::isnanx(x[0]) || edd::isnanx(x[1]);
}

inline bool edd_real::isfinite() const {
  return edd::isfinitex(x[0]);
}

inline bool edd_real::isinf() const {
  return edd::isinfx(x[0]);
}

inline edd_real edd_real::add(edd_word a, edd_word b) {
  edd_word s, e;
  s = edd::two_sum(a, b, e);
  return edd_real(s, e);
}

inline edd_real::edd_real(const dd_real &dd) {
  x[0] = (edd_word) dd._hi();
  x[1] = (edd_word) dd._lo();
  edd::renorm2(x[0], x[1]);
}

inline edd_real::edd_real(const qd_real &qd) {
  *this = edd::from_qd_truncate(qd);
}

inline edd_real operator+(const edd_real &a, edd_word b) {
  edd_word s1, s2;
  s1 = edd::two_sum(a.x[0], b, s2);
  s2 += a.x[1];
  s1 = edd::quick_two_sum(s1, s2, s2);
  return edd_real(s1, s2);
}

inline edd_real operator+(edd_word a, const edd_real &b) {
  return b + a;
}

inline edd_real operator+(const edd_real &a, const edd_real &b) {
  edd_word s1, s2, t1, t2;
  s1 = edd::two_sum(a.x[0], b.x[0], s2);
  t1 = edd::two_sum(a.x[1], b.x[1], t2);
  s2 += t1;
  s1 = edd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = edd::quick_two_sum(s1, s2, s2);
  return edd_real(s1, s2);
}

inline edd_real &edd_real::operator+=(edd_word a) {
  edd_word s1, s2;
  s1 = edd::two_sum(x[0], a, s2);
  s2 += x[1];
  x[0] = edd::quick_two_sum(s1, s2, x[1]);
  return *this;
}

inline edd_real &edd_real::operator+=(const edd_real &a) {
  edd_word s1, s2, t1, t2;
  s1 = edd::two_sum(x[0], a.x[0], s2);
  t1 = edd::two_sum(x[1], a.x[1], t2);
  s2 += t1;
  s1 = edd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  x[0] = edd::quick_two_sum(s1, s2, x[1]);
  return *this;
}

inline edd_real edd_real::sub(edd_word a, edd_word b) {
  edd_word s, e;
  s = edd::two_diff(a, b, e);
  return edd_real(s, e);
}

inline edd_real operator-(const edd_real &a, edd_word b) {
  edd_word s1, s2;
  s1 = edd::two_diff(a.x[0], b, s2);
  s2 += a.x[1];
  s1 = edd::quick_two_sum(s1, s2, s2);
  return edd_real(s1, s2);
}

inline edd_real operator-(edd_word a, const edd_real &b) {
  edd_word s1, s2;
  s1 = edd::two_diff(a, b.x[0], s2);
  s2 -= b.x[1];
  s1 = edd::quick_two_sum(s1, s2, s2);
  return edd_real(s1, s2);
}

inline edd_real operator-(const edd_real &a, const edd_real &b) {
  edd_word s1, s2, t1, t2;
  s1 = edd::two_diff(a.x[0], b.x[0], s2);
  t1 = edd::two_diff(a.x[1], b.x[1], t2);
  s2 += t1;
  s1 = edd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = edd::quick_two_sum(s1, s2, s2);
  return edd_real(s1, s2);
}

inline edd_real &edd_real::operator-=(edd_word a) {
  edd_word s1, s2;
  s1 = edd::two_diff(x[0], a, s2);
  s2 += x[1];
  x[0] = edd::quick_two_sum(s1, s2, x[1]);
  return *this;
}

inline edd_real &edd_real::operator-=(const edd_real &a) {
  edd_word s1, s2, t1, t2;
  s1 = edd::two_diff(x[0], a.x[0], s2);
  t1 = edd::two_diff(x[1], a.x[1], t2);
  s2 += t1;
  s1 = edd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  x[0] = edd::quick_two_sum(s1, s2, x[1]);
  return *this;
}

inline edd_real edd_real::operator-() const {
  return edd_real(-x[0], -x[1]);
}

inline edd_real edd_real::mul(edd_word a, edd_word b) {
  edd_word p, e;
  p = edd::two_prod(a, b, e);
  return edd_real(p, e);
}

inline edd_real ldexp(const edd_real &a, int exp) {
  return edd_real(edd::ldexpx(a.x[0], exp), edd::ldexpx(a.x[1], exp));
}

inline edd_real mul_pwr2(const edd_real &a, edd_word b) {
  return edd_real(a.x[0] * b, a.x[1] * b);
}

inline edd_real operator*(const edd_real &a, edd_word b) {
  edd_word p1, p2;
  p1 = edd::two_prod(a.x[0], b, p2);
  p2 += a.x[1] * b;
  p1 = edd::quick_two_sum(p1, p2, p2);
  return edd_real(p1, p2);
}

inline edd_real operator*(edd_word a, const edd_real &b) {
  return b * a;
}

inline edd_real operator*(const edd_real &a, const edd_real &b) {
  edd_word p1, p2;
  p1 = edd::two_prod(a.x[0], b.x[0], p2);
  p2 += a.x[0] * b.x[1] + a.x[1] * b.x[0];
  p1 = edd::quick_two_sum(p1, p2, p2);
  return edd_real(p1, p2);
}

inline edd_real &edd_real::operator*=(edd_word a) {
  edd_word p1, p2;
  p1 = edd::two_prod(x[0], a, p2);
  p2 += x[1] * a;
  x[0] = edd::quick_two_sum(p1, p2, x[1]);
  return *this;
}

inline edd_real &edd_real::operator*=(const edd_real &a) {
  edd_word p1, p2;
  p1 = edd::two_prod(x[0], a.x[0], p2);
  p2 += x[0] * a.x[1] + x[1] * a.x[0];
  x[0] = edd::quick_two_sum(p1, p2, x[1]);
  return *this;
}

inline bool edd_real_div_needs_rescale(edd_word a_hi) {
  return edd::fabsx(a_hi) > edd::div_rescale_thresh();
}

inline edd_real edd_real::div(edd_word a, edd_word b) {
  edd_word q1, q2;
  edd_word p1, p2;
  edd_word s, e;

  const bool rescale = edd_real_div_needs_rescale(a);
  const edd_word aa = rescale ? edd::ldexpx(a, -QD_EDD_FLT64X_MANT_DIG) : a;

  q1 = aa / b;
  p1 = edd::two_prod(q1, b, p2);
  s = edd::two_diff(aa, p1, e);
  e -= p2;
  q2 = (s + e) / b;
  s = edd::quick_two_sum(q1, q2, e);

  edd_real r(s, e);
  return rescale ? mul_pwr2(r, edd::ldexpx((edd_word) 1.0, QD_EDD_FLT64X_MANT_DIG)) : r;
}

inline edd_real operator/(const edd_real &a, edd_word b) {
  edd_word q1, q2;
  edd_word p1, p2;
  edd_word s, e;
  edd_real r;

  const bool rescale = edd_real_div_needs_rescale(a.x[0]);
  const edd_real aa = rescale ? mul_pwr2(a, edd::ldexpx((edd_word) 1.0, -QD_EDD_FLT64X_MANT_DIG)) : a;

  q1 = aa.x[0] / b;
  p1 = edd::two_prod(q1, b, p2);
  s = edd::two_diff(aa.x[0], p1, e);
  e += aa.x[1];
  e -= p2;
  q2 = (s + e) / b;
  r.x[0] = edd::quick_two_sum(q1, q2, r.x[1]);

  return rescale ? mul_pwr2(r, edd::ldexpx((edd_word) 1.0, QD_EDD_FLT64X_MANT_DIG)) : r;
}

inline edd_real operator/(const edd_real &a, const edd_real &b) {
  edd_word q1, q2, q3;
  edd_real r;

  const bool rescale = edd_real_div_needs_rescale(a.x[0]);
  const edd_real aa = rescale ? mul_pwr2(a, edd::ldexpx((edd_word) 1.0, -QD_EDD_FLT64X_MANT_DIG)) : a;

  q1 = aa.x[0] / b.x[0];
  r = aa - q1 * b;
  q2 = r.x[0] / b.x[0];
  r -= q2 * b;
  q3 = r.x[0] / b.x[0];

  q1 = edd::quick_two_sum(q1, q2, q2);
  r = edd_real(q1, q2) + q3;
  return rescale ? mul_pwr2(r, edd::ldexpx((edd_word) 1.0, QD_EDD_FLT64X_MANT_DIG)) : r;
}

inline edd_real operator/(edd_word a, const edd_real &b) {
  return edd_real(a) / b;
}

inline edd_real &edd_real::operator/=(edd_word a) {
  *this = *this / a;
  return *this;
}

inline edd_real &edd_real::operator/=(const edd_real &a) {
  *this = *this / a;
  return *this;
}

inline edd_real sqr(const edd_real &a) {
  edd_word p1, p2, s2;
  p1 = edd::two_sqr(a.x[0], p2);
  p2 += ((edd_word) 2.0) * a.x[0] * a.x[1];
  p2 += a.x[1] * a.x[1];
  p1 = edd::quick_two_sum(p1, p2, s2);
  return edd_real(p1, s2);
}

inline edd_real edd_real::sqr(edd_word a) {
  edd_word p1, p2;
  p1 = edd::two_sqr(a, p2);
  return edd_real(p1, p2);
}

inline edd_real edd_real::operator^(int n) const {
  return npwr(*this, n);
}

inline edd_real &edd_real::operator=(edd_word a) {
  x[0] = a;
  x[1] = (edd_word) 0.0;
  return *this;
}

inline edd_real &edd_real::operator=(double a) {
  x[0] = (edd_word) a;
  x[1] = (edd_word) 0.0;
  return *this;
}

inline edd_real &edd_real::operator=(int a) {
  x[0] = (edd_word) a;
  x[1] = (edd_word) 0.0;
  return *this;
}

inline edd_real &edd_real::operator=(const dd_real &a) {
  x[0] = (edd_word) a._hi();
  x[1] = (edd_word) a._lo();
  edd::renorm2(x[0], x[1]);
  return *this;
}

inline edd_real &edd_real::operator=(const qd_real &a) {
  *this = edd::from_qd_truncate(a);
  return *this;
}

inline bool operator==(const edd_real &a, edd_word b) {
  return edd::comparison_eq(a, b);
}

inline bool operator==(edd_word a, const edd_real &b) {
  return edd::comparison_eq(a, b);
}

inline bool operator==(const edd_real &a, const edd_real &b) {
  return edd::comparison_eq(a, b);
}

inline bool operator!=(const edd_real &a, edd_word b) {
  return !(a == b);
}

inline bool operator!=(edd_word a, const edd_real &b) {
  return !(a == b);
}

inline bool operator!=(const edd_real &a, const edd_real &b) {
  return !(a == b);
}

inline bool operator<(const edd_real &a, edd_word b) {
  return edd::comparison_lt(a, b);
}

inline bool operator<(edd_word a, const edd_real &b) {
  return edd::comparison_lt(a, b);
}

inline bool operator<(const edd_real &a, const edd_real &b) {
  return edd::comparison_lt(a, b);
}

inline bool operator>(const edd_real &a, edd_word b) {
  return edd::comparison_gt(a, b);
}

inline bool operator>(edd_word a, const edd_real &b) {
  return edd::comparison_gt(a, b);
}

inline bool operator>(const edd_real &a, const edd_real &b) {
  return edd::comparison_gt(a, b);
}

inline bool operator<=(const edd_real &a, edd_word b) {
  return edd::comparison_le(a, b);
}

inline bool operator<=(edd_word a, const edd_real &b) {
  return edd::comparison_le(a, b);
}

inline bool operator<=(const edd_real &a, const edd_real &b) {
  return edd::comparison_le(a, b);
}

inline bool operator>=(const edd_real &a, edd_word b) {
  return edd::comparison_ge(a, b);
}

inline bool operator>=(edd_word a, const edd_real &b) {
  return edd::comparison_ge(a, b);
}

inline bool operator>=(const edd_real &a, const edd_real &b) {
  return edd::comparison_ge(a, b);
}

inline edd_real abs(const edd_real &a) {
  return (a[0] < (edd_word) 0.0) ? -a : a;
}

inline edd_real fabs(const edd_real &a) {
  return abs(a);
}

inline bool edd_real::is_zero() const {
  return x[0] == (edd_word) 0.0 && x[1] == (edd_word) 0.0;
}

inline bool edd_real::is_one() const {
  return x[0] == (edd_word) 1.0 && x[1] == (edd_word) 0.0;
}

inline bool edd_real::is_positive() const {
  return x[0] > (edd_word) 0.0 ||
      (x[0] == (edd_word) 0.0 && x[1] > (edd_word) 0.0);
}

inline bool edd_real::is_negative() const {
  return x[0] < (edd_word) 0.0 ||
      (x[0] == (edd_word) 0.0 && x[1] < (edd_word) 0.0);
}

inline edd_word to_float64x(const edd_real &a) {
  return a[0] + a[1];
}

inline dd_real to_dd_real(const edd_real &a) {
  qd_real q = edd::to_qd_conversion(a);
  return dd_real(q[0], q[1]);
}

inline edd_real to_edd_real(const dd_real &a) {
  return edd_real(a);
}

inline edd_real to_edd_real(const qd_real &a) {
  return edd_real(a);
}

inline qd_real to_qd_real(const edd_real &a) {
  return edd::to_qd_conversion(a);
}

inline double to_double(const edd_real &a) {
  return static_cast<double>(a[0] + a[1]);
}

inline int to_int(const edd_real &a) {
  return static_cast<int>(a[0] + a[1]);
}

#endif /* _QD_EDD_INLINE_H */
