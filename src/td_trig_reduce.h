/*
 * src/td_trig_reduce.h
 *
 * Private triple-double trigonometric argument reduction helpers.
 */
#ifndef QD_TD_TRIG_REDUCE_H
#define QD_TD_TRIG_REDUCE_H

#include <cmath>

#include <qd/td_real.h>
#include <qd/inline.h>

#ifndef QD_INLINE
#include <qd/td_inline.h>
#endif

namespace qd_internal {

static const td_real td_trig_pi16(1.963495408493620697e-01,
                                  7.654042494670957545e-18,
                                  -1.871731131073962291e-34);
static const double td_trig_2pi_limb3 = 2.224908441726730563e-49;
static const double td_trig_pi2_limb3 = 5.562271104316826408e-50;
static const double td_trig_pi16_limb3 = 6.952838880396033010e-51;

struct TdTrigExpansion4 {
  double x[4];
};

inline td_real td_trig_nint(const td_real &a) {
  double x0 = qd::nint(a[0]);
  double x1 = 0.0;
  double x2 = 0.0;

  if (x0 == a[0]) {
    x1 = qd::nint(a[1]);
    if (x1 == a[1]) {
      x2 = qd::nint(a[2]);
    } else if (std::abs(x1 - a[1]) == 0.5 && a[2] < 0.0) {
      x1 -= 1.0;
    }
  } else if (std::abs(x0 - a[0]) == 0.5 && a[1] < 0.0) {
    x0 -= 1.0;
  }

  td::renorm(x0, x1, x2);
  return td_real(x0, x1, x2);
}

inline TdTrigExpansion4 td_trig_make4(double x0, double x1, double x2,
                                      double x3) {
  TdTrigExpansion4 r;
  r.x[0] = x0;
  r.x[1] = x1;
  r.x[2] = x2;
  r.x[3] = x3;
  return r;
}

inline void td_trig_renorm4(double &c0, double &c1, double &c2, double &c3) {
  double s0, s1, s2 = 0.0, s3 = 0.0;

  if (QD_ISINF(c0)) return;

  s0 = qd::quick_two_sum(c2, c3, c3);
  s0 = qd::quick_two_sum(c1, s0, c2);
  c0 = qd::quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;
  if (s1 != 0.0) {
    s1 = qd::quick_two_sum(s1, c2, s2);
    if (s2 != 0.0)
      s2 = qd::quick_two_sum(s2, c3, s3);
    else
      s1 = qd::quick_two_sum(s1, c3, s2);
  } else {
    s0 = qd::quick_two_sum(s0, c2, s1);
    if (s1 != 0.0)
      s1 = qd::quick_two_sum(s1, c3, s2);
    else
      s0 = qd::quick_two_sum(s0, c3, s1);
  }

  c0 = s0;
  c1 = s1;
  c2 = s2;
  c3 = s3;
}

inline void td_trig_renorm5(double &c0, double &c1, double &c2, double &c3,
                            double &c4) {
  double s0, s1, s2 = 0.0, s3 = 0.0;

  if (QD_ISINF(c0)) return;

  s0 = qd::quick_two_sum(c3, c4, c4);
  s0 = qd::quick_two_sum(c2, s0, c3);
  s0 = qd::quick_two_sum(c1, s0, c2);
  c0 = qd::quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;

  if (s1 != 0.0) {
    s1 = qd::quick_two_sum(s1, c2, s2);
    if (s2 != 0.0) {
      s2 = qd::quick_two_sum(s2, c3, s3);
      if (s3 != 0.0)
        s3 += c4;
      else
        s2 = qd::quick_two_sum(s2, c4, s3);
    } else {
      s1 = qd::quick_two_sum(s1, c3, s2);
      if (s2 != 0.0)
        s2 = qd::quick_two_sum(s2, c4, s3);
      else
        s1 = qd::quick_two_sum(s1, c4, s2);
    }
  } else {
    s0 = qd::quick_two_sum(s0, c2, s1);
    if (s1 != 0.0) {
      s1 = qd::quick_two_sum(s1, c3, s2);
      if (s2 != 0.0)
        s2 = qd::quick_two_sum(s2, c4, s3);
      else
        s1 = qd::quick_two_sum(s1, c4, s2);
    } else {
      s0 = qd::quick_two_sum(s0, c3, s1);
      if (s1 != 0.0)
        s1 = qd::quick_two_sum(s1, c4, s2);
      else
        s0 = qd::quick_two_sum(s0, c4, s1);
    }
  }

  c0 = s0;
  c1 = s1;
  c2 = s2;
  c3 = s3;
}

inline TdTrigExpansion4 td_trig_from_td(const td_real &a) {
  return td_trig_make4(a[0], a[1], a[2], 0.0);
}

inline td_real td_trig_to_td(const TdTrigExpansion4 &a) {
  double x0 = a.x[0];
  double x1 = a.x[1];
  double x2 = a.x[2];
  double x3 = a.x[3];
  td_trig_renorm4(x0, x1, x2, x3);
  return td_real(x0, x1, x2);
}

inline TdTrigExpansion4 td_trig_neg(const TdTrigExpansion4 &a) {
  return td_trig_make4(-a.x[0], -a.x[1], -a.x[2], -a.x[3]);
}

inline TdTrigExpansion4 td_trig_add4(const TdTrigExpansion4 &a,
                                     const TdTrigExpansion4 &b) {
  int i = 0;
  int j = 0;
  int k = 0;
  double s;
  double t;
  double u;
  double v;
  double x[4] = {0.0, 0.0, 0.0, 0.0};

  if (std::abs(a.x[i]) > std::abs(b.x[j])) u = a.x[i++];
  else u = b.x[j++];

  if (std::abs(a.x[i]) > std::abs(b.x[j])) v = a.x[i++];
  else v = b.x[j++];

  u = qd::quick_two_sum(u, v, v);

  while (k < 4) {
    if (i >= 4 && j >= 4) {
      x[k] = u;
      if (k < 3) x[++k] = v;
      break;
    }

    if (i >= 4) t = b.x[j++];
    else if (j >= 4) t = a.x[i++];
    else if (std::abs(a.x[i]) > std::abs(b.x[j])) t = a.x[i++];
    else t = b.x[j++];

    s = td::quick_three_accum(u, v, t);
    if (s != 0.0) x[k++] = s;
  }

  for (k = i; k < 4; k++) x[3] += a.x[k];
  for (k = j; k < 4; k++) x[3] += b.x[k];

  td_trig_renorm4(x[0], x[1], x[2], x[3]);
  return td_trig_make4(x[0], x[1], x[2], x[3]);
}

inline TdTrigExpansion4 td_trig_mul4(const TdTrigExpansion4 &a,
                                     const TdTrigExpansion4 &b) {
  double p0, p1, p2, p3, p4, p5;
  double q0, q1, q2, q3, q4, q5;
  double p6, p7, p8, p9;
  double q6, q7, q8, q9;
  double r0, r1;
  double t0, t1;
  double s0, s1, s2;

  p0 = qd::two_prod(a.x[0], b.x[0], q0);

  p1 = qd::two_prod(a.x[0], b.x[1], q1);
  p2 = qd::two_prod(a.x[1], b.x[0], q2);

  p3 = qd::two_prod(a.x[0], b.x[2], q3);
  p4 = qd::two_prod(a.x[1], b.x[1], q4);
  p5 = qd::two_prod(a.x[2], b.x[0], q5);

  td::three_sum(p1, p2, q0);

  td::three_sum(p2, q1, q2);
  td::three_sum(p3, p4, p5);
  s0 = qd::two_sum(p2, p3, t0);
  s1 = qd::two_sum(q1, p4, t1);
  s2 = q2 + p5;
  s1 = qd::two_sum(s1, t0, t0);
  s2 += (t0 + t1);

  p6 = qd::two_prod(a.x[0], b.x[3], q6);
  p7 = qd::two_prod(a.x[1], b.x[2], q7);
  p8 = qd::two_prod(a.x[2], b.x[1], q8);
  p9 = qd::two_prod(a.x[3], b.x[0], q9);

  q0 = qd::two_sum(q0, q3, q3);
  q4 = qd::two_sum(q4, q5, q5);
  p6 = qd::two_sum(p6, p7, p7);
  p8 = qd::two_sum(p8, p9, p9);
  t0 = qd::two_sum(q0, q4, t1);
  t1 += (q3 + q5);
  r0 = qd::two_sum(p6, p8, r1);
  r1 += (p7 + p9);
  q3 = qd::two_sum(t0, r0, q4);
  q4 += (t1 + r1);
  t0 = qd::two_sum(q3, s1, t1);
  t1 += q4;

  t1 += a.x[1] * b.x[3] + a.x[2] * b.x[2] + a.x[3] * b.x[1] +
        q6 + q7 + q8 + q9 + s2;

  td_trig_renorm5(p0, p1, s0, t0, t1);
  return td_trig_make4(p0, p1, s0, t0);
}

inline TdTrigExpansion4 td_trig_sub_mul_const(const TdTrigExpansion4 &a,
                                              const td_real &n, double c0,
                                              double c1, double c2,
                                              double c3) {
  const TdTrigExpansion4 product =
      td_trig_mul4(td_trig_from_td(n), td_trig_make4(c0, c1, c2, c3));
  return td_trig_add4(a, td_trig_neg(product));
}

inline void reduce_td_trig_arg(const td_real &a, td_real &t, int &j, int &k) {
  td_real z = td_trig_nint(a / td_real::_2pi);
  TdTrigExpansion4 r = td_trig_sub_mul_const(
      td_trig_from_td(a), z, td_real::_2pi[0], td_real::_2pi[1],
      td_real::_2pi[2], td_trig_2pi_limb3);

  td_real q = td_trig_nint(td_trig_to_td(r) / td_real::_pi2);
  TdTrigExpansion4 tq = td_trig_sub_mul_const(
      r, q, td_real::_pi2[0], td_real::_pi2[1], td_real::_pi2[2],
      td_trig_pi2_limb3);
  j = static_cast<int>(q[0]);
  while (j > 2) j -= 4;
  while (j < -2) j += 4;

  q = td_trig_nint(td_trig_to_td(tq) / td_trig_pi16);
  tq = td_trig_sub_mul_const(tq, q, td_trig_pi16[0], td_trig_pi16[1],
                             td_trig_pi16[2], td_trig_pi16_limb3);
  k = static_cast<int>(q[0]);
  t = td_trig_to_td(tq);
}

}  // namespace qd_internal

#endif
