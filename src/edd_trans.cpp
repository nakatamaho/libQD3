/*
 * src/edd_trans.cpp
 *
 * Transcendental functions for edd_real.
 */
#include <algorithm>
#include <cmath>

#include "config.h"
#include <qd/edd_real.h>

#ifndef QD_INLINE
#include <qd/edd_inline.h>
#endif

namespace {

inline edd_real signed_zero(edd_word sign_source) {
  return edd_real(edd::copysignx((edd_word) 0.0, sign_source));
}

inline edd_word exp_overflow_limit() {
  return to_float64x(edd_real::_log2) * (edd_word) 16384.0;
}

inline edd_word exp_underflow_limit() {
  return to_float64x(edd_real::_log2) * (edd_word) -16445.0;
}

class qd_fpu_scope {
public:
  qd_fpu_scope() : old_cw_(0) { fpu_fix_start(&old_cw_); }
  ~qd_fpu_scope() { fpu_fix_end(&old_cw_); }

private:
  qd_fpu_scope(const qd_fpu_scope &);
  qd_fpu_scope &operator=(const qd_fpu_scope &);

  unsigned int old_cw_;
};

inline qd_real to_qd_mode(const edd_real &a) {
  qd_fpu_scope scope;
  return to_qd_real(a);
}

inline edd_real edd_nint_internal(const edd_real &a) {
  edd_word x0 = edd::nint(a[0]);
  edd_word x1 = (edd_word) 0.0;

  if (x0 == a[0]) {
    x1 = edd::nint(a[1]);
    if (x1 != a[1] && edd::fabsx(x1 - a[1]) == (edd_word) 0.5 && a[1] < (edd_word) 0.0) {
      x1 -= (edd_word) 1.0;
    }
  } else if (edd::fabsx(x0 - a[0]) == (edd_word) 0.5 && a[1] < (edd_word) 0.0) {
    x0 -= (edd_word) 1.0;
  }

  edd::renorm2(x0, x1);
  return edd_real(x0, x1);
}

static const int edd_exp_squares = 12;
static const edd_word edd_exp_k = (edd_word) 4096.0;
static const edd_word edd_exp_inv_k = (edd_word) 1.0 / edd_exp_k;
static const edd_real edd_pi16 = to_edd_real(qd_real(
    1.963495408493620697e-01, 7.654042494670957545e-18,
    -1.871731131073962291e-34, 8.553725411101285102e-51));

static const qd_real edd_sin_table_qd[4] = {
  qd_real(1.950903220161282758e-01, -7.991079068461731263e-18,
      6.184627002422071324e-34, -3.243618621526994199e-50),
  qd_real(3.826834323650897818e-01, -1.005077269646158761e-17,
      -2.060531630280669539e-34, 1.250254541928171821e-50),
  qd_real(5.555702330196021776e-01, 4.709410940561676821e-17,
      -2.064052038368292118e-33, -8.964987880773293235e-50),
  qd_real(7.071067811865475727e-01, -4.833646656726456726e-17,
      2.069337654349706787e-33, -2.467773495734175784e-50)
};

static const qd_real edd_cos_table_qd[4] = {
  qd_real(9.807852804032304306e-01, 1.854693999782500573e-17,
      -1.069656444553075966e-33, 4.044946249568903222e-50),
  qd_real(9.238795325112867385e-01, 1.764504708433667706e-17,
      -5.044253732158682026e-34, 1.503235349992942096e-50),
  qd_real(8.314696123025452357e-01, 1.407385698472802389e-18,
      4.695131538398083001e-35, -2.729711263611673679e-51),
  qd_real(7.071067811865475727e-01, -4.833646656726456726e-17,
      2.069337654349706787e-33, -2.467773495734175784e-50)
};

inline edd_real sin_table(int i) {
  return to_edd_real(edd_sin_table_qd[i]);
}

inline edd_real cos_table(int i) {
  return to_edd_real(edd_cos_table_qd[i]);
}

inline edd_real expm1_small(const edd_real &a) {
  edd_real term = a;
  edd_real sum = a;
  edd_word n = (edd_word) 1.0;
  const edd_real thresh(abs(a) * (edd_real::_eps * edd_exp_inv_k));

  do {
    n += (edd_word) 1.0;
    term *= a;
    term /= n;
    sum += term;
  } while (abs(term) > thresh);

  return sum;
}

inline edd_real sin_taylor(const edd_real &a) {
  if (a.is_zero())
    return signed_zero(a[0]);

  const edd_real thresh = abs(a) * (edd_real::_eps * (edd_word) 0.5);
  edd_real x = -sqr(a);
  edd_real term = a;
  edd_real sum = a;
  edd_word n = (edd_word) 1.0;

  do {
    term *= x;
    term /= (n + (edd_word) 1.0) * (n + (edd_word) 2.0);
    sum += term;
    n += (edd_word) 2.0;
  } while (abs(term) > thresh);

  return sum;
}

inline edd_real cos_taylor(const edd_real &a) {
  if (a.is_zero())
    return edd_real((edd_word) 1.0);

  const edd_real thresh(edd_real::_eps * (edd_word) 0.5);
  edd_real x = -sqr(a);
  edd_real term((edd_word) 1.0);
  edd_real sum((edd_word) 1.0);
  edd_word n = (edd_word) 0.0;

  do {
    term *= x;
    term /= (n + (edd_word) 1.0) * (n + (edd_word) 2.0);
    sum += term;
    n += (edd_word) 2.0;
  } while (abs(term) > thresh);

  return sum;
}

inline void sincos_taylor(const edd_real &a, edd_real &sin_a, edd_real &cos_a) {
  if (a.is_zero()) {
    sin_a = signed_zero(a[0]);
    cos_a = (edd_word) 1.0;
    return;
  }

  sin_a = sin_taylor(a);
  cos_a = cos_taylor(a);
}

inline void reduce_edd_trig_arg(const edd_real &a, edd_real &t, int &j, int &k) {
  const qd_real a_qd = to_qd_mode(a);

  if (!std::isfinite(to_double(a_qd))) {
    edd_real z = edd_nint_internal(a / edd_real::_2pi);
    edd_real r = a - edd_real::_2pi * z;

    edd_word q = edd::floorx(r[0] / edd_real::_pi2[0] + (edd_word) 0.5);
    t = r - edd_real::_pi2 * q;
    j = static_cast<int>(q);
    while (j > 2) j -= 4;
    while (j < -2) j += 4;

    q = edd::floorx(t[0] / edd_pi16[0] + (edd_word) 0.5);
    t -= edd_pi16 * q;
    k = static_cast<int>(q);
    return;
  }

  qd_real t_qd;
  {
    qd_fpu_scope scope;
    qd_real z = nint(a_qd / qd_real::_2pi);
    qd_real r = a_qd - qd_real::_2pi * z;

    double q = std::floor(r[0] / qd_real::_pi2[0] + 0.5);
    t_qd = r - qd_real::_pi2 * q;
    j = static_cast<int>(q);
    while (j > 2) j -= 4;
    while (j < -2) j += 4;

    q = std::floor(t_qd[0] / edd_pi16[0] + 0.5);
    t_qd -= to_qd_real(edd_pi16) * q;
    k = static_cast<int>(q);
  }
  t = to_edd_real(t_qd);
}

inline void sincos_native(const edd_real &a, edd_real &s, edd_real &c) {
  edd_real t;
  int j;
  int k;
  reduce_edd_trig_arg(a, t, j, k);

  int abs_j = std::abs(j);
  int abs_k = std::abs(k);
  if (abs_k > 4) {
    edd_real::error("(edd_real::sincos): Cannot reduce modulo pi/16.");
    s = c = edd_real::_nan;
    return;
  }

  edd_real sin_t;
  edd_real cos_t;
  edd_real ss;
  edd_real cc;
  sincos_taylor(t, sin_t, cos_t);

  if (abs_k == 0) {
    ss = sin_t;
    cc = cos_t;
  } else {
    edd_real u = cos_table(abs_k - 1);
    edd_real v = sin_table(abs_k - 1);

    if (k > 0) {
      ss = u * sin_t + v * cos_t;
      cc = u * cos_t - v * sin_t;
    } else {
      ss = u * sin_t - v * cos_t;
      cc = u * cos_t + v * sin_t;
    }
  }

  if (abs_j == 0) {
    s = ss;
    c = cc;
  } else if (j == 1) {
    s = cc;
    c = -ss;
  } else if (j == -1) {
    s = -cc;
    c = ss;
  } else {
    s = -ss;
    c = -cc;
  }
}

} // namespace

edd_real nroot(const edd_real &a, int n) {
  if (n <= 0) {
    edd_real::error("(edd_real::nroot): N must be positive.");
    return edd_real::_nan;
  }

  if (n % 2 == 0 && a.is_negative()) {
    edd_real::error("(edd_real::nroot): Negative argument.");
    return edd_real::_nan;
  }
  if (n == 1)
    return a;
  if (n == 2)
    return sqrt(a);
  if (a.is_zero())
    return edd_real((edd_word) 0.0);

  edd_real r = abs(a);
  edd_real x = edd_real(edd::expx(-edd::logx(r[0]) / (edd_word) n));
  edd_word inv_n = (edd_word) 1.0 / (edd_word) n;

  x += x * ((edd_word) 1.0 - r * npwr(x, n)) * inv_n;
  x += x * ((edd_word) 1.0 - r * npwr(x, n)) * inv_n;
  x += x * ((edd_word) 1.0 - r * npwr(x, n)) * inv_n;

  if (a[0] < (edd_word) 0.0)
    x = -x;

  return (edd_word) 1.0 / x;
}

edd_real pow(const edd_real &a, const edd_real &b) {
  return exp(b * log(a));
}

edd_real exp(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (a.is_zero())
    return (edd_word) 1.0;
  if (a.isinf())
    return a.is_positive() ? edd_real::_inf : edd_real((edd_word) 0.0);
  if (a[0] >= exp_overflow_limit())
    return edd_real::_inf;
  if (a[0] <= exp_underflow_limit())
    return edd_real((edd_word) 0.0);
  if (a.is_one())
    return edd_real::_e;

  edd_word m = edd::floorx(a[0] / edd_real::_log2[0] + (edd_word) 0.5);
  edd_real r = mul_pwr2(a - edd_real::_log2 * m, edd_exp_inv_k);
  edd_real s = expm1_small(r);

  for (int i = 0; i < edd_exp_squares; i++) {
    s = mul_pwr2(s, (edd_word) 2.0) + sqr(s);
  }

  return ldexp(s + (edd_word) 1.0, static_cast<int>(m));
}

edd_real log(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (a.is_zero())
    return -edd_real::_inf;
  if (a.is_negative()) {
    edd_real::error("(edd_real::log): Non-positive argument.");
    return edd_real::_nan;
  }
  if (a.isinf())
    return edd_real::_inf;
  if (a.is_one())
    return edd_real((edd_word) 0.0);

  int e;
  edd_word m = edd::frexpx(a[0], &e);
  edd_real x = edd_real(edd::logx(m)) + edd_real::_log2 * (edd_word) e;

  x = x + a * exp(-x) - (edd_word) 1.0;
  x = x + a * exp(-x) - (edd_word) 1.0;
  x = x + a * exp(-x) - (edd_word) 1.0;
  return x;
}

edd_real log10(const edd_real &a) {
  return log(a) / edd_real::_log10;
}

edd_real log2(const edd_real &a) {
  return log(a) / edd_real::_log2;
}

edd_real exp2(const edd_real &a) {
  return exp(edd_real::_log2 * a);
}

edd_real expm1(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (a.is_zero())
    return a;
  if (abs(a) < (edd_word) 0.125) {
    edd_real term = a;
    edd_real sum = a;
    edd_word n = (edd_word) 1.0;
    edd_real thresh = abs(a) * edd_real::_eps;
    if (thresh.is_zero())
      thresh = edd_real(edd_real::_eps);
    do {
      n += (edd_word) 1.0;
      term *= a;
      term /= n;
      sum += term;
    } while (abs(term) > thresh);
    return sum;
  }
  return exp(a) - (edd_word) 1.0;
}

edd_real log1p(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (a == (edd_word) -1.0)
    return -edd_real::_inf;
  if (a < (edd_word) -1.0) {
    edd_real::error("(edd_real::log1p): Argument out of domain.");
    return edd_real::_nan;
  }
  if (a.is_zero())
    return a;
  if (abs(a) < (edd_word) 0.125) {
    edd_real term = a;
    edd_real sum = a;
    edd_word n = (edd_word) 1.0;
    edd_real thresh = abs(a) * edd_real::_eps;
    if (thresh.is_zero())
      thresh = edd_real(edd_real::_eps);
    do {
      n += (edd_word) 1.0;
      term *= -a;
      edd_real add = term / n;
      sum += add;
      if (abs(add) <= thresh)
        break;
    } while (true);
    return sum;
  }
  return log((edd_word) 1.0 + a);
}

void sincos(const edd_real &a, edd_real &s, edd_real &c) {
  if (a.isnan()) {
    s = c = edd_real::_nan;
    return;
  }
  if (a.isinf()) {
    edd_real::error("(edd_real::sincos): Infinite argument.");
    s = c = edd_real::_nan;
    return;
  }
  if (a.is_zero()) {
    s = signed_zero(a[0]);
    c = (edd_word) 1.0;
    return;
  }

  sincos_native(a, s, c);
}

edd_real sin(const edd_real &a) {
  edd_real s;
  edd_real c;
  sincos(a, s, c);
  return s;
}

edd_real cos(const edd_real &a) {
  edd_real s;
  edd_real c;
  sincos(a, s, c);
  return c;
}

edd_real tan(const edd_real &a) {
  edd_real s;
  edd_real c;
  sincos(a, s, c);
  return s / c;
}

edd_real asin(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (abs(a) > (edd_word) 1.0) {
    edd_real::error("(edd_real::asin): Argument out of domain.");
    return edd_real::_nan;
  }
  return atan2(a, sqrt((edd_word) 1.0 - sqr(a)));
}

edd_real acos(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (abs(a) > (edd_word) 1.0) {
    edd_real::error("(edd_real::acos): Argument out of domain.");
    return edd_real::_nan;
  }
  return atan2(sqrt((edd_word) 1.0 - sqr(a)), a);
}

edd_real atan(const edd_real &a) {
  return atan2(a, (edd_word) 1.0);
}

edd_real atan2(const edd_real &y, const edd_real &x) {
  if (x.isnan() || y.isnan())
    return edd_real::_nan;
  if (x.is_zero()) {
    if (y.is_zero()) {
      edd_real::error("(edd_real::atan2): Both arguments zero.");
      return edd_real::_nan;
    }
    return edd::signbitx(y[0]) ? -edd_real::_pi2 : edd_real::_pi2;
  }

  if (y.is_zero()) {
    if (x.is_positive())
      return signed_zero(y[0]);
    return edd::signbitx(y[0]) ? -edd_real::_pi : edd_real::_pi;
  }

  if (x == y)
    return y.is_positive() ? edd_real::_pi4 : -edd_real::_3pi4;
  if (x == -y)
    return y.is_positive() ? edd_real::_3pi4 : -edd_real::_pi4;

  edd_real r = sqrt(sqr(x) + sqr(y));
  edd_real xx = x / r;
  edd_real yy = y / r;
  edd_real z(edd_real(edd::atan2x(to_float64x(y), to_float64x(x))));
  edd_real sin_z;
  edd_real cos_z;

  if (edd::fabsx(xx[0]) > edd::fabsx(yy[0])) {
    for (int i = 0; i < 3; i++) {
      sincos(z, sin_z, cos_z);
      z += (yy - sin_z) / cos_z;
    }
  } else {
    for (int i = 0; i < 3; i++) {
      sincos(z, sin_z, cos_z);
      z -= (xx - cos_z) / sin_z;
    }
  }

  return z;
}

edd_real sinh(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (a.is_zero())
    return signed_zero(a[0]);
  if (a.isinf())
    return a;

  if (abs(a) > (edd_word) 0.05) {
    edd_real ea = exp(a);
    return mul_pwr2(ea - ((edd_word) 1.0 / ea), (edd_word) 0.5);
  }

  edd_real s = a;
  edd_real t = a;
  edd_real r = sqr(a);
  edd_word m = (edd_word) 1.0;
  edd_real thresh = abs(a) * edd_real::_eps;

  do {
    m += (edd_word) 2.0;
    t *= r;
    t /= (m - (edd_word) 1.0) * m;
    s += t;
  } while (abs(t) > thresh);

  return s;
}

edd_real cosh(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (a.is_zero())
    return (edd_word) 1.0;
  if (a.isinf())
    return edd_real::_inf;

  edd_real ea = exp(a);
  return mul_pwr2(ea + ((edd_word) 1.0 / ea), (edd_word) 0.5);
}

edd_real tanh(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (a.is_zero())
    return signed_zero(a[0]);
  if (a.isinf())
    return a.is_positive() ? edd_real((edd_word) 1.0) : edd_real((edd_word) -1.0);

  if (abs(a) > (edd_word) 0.05) {
    edd_real ea = exp(a);
    edd_real inv_ea = (edd_word) 1.0 / ea;
    return (ea - inv_ea) / (ea + inv_ea);
  }

  edd_real s = sinh(a);
  edd_real c = sqrt((edd_word) 1.0 + sqr(s));
  return s / c;
}

void sincosh(const edd_real &a, edd_real &s, edd_real &c) {
  if (abs(a) <= (edd_word) 0.05) {
    s = sinh(a);
    c = sqrt((edd_word) 1.0 + sqr(s));
  } else {
    edd_real ea = exp(a);
    edd_real inv_ea = (edd_word) 1.0 / ea;
    s = mul_pwr2(ea - inv_ea, (edd_word) 0.5);
    c = mul_pwr2(ea + inv_ea, (edd_word) 0.5);
  }
}

edd_real asinh(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (a.isinf())
    return a;
  return log(a + sqrt(sqr(a) + (edd_word) 1.0));
}

edd_real acosh(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (a < (edd_word) 1.0) {
    edd_real::error("(edd_real::acosh): Argument out of domain.");
    return edd_real::_nan;
  }
  if (a.isinf())
    return edd_real::_inf;
  return log(a + sqrt(sqr(a) - (edd_word) 1.0));
}

edd_real atanh(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (abs(a) >= (edd_word) 1.0) {
    edd_real::error("(edd_real::atanh): Argument out of domain.");
    return edd_real::_nan;
  }
  return mul_pwr2(log(((edd_word) 1.0 + a) / ((edd_word) 1.0 - a)),
                  (edd_word) 0.5);
}

edd_real nint(const edd_real &a) {
  return edd_nint_internal(a);
}

edd_real fmod(const edd_real &a, const edd_real &b) {
  edd_real n = aint(a / b);
  return a - b * n;
}

edd_real hypot(const edd_real &a, const edd_real &b) {
  if (a.isnan() || b.isnan())
    return edd_real::_nan;
  if (a.isinf() || b.isinf())
    return edd_real::_inf;
  edd_real x = abs(a);
  edd_real y = abs(b);
  if (x < y)
    std::swap(x, y);
  if (x.is_zero())
    return edd_real((edd_word) 0.0);
  edd_real r = y / x;
  return x * sqrt((edd_word) 1.0 + sqr(r));
}

edd_real cbrt(const edd_real &a) {
  if (a.isnan())
    return edd_real::_nan;
  if (a.is_zero())
    return a;
  return nroot(a, 3);
}

edd_real trunc(const edd_real &a) {
  return aint(a);
}

edd_real round(const edd_real &a) {
  if (a.isnan() || a.isinf() || a.is_zero())
    return a;
  return a.is_positive() ? floor(a + (edd_word) 0.5)
                         : ceil(a - (edd_word) 0.5);
}
