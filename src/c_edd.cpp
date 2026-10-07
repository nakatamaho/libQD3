/*
 * src/c_edd.cpp
 *
 * C wrapper functions for edd_real.
 */
#include <cstdlib>
#include <cstring>

#include "config.h"
#include <qd/edd_real.h>
#include <qd/c_edd.h>

namespace {

inline edd_real edd_from_c(const _Float64x *a) {
  return edd_real(static_cast<edd_word>(a[0]), static_cast<edd_word>(a[1]));
}

inline void edd_to_c(const edd_real &a, _Float64x *ptr) {
  ptr[0] = static_cast<_Float64x>(a[0]);
  ptr[1] = static_cast<_Float64x>(a[1]);
}

} // namespace

extern "C" {

void c_edd_copy(const _Float64x *a, _Float64x *b) {
  b[0] = a[0];
  b[1] = a[1];
}

void c_edd_copy_d(double a, _Float64x *b) {
  b[0] = static_cast<_Float64x>((edd_word) a);
  b[1] = static_cast<_Float64x>((edd_word) 0.0);
}

void c_edd_add(const _Float64x *a, const _Float64x *b, _Float64x *c) {
  edd_real cc = edd_from_c(a) + edd_from_c(b);
  edd_to_c(cc, c);
}

void c_edd_sub(const _Float64x *a, const _Float64x *b, _Float64x *c) {
  edd_real cc = edd_from_c(a) - edd_from_c(b);
  edd_to_c(cc, c);
}

void c_edd_mul(const _Float64x *a, const _Float64x *b, _Float64x *c) {
  edd_real cc = edd_from_c(a) * edd_from_c(b);
  edd_to_c(cc, c);
}

void c_edd_div(const _Float64x *a, const _Float64x *b, _Float64x *c) {
  edd_real cc = edd_from_c(a) / edd_from_c(b);
  edd_to_c(cc, c);
}

void c_edd_sqrt(const _Float64x *a, _Float64x *b) {
  edd_real bb = sqrt(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_sqr(const _Float64x *a, _Float64x *b) {
  edd_real bb = sqr(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_abs(const _Float64x *a, _Float64x *b) {
  edd_real bb = abs(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_exp(const _Float64x *a, _Float64x *b) {
  edd_real bb = exp(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_log(const _Float64x *a, _Float64x *b) {
  edd_real bb = log(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_log10(const _Float64x *a, _Float64x *b) {
  edd_real bb = log10(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_sin(const _Float64x *a, _Float64x *b) {
  edd_real bb = sin(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_cos(const _Float64x *a, _Float64x *b) {
  edd_real bb = cos(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_tan(const _Float64x *a, _Float64x *b) {
  edd_real bb = tan(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_atan2(const _Float64x *a, const _Float64x *b, _Float64x *c) {
  edd_real cc = atan2(edd_from_c(a), edd_from_c(b));
  edd_to_c(cc, c);
}

void c_edd_sinh(const _Float64x *a, _Float64x *b) {
  edd_real bb = sinh(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_cosh(const _Float64x *a, _Float64x *b) {
  edd_real bb = cosh(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_tanh(const _Float64x *a, _Float64x *b) {
  edd_real bb = tanh(edd_from_c(a));
  edd_to_c(bb, b);
}

void c_edd_read(const char *s, _Float64x *a) {
  edd_real aa(s);
  edd_to_c(aa, a);
}

void c_edd_swrite(const _Float64x *a, int precision, char *s, int len) {
  edd_from_c(a).write(s, len, precision);
}

void c_edd_neg(const _Float64x *a, _Float64x *b) {
  edd_real bb = -edd_from_c(a);
  edd_to_c(bb, b);
}

void c_edd_comp(const _Float64x *a, const _Float64x *b, int *result) {
  edd_real aa = edd_from_c(a);
  edd_real bb = edd_from_c(b);
  *result = (aa < bb) ? -1 : ((aa > bb) ? 1 : 0);
}

void c_edd_comp_edd_d(const _Float64x *a, double b, int *result) {
  edd_real aa = edd_from_c(a);
  edd_real bb((edd_word) b);
  *result = (aa < bb) ? -1 : ((aa > bb) ? 1 : 0);
}

void c_edd_pi(_Float64x *a) {
  edd_to_c(edd_real::_pi, a);
}

void c_edd_2pi(_Float64x *a) {
  edd_to_c(edd_real::_2pi, a);
}

_Float64x c_edd_epsilon(void) {
  return static_cast<_Float64x>(edd_real::_eps);
}


void c_edd_add_edd_d(const _Float64x *a, double b, _Float64x *c) {
  edd_to_c(edd_from_c(a) + edd_real((edd_word) b), c);
}

void c_edd_add_d_edd(double a, const _Float64x *b, _Float64x *c) {
  edd_to_c(edd_real((edd_word) a) + edd_from_c(b), c);
}

void c_edd_sub_edd_d(const _Float64x *a, double b, _Float64x *c) {
  edd_to_c(edd_from_c(a) - edd_real((edd_word) b), c);
}

void c_edd_sub_d_edd(double a, const _Float64x *b, _Float64x *c) {
  edd_to_c(edd_real((edd_word) a) - edd_from_c(b), c);
}

void c_edd_mul_edd_d(const _Float64x *a, double b, _Float64x *c) {
  edd_to_c(edd_from_c(a) * edd_real((edd_word) b), c);
}

void c_edd_mul_d_edd(double a, const _Float64x *b, _Float64x *c) {
  edd_to_c(edd_real((edd_word) a) * edd_from_c(b), c);
}

void c_edd_div_edd_d(const _Float64x *a, double b, _Float64x *c) {
  edd_to_c(edd_from_c(a) / edd_real((edd_word) b), c);
}

void c_edd_div_d_edd(double a, const _Float64x *b, _Float64x *c) {
  edd_to_c(edd_real((edd_word) a) / edd_from_c(b), c);
}

void c_edd_selfadd(const _Float64x *a, _Float64x *b) {
  edd_to_c(edd_from_c(b) + edd_from_c(a), b);
}

void c_edd_selfadd_d(double a, _Float64x *b) {
  edd_to_c(edd_from_c(b) + edd_real((edd_word) a), b);
}

void c_edd_selfsub(const _Float64x *a, _Float64x *b) {
  edd_to_c(edd_from_c(b) - edd_from_c(a), b);
}

void c_edd_selfsub_d(double a, _Float64x *b) {
  edd_to_c(edd_from_c(b) - edd_real((edd_word) a), b);
}

void c_edd_selfmul(const _Float64x *a, _Float64x *b) {
  edd_to_c(edd_from_c(b) * edd_from_c(a), b);
}

void c_edd_selfmul_d(double a, _Float64x *b) {
  edd_to_c(edd_from_c(b) * edd_real((edd_word) a), b);
}

void c_edd_selfdiv(const _Float64x *a, _Float64x *b) {
  edd_to_c(edd_from_c(b) / edd_from_c(a), b);
}

void c_edd_selfdiv_d(double a, _Float64x *b) {
  edd_to_c(edd_from_c(b) / edd_real((edd_word) a), b);
}

void c_edd_npwr(const _Float64x *a, int n, _Float64x *b) {
  edd_to_c(npwr(edd_from_c(a), n), b);
}

void c_edd_nroot(const _Float64x *a, int n, _Float64x *b) {
  edd_to_c(nroot(edd_from_c(a), n), b);
}

void c_edd_nint(const _Float64x *a, _Float64x *b) {
  edd_to_c(nint(edd_from_c(a)), b);
}

void c_edd_aint(const _Float64x *a, _Float64x *b) {
  edd_to_c(aint(edd_from_c(a)), b);
}

void c_edd_floor(const _Float64x *a, _Float64x *b) {
  edd_to_c(floor(edd_from_c(a)), b);
}

void c_edd_ceil(const _Float64x *a, _Float64x *b) {
  edd_to_c(ceil(edd_from_c(a)), b);
}

void c_edd_asin(const _Float64x *a, _Float64x *b) {
  edd_to_c(asin(edd_from_c(a)), b);
}

void c_edd_acos(const _Float64x *a, _Float64x *b) {
  edd_to_c(acos(edd_from_c(a)), b);
}

void c_edd_atan(const _Float64x *a, _Float64x *b) {
  edd_to_c(atan(edd_from_c(a)), b);
}

void c_edd_asinh(const _Float64x *a, _Float64x *b) {
  edd_to_c(asinh(edd_from_c(a)), b);
}

void c_edd_acosh(const _Float64x *a, _Float64x *b) {
  edd_to_c(acosh(edd_from_c(a)), b);
}

void c_edd_atanh(const _Float64x *a, _Float64x *b) {
  edd_to_c(atanh(edd_from_c(a)), b);
}

void c_edd_sincos(const _Float64x *a, _Float64x *s, _Float64x *c) {
  edd_real ss, cc;
  sincos(edd_from_c(a), ss, cc);
  edd_to_c(ss, s);
  edd_to_c(cc, c);
}

void c_edd_sincosh(const _Float64x *a, _Float64x *s, _Float64x *c) {
  edd_real ss, cc;
  sincosh(edd_from_c(a), ss, cc);
  edd_to_c(ss, s);
  edd_to_c(cc, c);
}

void c_edd_comp_d_edd(double a, const _Float64x *b, int *result) {
  edd_real aa((edd_word) a);
  edd_real bb = edd_from_c(b);
  *result = (aa < bb) ? -1 : ((aa > bb) ? 1 : 0);
}

/* Uniform in [0, 1): 31 random bits per step, as ddrand/qdrand do. */
void c_edd_rand(_Float64x *a) {
  static const double m_const = 4.6566128730773926e-10; /* 2^-31 */
  double m = m_const;
  edd_real r((edd_word) 0.0);
  for (int i = 0; i < 5; i++, m *= m_const) {
    r += edd_real((edd_word) (std::rand() * m));
  }
  edd_to_c(r, a);
}

} // extern "C"
