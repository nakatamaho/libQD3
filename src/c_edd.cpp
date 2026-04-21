/*
 * src/c_edd.cpp
 *
 * C wrapper functions for edd_real.
 */
#include <cstring>

#include "config.h"
#include <qd/edd_real.h>
#include <qd/c_edd.h>

#define TO_EDD_PTR(a, ptr) do { (ptr)[0] = (a)[0]; (ptr)[1] = (a)[1]; } while (0)

extern "C" {

void c_edd_copy(const _Float64x *a, _Float64x *b) {
  b[0] = a[0];
  b[1] = a[1];
}

void c_edd_copy_d(double a, _Float64x *b) {
  b[0] = (edd_word) a;
  b[1] = (edd_word) 0.0;
}

void c_edd_add(const _Float64x *a, const _Float64x *b, _Float64x *c) {
  edd_real cc = edd_real(a) + edd_real(b);
  TO_EDD_PTR(cc, c);
}

void c_edd_sub(const _Float64x *a, const _Float64x *b, _Float64x *c) {
  edd_real cc = edd_real(a) - edd_real(b);
  TO_EDD_PTR(cc, c);
}

void c_edd_mul(const _Float64x *a, const _Float64x *b, _Float64x *c) {
  edd_real cc = edd_real(a) * edd_real(b);
  TO_EDD_PTR(cc, c);
}

void c_edd_div(const _Float64x *a, const _Float64x *b, _Float64x *c) {
  edd_real cc = edd_real(a) / edd_real(b);
  TO_EDD_PTR(cc, c);
}

void c_edd_sqrt(const _Float64x *a, _Float64x *b) {
  edd_real bb = sqrt(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_sqr(const _Float64x *a, _Float64x *b) {
  edd_real bb = sqr(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_abs(const _Float64x *a, _Float64x *b) {
  edd_real bb = abs(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_exp(const _Float64x *a, _Float64x *b) {
  edd_real bb = exp(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_log(const _Float64x *a, _Float64x *b) {
  edd_real bb = log(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_log10(const _Float64x *a, _Float64x *b) {
  edd_real bb = log10(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_sin(const _Float64x *a, _Float64x *b) {
  edd_real bb = sin(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_cos(const _Float64x *a, _Float64x *b) {
  edd_real bb = cos(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_tan(const _Float64x *a, _Float64x *b) {
  edd_real bb = tan(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_atan2(const _Float64x *a, const _Float64x *b, _Float64x *c) {
  edd_real cc = atan2(edd_real(a), edd_real(b));
  TO_EDD_PTR(cc, c);
}

void c_edd_sinh(const _Float64x *a, _Float64x *b) {
  edd_real bb = sinh(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_cosh(const _Float64x *a, _Float64x *b) {
  edd_real bb = cosh(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_tanh(const _Float64x *a, _Float64x *b) {
  edd_real bb = tanh(edd_real(a));
  TO_EDD_PTR(bb, b);
}

void c_edd_read(const char *s, _Float64x *a) {
  edd_real aa(s);
  TO_EDD_PTR(aa, a);
}

void c_edd_swrite(const _Float64x *a, int precision, char *s, int len) {
  edd_real(a).write(s, len, precision);
}

void c_edd_neg(const _Float64x *a, _Float64x *b) {
  edd_real bb = -edd_real(a);
  TO_EDD_PTR(bb, b);
}

void c_edd_comp(const _Float64x *a, const _Float64x *b, int *result) {
  edd_real aa(a);
  edd_real bb(b);
  *result = (aa < bb) ? -1 : ((aa > bb) ? 1 : 0);
}

void c_edd_comp_edd_d(const _Float64x *a, double b, int *result) {
  edd_real aa(a);
  edd_real bb((edd_word) b);
  *result = (aa < bb) ? -1 : ((aa > bb) ? 1 : 0);
}

void c_edd_pi(_Float64x *a) {
  TO_EDD_PTR(edd_real::_pi, a);
}

void c_edd_2pi(_Float64x *a) {
  TO_EDD_PTR(edd_real::_2pi, a);
}

_Float64x c_edd_epsilon(void) {
  return edd_real::_eps;
}

} // extern "C"
