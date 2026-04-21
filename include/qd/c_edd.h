/*
 * include/c_edd.h
 *
 * C wrapper function prototypes for edd_real.
 */
#ifndef _QD_C_EDD_H
#define _QD_C_EDD_H

#include <qd/qd_config.h>
#include <qd/fpu.h>

#ifndef QD_HAVE_EDD_REAL
#error "c_edd.h requires QD_HAVE_EDD_REAL."
#endif

#ifdef __cplusplus
extern "C" {
#endif

void c_edd_copy(const _Float64x *a, _Float64x *b);
void c_edd_copy_d(double a, _Float64x *b);

void c_edd_add(const _Float64x *a, const _Float64x *b, _Float64x *c);
void c_edd_sub(const _Float64x *a, const _Float64x *b, _Float64x *c);
void c_edd_mul(const _Float64x *a, const _Float64x *b, _Float64x *c);
void c_edd_div(const _Float64x *a, const _Float64x *b, _Float64x *c);

void c_edd_sqrt(const _Float64x *a, _Float64x *b);
void c_edd_sqr(const _Float64x *a, _Float64x *b);
void c_edd_abs(const _Float64x *a, _Float64x *b);

void c_edd_exp(const _Float64x *a, _Float64x *b);
void c_edd_log(const _Float64x *a, _Float64x *b);
void c_edd_log10(const _Float64x *a, _Float64x *b);

void c_edd_sin(const _Float64x *a, _Float64x *b);
void c_edd_cos(const _Float64x *a, _Float64x *b);
void c_edd_tan(const _Float64x *a, _Float64x *b);
void c_edd_atan2(const _Float64x *a, const _Float64x *b, _Float64x *c);

void c_edd_sinh(const _Float64x *a, _Float64x *b);
void c_edd_cosh(const _Float64x *a, _Float64x *b);
void c_edd_tanh(const _Float64x *a, _Float64x *b);

void c_edd_read(const char *s, _Float64x *a);
void c_edd_swrite(const _Float64x *a, int precision, char *s, int len);
void c_edd_neg(const _Float64x *a, _Float64x *b);
void c_edd_comp(const _Float64x *a, const _Float64x *b, int *result);
void c_edd_comp_edd_d(const _Float64x *a, double b, int *result);
void c_edd_pi(_Float64x *a);
void c_edd_2pi(_Float64x *a);
_Float64x c_edd_epsilon(void);

#ifdef __cplusplus
}
#endif

#endif /* _QD_C_EDD_H */
