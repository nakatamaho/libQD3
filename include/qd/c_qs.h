/*
 * include/c_qs.h
 *
 * C wrapper function prototypes for quad-single precision arithmetic.
 * A qs_real value is passed as an array of 4 floats; double scalar
 * operands are converted exactly where the qs_real format allows.
 *
 *
 * Copyright (c) 2026, Nakata Maho
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef _QD_C_QS_H
#define _QD_C_QS_H

#include <qd/qd_config.h>

#ifdef __cplusplus
extern "C" {
#endif

/* add */
void c_qs_add(const float *a, const float *b, float *c);
void c_qs_add_qs_ds(const float *a, const float *b, float *c);
void c_qs_add_ds_qs(const float *a, const float *b, float *c);
void c_qs_add_qs_ts(const float *a, const float *b, float *c);
void c_qs_add_ts_qs(const float *a, const float *b, float *c);
void c_qs_add_d_qs(double a, const float *b, float *c);
void c_qs_add_qs_d(const float *a, double b, float *c);

/* sub */
void c_qs_sub(const float *a, const float *b, float *c);
void c_qs_sub_qs_ds(const float *a, const float *b, float *c);
void c_qs_sub_ds_qs(const float *a, const float *b, float *c);
void c_qs_sub_qs_ts(const float *a, const float *b, float *c);
void c_qs_sub_ts_qs(const float *a, const float *b, float *c);
void c_qs_sub_d_qs(double a, const float *b, float *c);
void c_qs_sub_qs_d(const float *a, double b, float *c);

/* mul */
void c_qs_mul(const float *a, const float *b, float *c);
void c_qs_mul_qs_ds(const float *a, const float *b, float *c);
void c_qs_mul_ds_qs(const float *a, const float *b, float *c);
void c_qs_mul_qs_ts(const float *a, const float *b, float *c);
void c_qs_mul_ts_qs(const float *a, const float *b, float *c);
void c_qs_mul_d_qs(double a, const float *b, float *c);
void c_qs_mul_qs_d(const float *a, double b, float *c);

/* div */
void c_qs_div(const float *a, const float *b, float *c);
void c_qs_div_qs_ds(const float *a, const float *b, float *c);
void c_qs_div_ds_qs(const float *a, const float *b, float *c);
void c_qs_div_qs_ts(const float *a, const float *b, float *c);
void c_qs_div_ts_qs(const float *a, const float *b, float *c);
void c_qs_div_d_qs(double a, const float *b, float *c);
void c_qs_div_qs_d(const float *a, double b, float *c);

/* copy */
void c_qs_copy(const float *a, float *b);
void c_qs_copy_ds(const float *a, float *b);
void c_qs_copy_ts(const float *a, float *b);
void c_qs_copy_d(double a, float *b);

void c_qs_selfadd(const float *a, float *b);
void c_qs_selfadd_ds(const float *a, float *b);
void c_qs_selfadd_ts(const float *a, float *b);
void c_qs_selfadd_d(double a, float *b);
void c_qs_selfsub(const float *a, float *b);
void c_qs_selfsub_ds(const float *a, float *b);
void c_qs_selfsub_ts(const float *a, float *b);
void c_qs_selfsub_d(double a, float *b);
void c_qs_selfmul(const float *a, float *b);
void c_qs_selfmul_ds(const float *a, float *b);
void c_qs_selfmul_ts(const float *a, float *b);
void c_qs_selfmul_d(double a, float *b);
void c_qs_selfdiv(const float *a, float *b);
void c_qs_selfdiv_ds(const float *a, float *b);
void c_qs_selfdiv_ts(const float *a, float *b);
void c_qs_selfdiv_d(double a, float *b);

void c_qs_sqrt(const float *a, float *b);
void c_qs_sqr(const float *a, float *b);
void c_qs_abs(const float *a, float *b);
void c_qs_nint(const float *a, float *b);
void c_qs_aint(const float *a, float *b);
void c_qs_floor(const float *a, float *b);
void c_qs_ceil(const float *a, float *b);
void c_qs_exp(const float *a, float *b);
void c_qs_log(const float *a, float *b);
void c_qs_log10(const float *a, float *b);
void c_qs_sin(const float *a, float *b);
void c_qs_cos(const float *a, float *b);
void c_qs_tan(const float *a, float *b);
void c_qs_asin(const float *a, float *b);
void c_qs_acos(const float *a, float *b);
void c_qs_atan(const float *a, float *b);
void c_qs_sinh(const float *a, float *b);
void c_qs_cosh(const float *a, float *b);
void c_qs_tanh(const float *a, float *b);
void c_qs_asinh(const float *a, float *b);
void c_qs_acosh(const float *a, float *b);
void c_qs_atanh(const float *a, float *b);
void c_qs_npwr(const float *a, int n, float *b);
void c_qs_nroot(const float *a, int n, float *b);
void c_qs_atan2(const float *a, const float *b, float *c);
void c_qs_sincos(const float *a, float *s, float *c);
void c_qs_sincosh(const float *a, float *s, float *c);

void c_qs_read(const char *s, float *a);
void c_qs_swrite(const float *a, int precision, char *s, int len);
void c_qs_write(const float *a);
void c_qs_neg(const float *a, float *b);
void c_qs_rand(float *a);
void c_qs_comp(const float *a, const float *b, int *result);
void c_qs_comp_qs_d(const float *a, double b, int *result);
void c_qs_comp_d_qs(double a, const float *b, int *result);
void c_qs_comp_qs_ds(const float *a, const float *b, int *result);
void c_qs_comp_ds_qs(const float *a, const float *b, int *result);
void c_qs_comp_qs_ts(const float *a, const float *b, int *result);
void c_qs_comp_ts_qs(const float *a, const float *b, int *result);
void c_qs_pi(float *a);
void c_qs_2pi(float *a);
double c_qs_epsilon(void);

#ifdef __cplusplus
}
#endif

#endif /* _QD_C_QS_H */
