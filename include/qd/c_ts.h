/*
 * include/c_ts.h
 *
 * C wrapper function prototypes for triple-single precision arithmetic.
 * A ts_real value is passed as an array of 3 floats; double scalar
 * operands are converted exactly where the ts_real format allows.
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
#ifndef _QD_C_TS_H
#define _QD_C_TS_H

#include <qd/qd_config.h>

#ifdef __cplusplus
extern "C" {
#endif

/* add */
void c_ts_add(const float *a, const float *b, float *c);
void c_ts_add_ts_ds(const float *a, const float *b, float *c);
void c_ts_add_ds_ts(const float *a, const float *b, float *c);
void c_ts_add_ts_qs(const float *a, const float *b, float *c);
void c_ts_add_qs_ts(const float *a, const float *b, float *c);
void c_ts_add_d_ts(double a, const float *b, float *c);
void c_ts_add_ts_d(const float *a, double b, float *c);

/* sub */
void c_ts_sub(const float *a, const float *b, float *c);
void c_ts_sub_ts_ds(const float *a, const float *b, float *c);
void c_ts_sub_ds_ts(const float *a, const float *b, float *c);
void c_ts_sub_ts_qs(const float *a, const float *b, float *c);
void c_ts_sub_qs_ts(const float *a, const float *b, float *c);
void c_ts_sub_d_ts(double a, const float *b, float *c);
void c_ts_sub_ts_d(const float *a, double b, float *c);

/* mul */
void c_ts_mul(const float *a, const float *b, float *c);
void c_ts_mul_ts_ds(const float *a, const float *b, float *c);
void c_ts_mul_ds_ts(const float *a, const float *b, float *c);
void c_ts_mul_ts_qs(const float *a, const float *b, float *c);
void c_ts_mul_qs_ts(const float *a, const float *b, float *c);
void c_ts_mul_d_ts(double a, const float *b, float *c);
void c_ts_mul_ts_d(const float *a, double b, float *c);

/* div */
void c_ts_div(const float *a, const float *b, float *c);
void c_ts_div_ts_ds(const float *a, const float *b, float *c);
void c_ts_div_ds_ts(const float *a, const float *b, float *c);
void c_ts_div_ts_qs(const float *a, const float *b, float *c);
void c_ts_div_qs_ts(const float *a, const float *b, float *c);
void c_ts_div_d_ts(double a, const float *b, float *c);
void c_ts_div_ts_d(const float *a, double b, float *c);

/* copy */
void c_ts_copy(const float *a, float *b);
void c_ts_copy_ds(const float *a, float *b);
void c_ts_copy_qs(const float *a, float *b);
void c_ts_copy_d(double a, float *b);

void c_ts_selfadd(const float *a, float *b);
void c_ts_selfadd_ds(const float *a, float *b);
void c_ts_selfadd_qs(const float *a, float *b);
void c_ts_selfadd_d(double a, float *b);
void c_ts_selfsub(const float *a, float *b);
void c_ts_selfsub_ds(const float *a, float *b);
void c_ts_selfsub_qs(const float *a, float *b);
void c_ts_selfsub_d(double a, float *b);
void c_ts_selfmul(const float *a, float *b);
void c_ts_selfmul_ds(const float *a, float *b);
void c_ts_selfmul_qs(const float *a, float *b);
void c_ts_selfmul_d(double a, float *b);
void c_ts_selfdiv(const float *a, float *b);
void c_ts_selfdiv_ds(const float *a, float *b);
void c_ts_selfdiv_qs(const float *a, float *b);
void c_ts_selfdiv_d(double a, float *b);

void c_ts_sqrt(const float *a, float *b);
void c_ts_sqr(const float *a, float *b);
void c_ts_abs(const float *a, float *b);
void c_ts_nint(const float *a, float *b);
void c_ts_aint(const float *a, float *b);
void c_ts_floor(const float *a, float *b);
void c_ts_ceil(const float *a, float *b);
void c_ts_exp(const float *a, float *b);
void c_ts_log(const float *a, float *b);
void c_ts_log10(const float *a, float *b);
void c_ts_sin(const float *a, float *b);
void c_ts_cos(const float *a, float *b);
void c_ts_tan(const float *a, float *b);
void c_ts_asin(const float *a, float *b);
void c_ts_acos(const float *a, float *b);
void c_ts_atan(const float *a, float *b);
void c_ts_sinh(const float *a, float *b);
void c_ts_cosh(const float *a, float *b);
void c_ts_tanh(const float *a, float *b);
void c_ts_asinh(const float *a, float *b);
void c_ts_acosh(const float *a, float *b);
void c_ts_atanh(const float *a, float *b);
void c_ts_npwr(const float *a, int n, float *b);
void c_ts_nroot(const float *a, int n, float *b);
void c_ts_atan2(const float *a, const float *b, float *c);
void c_ts_sincos(const float *a, float *s, float *c);
void c_ts_sincosh(const float *a, float *s, float *c);

void c_ts_read(const char *s, float *a);
void c_ts_swrite(const float *a, int precision, char *s, int len);
void c_ts_write(const float *a);
void c_ts_neg(const float *a, float *b);
void c_ts_rand(float *a);
void c_ts_comp(const float *a, const float *b, int *result);
void c_ts_comp_ts_d(const float *a, double b, int *result);
void c_ts_comp_d_ts(double a, const float *b, int *result);
void c_ts_comp_ts_ds(const float *a, const float *b, int *result);
void c_ts_comp_ds_ts(const float *a, const float *b, int *result);
void c_ts_comp_ts_qs(const float *a, const float *b, int *result);
void c_ts_comp_qs_ts(const float *a, const float *b, int *result);
void c_ts_pi(float *a);
void c_ts_2pi(float *a);
double c_ts_epsilon(void);

#ifdef __cplusplus
}
#endif

#endif /* _QD_C_TS_H */
