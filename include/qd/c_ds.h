/*
 * include/c_ds.h
 *
 * C wrapper function prototypes for double-single precision arithmetic.
 * A ds_real value is passed as an array of 2 floats; double scalar
 * operands are converted exactly where the ds_real format allows.
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
#ifndef _QD_C_DS_H
#define _QD_C_DS_H

#include <qd/qd_config.h>

#ifdef __cplusplus
extern "C" {
#endif

/* add */
void c_ds_add(const float *a, const float *b, float *c);
void c_ds_add_ds_ts(const float *a, const float *b, float *c);
void c_ds_add_ts_ds(const float *a, const float *b, float *c);
void c_ds_add_ds_qs(const float *a, const float *b, float *c);
void c_ds_add_qs_ds(const float *a, const float *b, float *c);
void c_ds_add_d_ds(double a, const float *b, float *c);
void c_ds_add_ds_d(const float *a, double b, float *c);

/* sub */
void c_ds_sub(const float *a, const float *b, float *c);
void c_ds_sub_ds_ts(const float *a, const float *b, float *c);
void c_ds_sub_ts_ds(const float *a, const float *b, float *c);
void c_ds_sub_ds_qs(const float *a, const float *b, float *c);
void c_ds_sub_qs_ds(const float *a, const float *b, float *c);
void c_ds_sub_d_ds(double a, const float *b, float *c);
void c_ds_sub_ds_d(const float *a, double b, float *c);

/* mul */
void c_ds_mul(const float *a, const float *b, float *c);
void c_ds_mul_ds_ts(const float *a, const float *b, float *c);
void c_ds_mul_ts_ds(const float *a, const float *b, float *c);
void c_ds_mul_ds_qs(const float *a, const float *b, float *c);
void c_ds_mul_qs_ds(const float *a, const float *b, float *c);
void c_ds_mul_d_ds(double a, const float *b, float *c);
void c_ds_mul_ds_d(const float *a, double b, float *c);

/* div */
void c_ds_div(const float *a, const float *b, float *c);
void c_ds_div_ds_ts(const float *a, const float *b, float *c);
void c_ds_div_ts_ds(const float *a, const float *b, float *c);
void c_ds_div_ds_qs(const float *a, const float *b, float *c);
void c_ds_div_qs_ds(const float *a, const float *b, float *c);
void c_ds_div_d_ds(double a, const float *b, float *c);
void c_ds_div_ds_d(const float *a, double b, float *c);

/* copy */
void c_ds_copy(const float *a, float *b);
void c_ds_copy_ts(const float *a, float *b);
void c_ds_copy_qs(const float *a, float *b);
void c_ds_copy_d(double a, float *b);

void c_ds_selfadd(const float *a, float *b);
void c_ds_selfadd_ts(const float *a, float *b);
void c_ds_selfadd_qs(const float *a, float *b);
void c_ds_selfadd_d(double a, float *b);
void c_ds_selfsub(const float *a, float *b);
void c_ds_selfsub_ts(const float *a, float *b);
void c_ds_selfsub_qs(const float *a, float *b);
void c_ds_selfsub_d(double a, float *b);
void c_ds_selfmul(const float *a, float *b);
void c_ds_selfmul_ts(const float *a, float *b);
void c_ds_selfmul_qs(const float *a, float *b);
void c_ds_selfmul_d(double a, float *b);
void c_ds_selfdiv(const float *a, float *b);
void c_ds_selfdiv_ts(const float *a, float *b);
void c_ds_selfdiv_qs(const float *a, float *b);
void c_ds_selfdiv_d(double a, float *b);

void c_ds_sqrt(const float *a, float *b);
void c_ds_sqr(const float *a, float *b);
void c_ds_abs(const float *a, float *b);
void c_ds_nint(const float *a, float *b);
void c_ds_aint(const float *a, float *b);
void c_ds_floor(const float *a, float *b);
void c_ds_ceil(const float *a, float *b);
void c_ds_exp(const float *a, float *b);
void c_ds_log(const float *a, float *b);
void c_ds_log10(const float *a, float *b);
void c_ds_sin(const float *a, float *b);
void c_ds_cos(const float *a, float *b);
void c_ds_tan(const float *a, float *b);
void c_ds_asin(const float *a, float *b);
void c_ds_acos(const float *a, float *b);
void c_ds_atan(const float *a, float *b);
void c_ds_sinh(const float *a, float *b);
void c_ds_cosh(const float *a, float *b);
void c_ds_tanh(const float *a, float *b);
void c_ds_asinh(const float *a, float *b);
void c_ds_acosh(const float *a, float *b);
void c_ds_atanh(const float *a, float *b);
void c_ds_npwr(const float *a, int n, float *b);
void c_ds_nroot(const float *a, int n, float *b);
void c_ds_atan2(const float *a, const float *b, float *c);
void c_ds_sincos(const float *a, float *s, float *c);
void c_ds_sincosh(const float *a, float *s, float *c);

void c_ds_read(const char *s, float *a);
void c_ds_swrite(const float *a, int precision, char *s, int len);
void c_ds_write(const float *a);
void c_ds_neg(const float *a, float *b);
void c_ds_rand(float *a);
void c_ds_comp(const float *a, const float *b, int *result);
void c_ds_comp_ds_d(const float *a, double b, int *result);
void c_ds_comp_d_ds(double a, const float *b, int *result);
void c_ds_comp_ds_ts(const float *a, const float *b, int *result);
void c_ds_comp_ts_ds(const float *a, const float *b, int *result);
void c_ds_comp_ds_qs(const float *a, const float *b, int *result);
void c_ds_comp_qs_ds(const float *a, const float *b, int *result);
void c_ds_pi(float *a);
void c_ds_2pi(float *a);
double c_ds_epsilon(void);

#ifdef __cplusplus
}
#endif

#endif /* _QD_C_DS_H */
