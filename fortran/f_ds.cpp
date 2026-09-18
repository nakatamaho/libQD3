/*
 * fortran/f_ds.cpp
 *
 * C++ wrapper functions for double-single precision arithmetic.
 * This can be used from Fortran code.
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
#include "config.h"
#ifdef HAVE_FORTRAN

#include <cstring>
#include <iostream>
#include <cstdlib>

#include <qd/ds_real.h>
#include <qd/inline.h>

#define f_ds_add          FC_FUNC_(f_ds_add, F_DS_ADD)
#define f_ds_add_ds_d     FC_FUNC_(f_ds_add_ds_d, F_DS_ADD_DS_D)

#define f_ds_sub          FC_FUNC_(f_ds_sub, F_DS_SUB)
#define f_ds_sub_ds_d     FC_FUNC_(f_ds_sub_ds_d, F_DS_SUB_DS_D)
#define f_ds_sub_d_ds     FC_FUNC_(f_ds_sub_d_ds, F_DS_SUB_D_DS)

#define f_ds_mul          FC_FUNC_(f_ds_mul, F_DS_MUL)
#define f_ds_mul_ds_d     FC_FUNC_(f_ds_mul_ds_d, F_DS_MUL_DS_D)

#define f_ds_div          FC_FUNC_(f_ds_div, F_DS_DIV)
#define f_ds_div_ds_d     FC_FUNC_(f_ds_div_ds_d, F_DS_DIV_DS_D)
#define f_ds_div_d_ds     FC_FUNC_(f_ds_div_d_ds, F_DS_DIV_D_DS)

#define f_ds_sqrt         FC_FUNC_(f_ds_sqrt, F_DS_SQRT)
#define f_ds_sqr          FC_FUNC_(f_ds_sqr, F_DS_SQR)
#define f_ds_abs          FC_FUNC_(f_ds_abs, F_DS_ABS)

#define f_ds_npwr         FC_FUNC_(f_ds_npwr, F_DS_NPWR)
#define f_ds_nroot        FC_FUNC_(f_ds_nroot, F_DS_NROOT)
#define f_ds_nint         FC_FUNC_(f_ds_nint, F_DS_NINT)
#define f_ds_aint         FC_FUNC_(f_ds_aint, F_DS_AINT)
#define f_ds_floor        FC_FUNC_(f_ds_floor, F_DS_FLOOR)
#define f_ds_ceil         FC_FUNC_(f_ds_ceil, F_DS_CEIL)

#define f_ds_exp          FC_FUNC_(f_ds_exp, F_DS_EXP)
#define f_ds_log          FC_FUNC_(f_ds_log, F_DS_LOG)
#define f_ds_log10        FC_FUNC_(f_ds_log10, F_DS_LOG10)

#define f_ds_sin          FC_FUNC_(f_ds_sin, F_DS_SIN)
#define f_ds_cos          FC_FUNC_(f_ds_cos, F_DS_COS)
#define f_ds_tan          FC_FUNC_(f_ds_tan, F_DS_TAN)
#define f_ds_sincos       FC_FUNC_(f_ds_sincos, F_DS_SINCOS)

#define f_ds_asin         FC_FUNC_(f_ds_asin, F_DS_ASIN)
#define f_ds_acos         FC_FUNC_(f_ds_acos, F_DS_ACOS)
#define f_ds_atan         FC_FUNC_(f_ds_atan, F_DS_ATAN)
#define f_ds_atan2        FC_FUNC_(f_ds_atan2, F_DS_ATAN2)

#define f_ds_sinh         FC_FUNC_(f_ds_sinh, F_DS_SINH)
#define f_ds_cosh         FC_FUNC_(f_ds_cosh, F_DS_COSH)
#define f_ds_tanh         FC_FUNC_(f_ds_tanh, F_DS_TANH)
#define f_ds_sincosh      FC_FUNC_(f_ds_sincosh, F_DS_SINCOSH)

#define f_ds_asinh        FC_FUNC_(f_ds_asinh, F_DS_ASINH)
#define f_ds_acosh        FC_FUNC_(f_ds_acosh, F_DS_ACOSH)
#define f_ds_atanh        FC_FUNC_(f_ds_atanh, F_DS_ATANH)

#define f_ds_swrite       FC_FUNC_(f_ds_swrite, F_DS_SWRITE)
#define f_ds_write        FC_FUNC_(f_ds_write, F_DS_WRITE)
#define f_ds_neg          FC_FUNC_(f_ds_neg, F_DS_NEG)
#define f_ds_rand         FC_FUNC_(f_ds_rand, F_DS_RAND)
#define f_ds_comp         FC_FUNC_(f_ds_comp, F_DS_COMP)
#define f_ds_comp_ds_d    FC_FUNC_(f_ds_comp_ds_d, F_DS_COMP_DS_D)
#define f_ds_comp_d_ds    FC_FUNC_(f_ds_comp_d_ds, F_DS_COMP_D_DS)
#define f_ds_pi           FC_FUNC_(f_ds_pi, F_DS_PI)
#define f_ds_nan          FC_FUNC_(f_ds_nan, F_DS_NAN)

#define TO_FLOAT_PTR(a, ptr) \
  ptr[0] = (a)[0]; \
  ptr[1] = (a)[1];

extern "C" {

static ds_real ds_floor_local(const ds_real &a) {
  return ::floor(a);
}

static ds_real ds_ceil_local(const ds_real &a) {
  return ::ceil(a);
}

static ds_real ds_aint_local(const ds_real &a) {
  return (a[0] >= 0.0) ? ds_floor_local(a) : ds_ceil_local(a);
}

static ds_real dsrand_local() {
  return ds_real::rand();
}

void f_ds_add(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ds_real(a) + ds_real(b), c);
}

void f_ds_add_ds_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ds_real(a) + *b, c);
}

void f_ds_sub(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ds_real(a) - ds_real(b), c);
}

void f_ds_sub_ds_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ds_real(a) - *b, c);
}

void f_ds_sub_d_ds(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(*a - ds_real(b), c);
}

void f_ds_mul(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ds_real(a) * ds_real(b), c);
}

void f_ds_mul_ds_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ds_real(a) * *b, c);
}

void f_ds_div(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ds_real(a) / ds_real(b), c);
}

void f_ds_div_ds_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ds_real(a) / *b, c);
}

void f_ds_div_d_ds(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(*a / ds_real(b), c);
}

void f_ds_sqrt(const float *a, float *b) {
  TO_FLOAT_PTR(sqrt(ds_real(a)), b);
}

void f_ds_sqr(const float *a, float *b) {
  TO_FLOAT_PTR(sqr(ds_real(a)), b);
}

void f_ds_abs(const float *a, float *b) {
  TO_FLOAT_PTR(abs(ds_real(a)), b);
}

void f_ds_npwr(const float *a, const int *n, float *b) {
  TO_FLOAT_PTR(npwr(ds_real(a), *n), b);
}

void f_ds_nroot(const float *a, const int *n, float *b) {
  TO_FLOAT_PTR(nroot(ds_real(a), *n), b);
}

void f_ds_nint(const float *a, float *b) {
  TO_FLOAT_PTR(nint(ds_real(a)), b);
}

void f_ds_aint(const float *a, float *b) {
  TO_FLOAT_PTR(ds_aint_local(ds_real(a)), b);
}

void f_ds_floor(const float *a, float *b) {
  TO_FLOAT_PTR(ds_floor_local(ds_real(a)), b);
}

void f_ds_ceil(const float *a, float *b) {
  TO_FLOAT_PTR(ds_ceil_local(ds_real(a)), b);
}

void f_ds_exp(const float *a, float *b) {
  TO_FLOAT_PTR(exp(ds_real(a)), b);
}

void f_ds_log(const float *a, float *b) {
  TO_FLOAT_PTR(log(ds_real(a)), b);
}

void f_ds_log10(const float *a, float *b) {
  TO_FLOAT_PTR(log10(ds_real(a)), b);
}

void f_ds_sin(const float *a, float *b) {
  TO_FLOAT_PTR(sin(ds_real(a)), b);
}

void f_ds_cos(const float *a, float *b) {
  TO_FLOAT_PTR(cos(ds_real(a)), b);
}

void f_ds_tan(const float *a, float *b) {
  TO_FLOAT_PTR(tan(ds_real(a)), b);
}

void f_ds_sincos(const float *a, float *s, float *c) {
  ds_real ss, cc;
  sincos(ds_real(a), ss, cc);
  TO_FLOAT_PTR(ss, s);
  TO_FLOAT_PTR(cc, c);
}

void f_ds_asin(const float *a, float *b) {
  TO_FLOAT_PTR(asin(ds_real(a)), b);
}

void f_ds_acos(const float *a, float *b) {
  TO_FLOAT_PTR(acos(ds_real(a)), b);
}

void f_ds_atan(const float *a, float *b) {
  TO_FLOAT_PTR(atan(ds_real(a)), b);
}

void f_ds_atan2(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(atan2(ds_real(a), ds_real(b)), c);
}

void f_ds_sinh(const float *a, float *b) {
  TO_FLOAT_PTR(sinh(ds_real(a)), b);
}

void f_ds_cosh(const float *a, float *b) {
  TO_FLOAT_PTR(cosh(ds_real(a)), b);
}

void f_ds_tanh(const float *a, float *b) {
  TO_FLOAT_PTR(tanh(ds_real(a)), b);
}

void f_ds_sincosh(const float *a, float *s, float *c) {
  ds_real ss, cc;
  sincosh(ds_real(a), ss, cc);
  TO_FLOAT_PTR(ss, s);
  TO_FLOAT_PTR(cc, c);
}

void f_ds_asinh(const float *a, float *b) {
  TO_FLOAT_PTR(asinh(ds_real(a)), b);
}

void f_ds_acosh(const float *a, float *b) {
  TO_FLOAT_PTR(acosh(ds_real(a)), b);
}

void f_ds_atanh(const float *a, float *b) {
  TO_FLOAT_PTR(atanh(ds_real(a)), b);
}

void f_ds_swrite(const float *a, int *precision, char *s, int *maxlen) {
  int prec = *precision;
  if (prec <= 0 || prec > ds_real::_ndigits) prec = ds_real::_ndigits;
  std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0);
  std::string str = ds_real(a).to_string(prec, 0, fmt, false, true);

  int len = 0;
  if (a[0] < 0.0) {
    strncpy(&s[len], str.c_str(), *maxlen - len);
  } else {
    s[len++] = ' ';
    strncpy(&s[len], str.c_str(), *maxlen - len);
  }

  len += str.length();
  for (int i = len; i < *maxlen; i++) s[i] = ' ';
}

void f_ds_write(const float *a) {
  std::cout << ds_real(a) << std::endl;
}

void f_ds_neg(const float *a, float *b) {
  b[0] = -a[0];
  b[1] = -a[1];
}

void f_ds_rand(float *a) {
  TO_FLOAT_PTR(dsrand_local(), a);
}

void f_ds_comp(const float *a, const float *b, int *result) {
  ds_real aa(a), bb(b);
  if (aa < bb) {
    *result = -1;
  } else if (aa > bb) {
    *result = 1;
  } else {
    *result = 0;
  }
}

void f_ds_comp_ds_d(const float *a, const float *b, int *result) {
  ds_real aa(a);
  if (aa < *b) {
    *result = -1;
  } else if (aa > *b) {
    *result = 1;
  } else {
    *result = 0;
  }
}

void f_ds_comp_d_ds(const float *a, const float *b, int *result) {
  ds_real bb(b);
  if (*a < bb) {
    *result = -1;
  } else if (*a > bb) {
    *result = 1;
  } else {
    *result = 0;
  }
}

void f_ds_pi(float *a) {
  TO_FLOAT_PTR(ds_real::_pi, a);
}

void f_ds_nan(float *a) {
  TO_FLOAT_PTR(ds_real::_nan, a);
}

}

#endif /* HAVE_FORTRAN */
