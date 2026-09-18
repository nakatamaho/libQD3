/*
 * fortran/f_ts.cpp
 *
 * C++ wrapper functions for triple-float precision arithmetic.
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

#include <qd/ts_real.h>
#include <qd/inline.h>

#define f_ts_add          FC_FUNC_(f_ts_add, F_TS_ADD)
#define f_ts_add_ts_d     FC_FUNC_(f_ts_add_ts_d, F_TS_ADD_TS_D)

#define f_ts_sub          FC_FUNC_(f_ts_sub, F_TS_SUB)
#define f_ts_sub_ts_d     FC_FUNC_(f_ts_sub_ts_d, F_TS_SUB_TS_D)
#define f_ts_sub_d_ts     FC_FUNC_(f_ts_sub_d_ts, F_TS_SUB_D_TS)

#define f_ts_mul          FC_FUNC_(f_ts_mul, F_TS_MUL)
#define f_ts_mul_ts_d     FC_FUNC_(f_ts_mul_ts_d, F_TS_MUL_TS_D)

#define f_ts_div          FC_FUNC_(f_ts_div, F_TS_DIV)
#define f_ts_div_ts_d     FC_FUNC_(f_ts_div_ts_d, F_TS_DIV_TS_D)
#define f_ts_div_d_ts     FC_FUNC_(f_ts_div_d_ts, F_TS_DIV_D_TS)

#define f_ts_sqrt         FC_FUNC_(f_ts_sqrt, F_TS_SQRT)
#define f_ts_sqr          FC_FUNC_(f_ts_sqr, F_TS_SQR)
#define f_ts_abs          FC_FUNC_(f_ts_abs, F_TS_ABS)

#define f_ts_npwr         FC_FUNC_(f_ts_npwr, F_TS_NPWR)
#define f_ts_nroot        FC_FUNC_(f_ts_nroot, F_TS_NROOT)
#define f_ts_nint         FC_FUNC_(f_ts_nint, F_TS_NINT)
#define f_ts_aint         FC_FUNC_(f_ts_aint, F_TS_AINT)
#define f_ts_floor        FC_FUNC_(f_ts_floor, F_TS_FLOOR)
#define f_ts_ceil         FC_FUNC_(f_ts_ceil, F_TS_CEIL)

#define f_ts_exp          FC_FUNC_(f_ts_exp, F_TS_EXP)
#define f_ts_log          FC_FUNC_(f_ts_log, F_TS_LOG)
#define f_ts_log10        FC_FUNC_(f_ts_log10, F_TS_LOG10)

#define f_ts_sin          FC_FUNC_(f_ts_sin, F_TS_SIN)
#define f_ts_cos          FC_FUNC_(f_ts_cos, F_TS_COS)
#define f_ts_tan          FC_FUNC_(f_ts_tan, F_TS_TAN)
#define f_ts_sincos       FC_FUNC_(f_ts_sincos, F_TS_SINCOS)

#define f_ts_asin         FC_FUNC_(f_ts_asin, F_TS_ASIN)
#define f_ts_acos         FC_FUNC_(f_ts_acos, F_TS_ACOS)
#define f_ts_atan         FC_FUNC_(f_ts_atan, F_TS_ATAN)
#define f_ts_atan2        FC_FUNC_(f_ts_atan2, F_TS_ATAN2)

#define f_ts_sinh         FC_FUNC_(f_ts_sinh, F_TS_SINH)
#define f_ts_cosh         FC_FUNC_(f_ts_cosh, F_TS_COSH)
#define f_ts_tanh         FC_FUNC_(f_ts_tanh, F_TS_TANH)
#define f_ts_sincosh      FC_FUNC_(f_ts_sincosh, F_TS_SINCOSH)

#define f_ts_asinh        FC_FUNC_(f_ts_asinh, F_TS_ASINH)
#define f_ts_acosh        FC_FUNC_(f_ts_acosh, F_TS_ACOSH)
#define f_ts_atanh        FC_FUNC_(f_ts_atanh, F_TS_ATANH)

#define f_ts_swrite       FC_FUNC_(f_ts_swrite, F_TS_SWRITE)
#define f_ts_write        FC_FUNC_(f_ts_write, F_TS_WRITE)
#define f_ts_neg          FC_FUNC_(f_ts_neg, F_TS_NEG)
#define f_ts_rand         FC_FUNC_(f_ts_rand, F_TS_RAND)
#define f_ts_comp         FC_FUNC_(f_ts_comp, F_TS_COMP)
#define f_ts_comp_ts_d    FC_FUNC_(f_ts_comp_ts_d, F_TS_COMP_TS_D)
#define f_ts_comp_d_ts    FC_FUNC_(f_ts_comp_d_ts, F_TS_COMP_D_TS)
#define f_ts_pi           FC_FUNC_(f_ts_pi, F_TS_PI)
#define f_ts_nan          FC_FUNC_(f_ts_nan, F_TS_NAN)

#define TO_FLOAT_PTR(a, ptr) \
  ptr[0] = (a)[0]; \
  ptr[1] = (a)[1]; \
  ptr[2] = (a)[2];

extern "C" {

static ts_real ts_floor_local(const ts_real &a) {
  return ::floor(a);
}

static ts_real ts_ceil_local(const ts_real &a) {
  return ::ceil(a);
}

static ts_real ts_aint_local(const ts_real &a) {
  return (a[0] >= 0.0) ? ts_floor_local(a) : ts_ceil_local(a);
}

static ts_real tsrand_local() {
  return ts_real::rand();
}

void f_ts_add(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ts_real(a) + ts_real(b), c);
}

void f_ts_add_ts_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ts_real(a) + *b, c);
}

void f_ts_sub(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ts_real(a) - ts_real(b), c);
}

void f_ts_sub_ts_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ts_real(a) - *b, c);
}

void f_ts_sub_d_ts(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(*a - ts_real(b), c);
}

void f_ts_mul(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ts_real(a) * ts_real(b), c);
}

void f_ts_mul_ts_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ts_real(a) * *b, c);
}

void f_ts_div(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ts_real(a) / ts_real(b), c);
}

void f_ts_div_ts_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(ts_real(a) / *b, c);
}

void f_ts_div_d_ts(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(*a / ts_real(b), c);
}

void f_ts_sqrt(const float *a, float *b) {
  TO_FLOAT_PTR(sqrt(ts_real(a)), b);
}

void f_ts_sqr(const float *a, float *b) {
  TO_FLOAT_PTR(sqr(ts_real(a)), b);
}

void f_ts_abs(const float *a, float *b) {
  TO_FLOAT_PTR(abs(ts_real(a)), b);
}

void f_ts_npwr(const float *a, const int *n, float *b) {
  TO_FLOAT_PTR(npwr(ts_real(a), *n), b);
}

void f_ts_nroot(const float *a, const int *n, float *b) {
  TO_FLOAT_PTR(nroot(ts_real(a), *n), b);
}

void f_ts_nint(const float *a, float *b) {
  TO_FLOAT_PTR(nint(ts_real(a)), b);
}

void f_ts_aint(const float *a, float *b) {
  TO_FLOAT_PTR(ts_aint_local(ts_real(a)), b);
}

void f_ts_floor(const float *a, float *b) {
  TO_FLOAT_PTR(ts_floor_local(ts_real(a)), b);
}

void f_ts_ceil(const float *a, float *b) {
  TO_FLOAT_PTR(ts_ceil_local(ts_real(a)), b);
}

void f_ts_exp(const float *a, float *b) {
  TO_FLOAT_PTR(exp(ts_real(a)), b);
}

void f_ts_log(const float *a, float *b) {
  TO_FLOAT_PTR(log(ts_real(a)), b);
}

void f_ts_log10(const float *a, float *b) {
  TO_FLOAT_PTR(log10(ts_real(a)), b);
}

void f_ts_sin(const float *a, float *b) {
  TO_FLOAT_PTR(sin(ts_real(a)), b);
}

void f_ts_cos(const float *a, float *b) {
  TO_FLOAT_PTR(cos(ts_real(a)), b);
}

void f_ts_tan(const float *a, float *b) {
  TO_FLOAT_PTR(tan(ts_real(a)), b);
}

void f_ts_sincos(const float *a, float *s, float *c) {
  ts_real ss, cc;
  sincos(ts_real(a), ss, cc);
  TO_FLOAT_PTR(ss, s);
  TO_FLOAT_PTR(cc, c);
}

void f_ts_asin(const float *a, float *b) {
  TO_FLOAT_PTR(asin(ts_real(a)), b);
}

void f_ts_acos(const float *a, float *b) {
  TO_FLOAT_PTR(acos(ts_real(a)), b);
}

void f_ts_atan(const float *a, float *b) {
  TO_FLOAT_PTR(atan(ts_real(a)), b);
}

void f_ts_atan2(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(atan2(ts_real(a), ts_real(b)), c);
}

void f_ts_sinh(const float *a, float *b) {
  TO_FLOAT_PTR(sinh(ts_real(a)), b);
}

void f_ts_cosh(const float *a, float *b) {
  TO_FLOAT_PTR(cosh(ts_real(a)), b);
}

void f_ts_tanh(const float *a, float *b) {
  TO_FLOAT_PTR(tanh(ts_real(a)), b);
}

void f_ts_sincosh(const float *a, float *s, float *c) {
  ts_real ss, cc;
  sincosh(ts_real(a), ss, cc);
  TO_FLOAT_PTR(ss, s);
  TO_FLOAT_PTR(cc, c);
}

void f_ts_asinh(const float *a, float *b) {
  TO_FLOAT_PTR(asinh(ts_real(a)), b);
}

void f_ts_acosh(const float *a, float *b) {
  TO_FLOAT_PTR(acosh(ts_real(a)), b);
}

void f_ts_atanh(const float *a, float *b) {
  TO_FLOAT_PTR(atanh(ts_real(a)), b);
}

void f_ts_swrite(const float *a, int *precision, char *s, int *maxlen) {
  int prec = *precision;
  if (prec <= 0 || prec > ts_real::_ndigits) prec = ts_real::_ndigits;
  std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0);
  std::string str = ts_real(a).to_string(prec, 0, fmt, false, true);

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

void f_ts_write(const float *a) {
  std::cout << ts_real(a) << std::endl;
}

void f_ts_neg(const float *a, float *b) {
  b[0] = -a[0];
  b[1] = -a[1];
  b[2] = -a[2];
}

void f_ts_rand(float *a) {
  TO_FLOAT_PTR(tsrand_local(), a);
}

void f_ts_comp(const float *a, const float *b, int *result) {
  ts_real aa(a), bb(b);
  if (aa < bb) {
    *result = -1;
  } else if (aa > bb) {
    *result = 1;
  } else {
    *result = 0;
  }
}

void f_ts_comp_ts_d(const float *a, const float *b, int *result) {
  ts_real aa(a);
  if (aa < *b) {
    *result = -1;
  } else if (aa > *b) {
    *result = 1;
  } else {
    *result = 0;
  }
}

void f_ts_comp_d_ts(const float *a, const float *b, int *result) {
  ts_real bb(b);
  if (*a < bb) {
    *result = -1;
  } else if (*a > bb) {
    *result = 1;
  } else {
    *result = 0;
  }
}

void f_ts_pi(float *a) {
  TO_FLOAT_PTR(ts_real::_pi, a);
}

void f_ts_nan(float *a) {
  TO_FLOAT_PTR(ts_real::_nan, a);
}

}

#endif /* HAVE_FORTRAN */
