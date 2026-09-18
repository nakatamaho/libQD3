/*
 * fortran/f_qs.cpp
 *
 * C++ wrapper functions for quad-single precision arithmetic.
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

#include <qd/qs_real.h>
#include <qd/inline.h>

#define f_qs_add          FC_FUNC_(f_qs_add, F_QS_ADD)
#define f_qs_add_qs_d     FC_FUNC_(f_qs_add_qs_d, F_QS_ADD_QS_D)

#define f_qs_sub          FC_FUNC_(f_qs_sub, F_QS_SUB)
#define f_qs_sub_qs_d     FC_FUNC_(f_qs_sub_qs_d, F_QS_SUB_QS_D)
#define f_qs_sub_d_qs     FC_FUNC_(f_qs_sub_d_qs, F_QS_SUB_D_QS)

#define f_qs_mul          FC_FUNC_(f_qs_mul, F_QS_MUL)
#define f_qs_mul_qs_d     FC_FUNC_(f_qs_mul_qs_d, F_QS_MUL_QS_D)

#define f_qs_div          FC_FUNC_(f_qs_div, F_QS_DIV)
#define f_qs_div_qs_d     FC_FUNC_(f_qs_div_qs_d, F_QS_DIV_QS_D)
#define f_qs_div_d_qs     FC_FUNC_(f_qs_div_d_qs, F_QS_DIV_D_QS)

#define f_qs_sqrt         FC_FUNC_(f_qs_sqrt, F_QS_SQRT)
#define f_qs_sqr          FC_FUNC_(f_qs_sqr, F_QS_SQR)
#define f_qs_abs          FC_FUNC_(f_qs_abs, F_QS_ABS)

#define f_qs_npwr         FC_FUNC_(f_qs_npwr, F_QS_NPWR)
#define f_qs_nroot        FC_FUNC_(f_qs_nroot, F_QS_NROOT)
#define f_qs_nint         FC_FUNC_(f_qs_nint, F_QS_NINT)
#define f_qs_aint         FC_FUNC_(f_qs_aint, F_QS_AINT)
#define f_qs_floor        FC_FUNC_(f_qs_floor, F_QS_FLOOR)
#define f_qs_ceil         FC_FUNC_(f_qs_ceil, F_QS_CEIL)

#define f_qs_exp          FC_FUNC_(f_qs_exp, F_QS_EXP)
#define f_qs_log          FC_FUNC_(f_qs_log, F_QS_LOG)
#define f_qs_log10        FC_FUNC_(f_qs_log10, F_QS_LOG10)

#define f_qs_sin          FC_FUNC_(f_qs_sin, F_QS_SIN)
#define f_qs_cos          FC_FUNC_(f_qs_cos, F_QS_COS)
#define f_qs_tan          FC_FUNC_(f_qs_tan, F_QS_TAN)
#define f_qs_sincos       FC_FUNC_(f_qs_sincos, F_QS_SINCOS)

#define f_qs_asin         FC_FUNC_(f_qs_asin, F_QS_ASIN)
#define f_qs_acos         FC_FUNC_(f_qs_acos, F_QS_ACOS)
#define f_qs_atan         FC_FUNC_(f_qs_atan, F_QS_ATAN)
#define f_qs_atan2        FC_FUNC_(f_qs_atan2, F_QS_ATAN2)

#define f_qs_sinh         FC_FUNC_(f_qs_sinh, F_QS_SINH)
#define f_qs_cosh         FC_FUNC_(f_qs_cosh, F_QS_COSH)
#define f_qs_tanh         FC_FUNC_(f_qs_tanh, F_QS_TANH)
#define f_qs_sincosh      FC_FUNC_(f_qs_sincosh, F_QS_SINCOSH)

#define f_qs_asinh        FC_FUNC_(f_qs_asinh, F_QS_ASINH)
#define f_qs_acosh        FC_FUNC_(f_qs_acosh, F_QS_ACOSH)
#define f_qs_atanh        FC_FUNC_(f_qs_atanh, F_QS_ATANH)

#define f_qs_swrite       FC_FUNC_(f_qs_swrite, F_QS_SWRITE)
#define f_qs_write        FC_FUNC_(f_qs_write, F_QS_WRITE)
#define f_qs_neg          FC_FUNC_(f_qs_neg, F_QS_NEG)
#define f_qs_rand         FC_FUNC_(f_qs_rand, F_QS_RAND)
#define f_qs_comp         FC_FUNC_(f_qs_comp, F_QS_COMP)
#define f_qs_comp_qs_d    FC_FUNC_(f_qs_comp_qs_d, F_QS_COMP_QS_D)
#define f_qs_comp_d_qs    FC_FUNC_(f_qs_comp_d_qs, F_QS_COMP_D_QS)
#define f_qs_pi           FC_FUNC_(f_qs_pi, F_QS_PI)
#define f_qs_nan          FC_FUNC_(f_qs_nan, F_QS_NAN)

#define TO_FLOAT_PTR(a, ptr) \
  ptr[0] = (a)[0]; \
  ptr[1] = (a)[1]; \
  ptr[2] = (a)[2]; \
  ptr[3] = (a)[3];

extern "C" {

static qs_real qs_floor_local(const qs_real &a) {
  return ::floor(a);
}

static qs_real qs_ceil_local(const qs_real &a) {
  return ::ceil(a);
}

static qs_real qs_aint_local(const qs_real &a) {
  return (a[0] >= 0.0) ? qs_floor_local(a) : qs_ceil_local(a);
}

static qs_real qsrand_local() {
  return qs_real::rand();
}

void f_qs_add(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(qs_real(a) + qs_real(b), c);
}

void f_qs_add_qs_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(qs_real(a) + *b, c);
}

void f_qs_sub(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(qs_real(a) - qs_real(b), c);
}

void f_qs_sub_qs_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(qs_real(a) - *b, c);
}

void f_qs_sub_d_qs(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(*a - qs_real(b), c);
}

void f_qs_mul(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(qs_real(a) * qs_real(b), c);
}

void f_qs_mul_qs_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(qs_real(a) * *b, c);
}

void f_qs_div(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(qs_real(a) / qs_real(b), c);
}

void f_qs_div_qs_d(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(qs_real(a) / *b, c);
}

void f_qs_div_d_qs(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(*a / qs_real(b), c);
}

void f_qs_sqrt(const float *a, float *b) {
  TO_FLOAT_PTR(sqrt(qs_real(a)), b);
}

void f_qs_sqr(const float *a, float *b) {
  TO_FLOAT_PTR(sqr(qs_real(a)), b);
}

void f_qs_abs(const float *a, float *b) {
  TO_FLOAT_PTR(abs(qs_real(a)), b);
}

void f_qs_npwr(const float *a, const int *n, float *b) {
  TO_FLOAT_PTR(npwr(qs_real(a), *n), b);
}

void f_qs_nroot(const float *a, const int *n, float *b) {
  TO_FLOAT_PTR(nroot(qs_real(a), *n), b);
}

void f_qs_nint(const float *a, float *b) {
  TO_FLOAT_PTR(nint(qs_real(a)), b);
}

void f_qs_aint(const float *a, float *b) {
  TO_FLOAT_PTR(qs_aint_local(qs_real(a)), b);
}

void f_qs_floor(const float *a, float *b) {
  TO_FLOAT_PTR(qs_floor_local(qs_real(a)), b);
}

void f_qs_ceil(const float *a, float *b) {
  TO_FLOAT_PTR(qs_ceil_local(qs_real(a)), b);
}

void f_qs_exp(const float *a, float *b) {
  TO_FLOAT_PTR(exp(qs_real(a)), b);
}

void f_qs_log(const float *a, float *b) {
  TO_FLOAT_PTR(log(qs_real(a)), b);
}

void f_qs_log10(const float *a, float *b) {
  TO_FLOAT_PTR(log10(qs_real(a)), b);
}

void f_qs_sin(const float *a, float *b) {
  TO_FLOAT_PTR(sin(qs_real(a)), b);
}

void f_qs_cos(const float *a, float *b) {
  TO_FLOAT_PTR(cos(qs_real(a)), b);
}

void f_qs_tan(const float *a, float *b) {
  TO_FLOAT_PTR(tan(qs_real(a)), b);
}

void f_qs_sincos(const float *a, float *s, float *c) {
  qs_real ss, cc;
  sincos(qs_real(a), ss, cc);
  TO_FLOAT_PTR(ss, s);
  TO_FLOAT_PTR(cc, c);
}

void f_qs_asin(const float *a, float *b) {
  TO_FLOAT_PTR(asin(qs_real(a)), b);
}

void f_qs_acos(const float *a, float *b) {
  TO_FLOAT_PTR(acos(qs_real(a)), b);
}

void f_qs_atan(const float *a, float *b) {
  TO_FLOAT_PTR(atan(qs_real(a)), b);
}

void f_qs_atan2(const float *a, const float *b, float *c) {
  TO_FLOAT_PTR(atan2(qs_real(a), qs_real(b)), c);
}

void f_qs_sinh(const float *a, float *b) {
  TO_FLOAT_PTR(sinh(qs_real(a)), b);
}

void f_qs_cosh(const float *a, float *b) {
  TO_FLOAT_PTR(cosh(qs_real(a)), b);
}

void f_qs_tanh(const float *a, float *b) {
  TO_FLOAT_PTR(tanh(qs_real(a)), b);
}

void f_qs_sincosh(const float *a, float *s, float *c) {
  qs_real ss, cc;
  sincosh(qs_real(a), ss, cc);
  TO_FLOAT_PTR(ss, s);
  TO_FLOAT_PTR(cc, c);
}

void f_qs_asinh(const float *a, float *b) {
  TO_FLOAT_PTR(asinh(qs_real(a)), b);
}

void f_qs_acosh(const float *a, float *b) {
  TO_FLOAT_PTR(acosh(qs_real(a)), b);
}

void f_qs_atanh(const float *a, float *b) {
  TO_FLOAT_PTR(atanh(qs_real(a)), b);
}

void f_qs_swrite(const float *a, int *precision, char *s, int *maxlen) {
  int prec = *precision;
  if (prec <= 0 || prec > qs_real::_ndigits) prec = qs_real::_ndigits;
  std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0);
  std::string str = qs_real(a).to_string(prec, 0, fmt, false, true);

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

void f_qs_write(const float *a) {
  std::cout << qs_real(a) << std::endl;
}

void f_qs_neg(const float *a, float *b) {
  b[0] = -a[0];
  b[1] = -a[1];
  b[2] = -a[2];
  b[3] = -a[3];
}

void f_qs_rand(float *a) {
  TO_FLOAT_PTR(qsrand_local(), a);
}

void f_qs_comp(const float *a, const float *b, int *result) {
  qs_real aa(a), bb(b);
  if (aa < bb) {
    *result = -1;
  } else if (aa > bb) {
    *result = 1;
  } else {
    *result = 0;
  }
}

void f_qs_comp_qs_d(const float *a, const float *b, int *result) {
  qs_real aa(a);
  if (aa < *b) {
    *result = -1;
  } else if (aa > *b) {
    *result = 1;
  } else {
    *result = 0;
  }
}

void f_qs_comp_d_qs(const float *a, const float *b, int *result) {
  qs_real bb(b);
  if (*a < bb) {
    *result = -1;
  } else if (*a > bb) {
    *result = 1;
  } else {
    *result = 0;
  }
}

void f_qs_pi(float *a) {
  TO_FLOAT_PTR(qs_real::_pi, a);
}

void f_qs_nan(float *a) {
  TO_FLOAT_PTR(qs_real::_nan, a);
}

}

#endif /* HAVE_FORTRAN */
