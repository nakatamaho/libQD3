/*
 * src/c_qs.cpp
 *
 * C wrapper functions for quad-single precision arithmetic.
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
#include <iostream>

#include "config.h"
#include <qd/ds_real.h>
#include <qd/ts_real.h>
#include <qd/qs_real.h>
#include <qd/c_qs.h>

namespace {

typedef single_real<4> real_t;

// Mixed operations are evaluated in the wider format and rounded once.
template <int A, int B> struct wider { typedef single_real<(A > B ? A : B)> type; };

inline real_t load(const float *a) { return real_t(a); }

template <int M> inline single_real<M> load_n(const float *a) {
  return single_real<M>(a);
}

inline void store(const real_t &v, float *p) {
  for (int i = 0; i < 4; ++i) p[i] = v.x[i];
}

// double has 53 bits, so ts_real and wider hold it exactly.
typedef single_real<(4 > 3 ? 4 : 3)> scalar_t;

inline int compare(int lt, int gt) { return lt ? -1 : (gt ? 1 : 0); }

} // namespace

extern "C" {

void c_qs_add(const float *a, const float *b, float *c) {
  store(load(a) + load(b), c);
}

void c_qs_add_qs_ds(const float *a, const float *b, float *c) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load(a)) + W(load_n<2>(b))), c);
}

void c_qs_add_ds_qs(const float *a, const float *b, float *c) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load_n<2>(a)) + W(load(b))), c);
}

void c_qs_add_qs_ts(const float *a, const float *b, float *c) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load(a)) + W(load_n<3>(b))), c);
}

void c_qs_add_ts_qs(const float *a, const float *b, float *c) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load_n<3>(a)) + W(load(b))), c);
}

void c_qs_add_d_qs(double a, const float *b, float *c) {
  store(real_t(scalar_t(a) + scalar_t(load(b))), c);
}

void c_qs_add_qs_d(const float *a, double b, float *c) {
  store(real_t(scalar_t(load(a)) + scalar_t(b)), c);
}

void c_qs_sub(const float *a, const float *b, float *c) {
  store(load(a) - load(b), c);
}

void c_qs_sub_qs_ds(const float *a, const float *b, float *c) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load(a)) - W(load_n<2>(b))), c);
}

void c_qs_sub_ds_qs(const float *a, const float *b, float *c) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load_n<2>(a)) - W(load(b))), c);
}

void c_qs_sub_qs_ts(const float *a, const float *b, float *c) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load(a)) - W(load_n<3>(b))), c);
}

void c_qs_sub_ts_qs(const float *a, const float *b, float *c) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load_n<3>(a)) - W(load(b))), c);
}

void c_qs_sub_d_qs(double a, const float *b, float *c) {
  store(real_t(scalar_t(a) - scalar_t(load(b))), c);
}

void c_qs_sub_qs_d(const float *a, double b, float *c) {
  store(real_t(scalar_t(load(a)) - scalar_t(b)), c);
}

void c_qs_mul(const float *a, const float *b, float *c) {
  store(load(a) * load(b), c);
}

void c_qs_mul_qs_ds(const float *a, const float *b, float *c) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load(a)) * W(load_n<2>(b))), c);
}

void c_qs_mul_ds_qs(const float *a, const float *b, float *c) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load_n<2>(a)) * W(load(b))), c);
}

void c_qs_mul_qs_ts(const float *a, const float *b, float *c) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load(a)) * W(load_n<3>(b))), c);
}

void c_qs_mul_ts_qs(const float *a, const float *b, float *c) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load_n<3>(a)) * W(load(b))), c);
}

void c_qs_mul_d_qs(double a, const float *b, float *c) {
  store(real_t(scalar_t(a) * scalar_t(load(b))), c);
}

void c_qs_mul_qs_d(const float *a, double b, float *c) {
  store(real_t(scalar_t(load(a)) * scalar_t(b)), c);
}

void c_qs_div(const float *a, const float *b, float *c) {
  store(load(a) / load(b), c);
}

void c_qs_div_qs_ds(const float *a, const float *b, float *c) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load(a)) / W(load_n<2>(b))), c);
}

void c_qs_div_ds_qs(const float *a, const float *b, float *c) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load_n<2>(a)) / W(load(b))), c);
}

void c_qs_div_qs_ts(const float *a, const float *b, float *c) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load(a)) / W(load_n<3>(b))), c);
}

void c_qs_div_ts_qs(const float *a, const float *b, float *c) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load_n<3>(a)) / W(load(b))), c);
}

void c_qs_div_d_qs(double a, const float *b, float *c) {
  store(real_t(scalar_t(a) / scalar_t(load(b))), c);
}

void c_qs_div_qs_d(const float *a, double b, float *c) {
  store(real_t(scalar_t(load(a)) / scalar_t(b)), c);
}

void c_qs_copy(const float *a, float *b) {
  store(load(a), b);
}

void c_qs_copy_ds(const float *a, float *b) {
  store(real_t(load_n<2>(a)), b);
}

void c_qs_copy_ts(const float *a, float *b) {
  store(real_t(load_n<3>(a)), b);
}

void c_qs_copy_d(double a, float *b) {
  store(real_t(a), b);
}

void c_qs_selfadd(const float *a, float *b) {
  store(load(b) + load(a), b);
}

void c_qs_selfadd_ds(const float *a, float *b) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load(b)) + W(load_n<2>(a))), b);
}

void c_qs_selfadd_ts(const float *a, float *b) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load(b)) + W(load_n<3>(a))), b);
}

void c_qs_selfadd_d(double a, float *b) {
  store(real_t(scalar_t(load(b)) + scalar_t(a)), b);
}

void c_qs_selfsub(const float *a, float *b) {
  store(load(b) - load(a), b);
}

void c_qs_selfsub_ds(const float *a, float *b) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load(b)) - W(load_n<2>(a))), b);
}

void c_qs_selfsub_ts(const float *a, float *b) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load(b)) - W(load_n<3>(a))), b);
}

void c_qs_selfsub_d(double a, float *b) {
  store(real_t(scalar_t(load(b)) - scalar_t(a)), b);
}

void c_qs_selfmul(const float *a, float *b) {
  store(load(b) * load(a), b);
}

void c_qs_selfmul_ds(const float *a, float *b) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load(b)) * W(load_n<2>(a))), b);
}

void c_qs_selfmul_ts(const float *a, float *b) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load(b)) * W(load_n<3>(a))), b);
}

void c_qs_selfmul_d(double a, float *b) {
  store(real_t(scalar_t(load(b)) * scalar_t(a)), b);
}

void c_qs_selfdiv(const float *a, float *b) {
  store(load(b) / load(a), b);
}

void c_qs_selfdiv_ds(const float *a, float *b) {
  typedef wider<4, 2>::type W;
  store(real_t(W(load(b)) / W(load_n<2>(a))), b);
}

void c_qs_selfdiv_ts(const float *a, float *b) {
  typedef wider<4, 3>::type W;
  store(real_t(W(load(b)) / W(load_n<3>(a))), b);
}

void c_qs_selfdiv_d(double a, float *b) {
  store(real_t(scalar_t(load(b)) / scalar_t(a)), b);
}

void c_qs_sqrt(const float *a, float *b) {
  store(sqrt(load(a)), b);
}

void c_qs_sqr(const float *a, float *b) {
  store(sqr(load(a)), b);
}

void c_qs_abs(const float *a, float *b) {
  store(abs(load(a)), b);
}

void c_qs_nint(const float *a, float *b) {
  store(nint(load(a)), b);
}

void c_qs_aint(const float *a, float *b) {
  store(aint(load(a)), b);
}

void c_qs_floor(const float *a, float *b) {
  store(floor(load(a)), b);
}

void c_qs_ceil(const float *a, float *b) {
  store(ceil(load(a)), b);
}

void c_qs_exp(const float *a, float *b) {
  store(exp(load(a)), b);
}

void c_qs_log(const float *a, float *b) {
  store(log(load(a)), b);
}

void c_qs_log10(const float *a, float *b) {
  store(log10(load(a)), b);
}

void c_qs_sin(const float *a, float *b) {
  store(sin(load(a)), b);
}

void c_qs_cos(const float *a, float *b) {
  store(cos(load(a)), b);
}

void c_qs_tan(const float *a, float *b) {
  store(tan(load(a)), b);
}

void c_qs_asin(const float *a, float *b) {
  store(asin(load(a)), b);
}

void c_qs_acos(const float *a, float *b) {
  store(acos(load(a)), b);
}

void c_qs_atan(const float *a, float *b) {
  store(atan(load(a)), b);
}

void c_qs_sinh(const float *a, float *b) {
  store(sinh(load(a)), b);
}

void c_qs_cosh(const float *a, float *b) {
  store(cosh(load(a)), b);
}

void c_qs_tanh(const float *a, float *b) {
  store(tanh(load(a)), b);
}

void c_qs_asinh(const float *a, float *b) {
  store(asinh(load(a)), b);
}

void c_qs_acosh(const float *a, float *b) {
  store(acosh(load(a)), b);
}

void c_qs_atanh(const float *a, float *b) {
  store(atanh(load(a)), b);
}

void c_qs_npwr(const float *a, int n, float *b) {
  store(npwr(load(a), n), b);
}

void c_qs_nroot(const float *a, int n, float *b) {
  store(nroot(load(a), n), b);
}

void c_qs_atan2(const float *a, const float *b, float *c) {
  store(atan2(load(a), load(b)), c);
}

void c_qs_sincos(const float *a, float *s, float *c) {
  real_t ss, cc;
  sincos(load(a), ss, cc);
  store(ss, s);
  store(cc, c);
}

void c_qs_sincosh(const float *a, float *s, float *c) {
  real_t ss, cc;
  sincosh(load(a), ss, cc);
  store(ss, s);
  store(cc, c);
}

void c_qs_read(const char *s, float *a) {
  store(real_t(s), a);
}

void c_qs_swrite(const float *a, int precision, char *s, int len) {
  load(a).write(s, len, precision);
}

void c_qs_write(const float *a) {
  std::cout << load(a).to_string(real_t::_ndigits) << std::endl;
}

void c_qs_neg(const float *a, float *b) {
  store(-load(a), b);
}

void c_qs_rand(float *a) {
  store(real_t::rand(), a);
}

void c_qs_comp(const float *a, const float *b, int *result) {
  const real_t aa = load(a), bb = load(b);
  *result = compare(aa < bb, aa > bb);
}

void c_qs_comp_qs_d(const float *a, double b, int *result) {
  const scalar_t aa(load(a)), bb(b);
  *result = compare(aa < bb, aa > bb);
}

void c_qs_comp_d_qs(double a, const float *b, int *result) {
  const scalar_t aa(a), bb(load(b));
  *result = compare(aa < bb, aa > bb);
}

void c_qs_comp_qs_ds(const float *a, const float *b, int *result) {
  typedef wider<4, 2>::type W;
  const W aa(load(a)), bb(load_n<2>(b));
  *result = compare(aa < bb, aa > bb);
}

void c_qs_comp_ds_qs(const float *a, const float *b, int *result) {
  typedef wider<4, 2>::type W;
  const W aa(load_n<2>(a)), bb(load(b));
  *result = compare(aa < bb, aa > bb);
}

void c_qs_comp_qs_ts(const float *a, const float *b, int *result) {
  typedef wider<4, 3>::type W;
  const W aa(load(a)), bb(load_n<3>(b));
  *result = compare(aa < bb, aa > bb);
}

void c_qs_comp_ts_qs(const float *a, const float *b, int *result) {
  typedef wider<4, 3>::type W;
  const W aa(load_n<3>(a)), bb(load(b));
  *result = compare(aa < bb, aa > bb);
}

void c_qs_pi(float *a) {
  store(real_t::_pi, a);
}

void c_qs_2pi(float *a) {
  store(real_t::_2pi, a);
}

double c_qs_epsilon(void) {
  return static_cast<double>(real_t::_eps);
}

} // extern "C"
