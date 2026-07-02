#ifndef QD_ORACLE_FN_REGISTRY_H
#define QD_ORACLE_FN_REGISTRY_H

#include <cstddef>

namespace qd_oracle {

struct ErrorBound {
  const char *name;
  double eps_multiplier;
  const char *justification;
};

enum ArithOp {
  arith_add,
  arith_sub,
  arith_mul,
  arith_div,
  arith_sqr
};

struct ArithRegistryEntry {
  ArithOp op;
  const char *name;
  ErrorBound bound;
};

enum UnaryOp {
  unary_sqrt,
  unary_sqr,
  unary_exp,
  unary_exp2,
  unary_expm1,
  unary_log,
  unary_log10,
  unary_log1p,
  unary_log2,
  unary_sin,
  unary_cos,
  unary_tan,
  unary_asin,
  unary_acos,
  unary_atan,
  unary_sinh,
  unary_cosh,
  unary_tanh,
  unary_cbrt,
  unary_trunc,
  unary_round
};

enum InputDomain {
  domain_all_moderate,
  domain_nonnegative,
  domain_positive,
  domain_unit,
  domain_exp_moderate,
  domain_trig_moderate,
  domain_trig_stable,
  domain_trig_conditioned,
  domain_hyperbolic_moderate,
  domain_log1p
};

struct UnaryRegistryEntry {
  UnaryOp op;
  const char *name;
  InputDomain domain;
  ErrorBound bound;
};


enum BinaryOp {
  binary_pow_int,
  binary_pow_real,
  binary_nroot3,
  binary_atan2,
  binary_ldexp,
  binary_fmod,
  binary_drem,
  binary_hypot
};

inline ErrorBound binary_algebraic_bound(const char *name) {
  ErrorBound bound = {
      name,
      128.0,
      "binary algebraic operations use the same arithmetic kernels plus "
      "Newton/refinement steps; 128 library eps covers the M3 random domain"};
  return bound;
}

inline ErrorBound binary_transcendental_bound(const char *name) {
  ErrorBound bound = {
      name,
      512.0,
      "binary transcendental operations use argument reduction and Newton/Taylor "
      "refinement; 512 library eps matches the M2 transcendental triage bound"};
  return bound;
}

inline ErrorBound exact_ldexp_bound() {
  ErrorBound bound = {
      "ldexp",
      0.0,
      "ldexp scales every expansion limb by an exact power of two, so the "
      "result must match the exact MPFR reference"};
  return bound;
}

inline ErrorBound identity_bound(const char *name) {
  ErrorBound bound = {
      name,
      1024.0,
      "M3 identity checks compose multiple library functions; 1024 library eps "
      "is a triage bound until per-identity conditioning is split further"};
  return bound;
}

inline ErrorBound exact_add_smoke_bound() {
  ErrorBound bound = {
      "add_exact_smoke",
      0.0,
      "M1 smoke uses exactly representable inputs, so a+b must be bit exact"};
  return bound;
}

inline ErrorBound add_sub_bound(const char *name) {
#if defined(QD_IEEE_ADD)
  ErrorBound bound = {
      name,
      8.0,
      "docs/qd.tex gives the IEEE-style expansion add bound as about 2 eps_qd; "
      "8 library eps allows renormalization and DD/TD/EDD analog margin"};
#else
  ErrorBound bound = {
      name,
      16.0,
      "docs/qd.tex gives the default Cray-style add bound as separate "
      "operand perturbations <= eps_qd; randomized M2 inputs avoid severe "
      "cancellation and allow 16 library eps"};
#endif
  return bound;
}

inline ErrorBound mul_bound(const char *name) {
#if defined(QD_SLOPPY_MUL)
  ErrorBound bound = {
      name,
      32.0,
      "docs/qd.tex describes sloppy QD multiplication as dropping low "
      "O(eps^4) terms; 32 library eps is the documented loose M2 bound"};
#else
  ErrorBound bound = {
      name,
      8.0,
      "docs/qd.tex states accurate multiplication satisfies an IEEE-style "
      "relative error near eps_qd; 8 library eps allows renormalization margin"};
#endif
  return bound;
}

inline ErrorBound div_bound(const char *name) {
#if defined(QD_SLOPPY_DIV)
  ErrorBound bound = {
      name,
      64.0,
      "docs/qd.tex and docs/td.tex describe sloppy division as one fewer "
      "quotient digit; 64 library eps is the loose variant bound"};
#else
  ErrorBound bound = {
      name,
      16.0,
      "docs/qd.tex and docs/td.tex describe accurate division as retaining "
      "one guard quotient digit; 16 library eps is the tight variant bound"};
#endif
  return bound;
}

inline const ArithRegistryEntry *arithmetic_registry(std::size_t *count) {
  static const ArithRegistryEntry entries[] = {
      {arith_add, "add", add_sub_bound("add")},
      {arith_sub, "sub", add_sub_bound("sub")},
      {arith_mul, "mul", mul_bound("mul")},
      {arith_div, "div", div_bound("div")},
      {arith_sqr, "sqr", mul_bound("sqr")}};
  *count = sizeof(entries) / sizeof(entries[0]);
  return entries;
}

inline ErrorBound algebraic_unary_bound(const char *name) {
  ErrorBound bound = {
      name,
      32.0,
      "algebraic unary operations use Newton/refinement plus one "
      "renormalization; 32 library eps covers the M2 random domain"};
  return bound;
}

inline ErrorBound transcendental_bound(const char *name) {
  ErrorBound bound = {
      name,
      512.0,
      "transcendentals use argument reduction plus Taylor/Newton refinement; "
      "512 library eps covers DD/TD/QD/EDD until per-function tight bounds are "
      "split in the corner-case phases"};
  return bound;
}

inline ErrorBound stable_trig_bound(const char *name) {
  ErrorBound bound = {
      name,
      64.0,
      "stable trig inputs require |sin(x)| and |cos(x)| >= 0.25, so relative "
      "conditioning is bounded and the M2 argument-reduction path should stay "
      "well below the broad 512 eps transcendental bound"};
  return bound;
}

inline ErrorBound conditioned_trig_bound(const char *name) {
  ErrorBound bound = {
      name,
      16.0,
      "conditioned trig inputs allow zeros/poles; the driver scales this base "
      "bound by max(1, relative condition number) and reports sin/cos/tan MPFR "
      "components for triage"};
  return bound;
}

inline const UnaryRegistryEntry *unary_registry(std::size_t *count) {
  static const UnaryRegistryEntry entries[] = {
      {unary_sqrt, "sqrt", domain_nonnegative,
       algebraic_unary_bound("sqrt")},
      {unary_sqr, "sqr", domain_all_moderate,
       algebraic_unary_bound("sqr")},
      {unary_exp, "exp", domain_exp_moderate,
       transcendental_bound("exp")},
      {unary_exp2, "exp2", domain_exp_moderate,
       transcendental_bound("exp2")},
      {unary_expm1, "expm1", domain_exp_moderate,
       transcendental_bound("expm1")},
      {unary_log, "log", domain_positive,
       transcendental_bound("log")},
      {unary_log10, "log10", domain_positive,
       transcendental_bound("log10")},
      {unary_log1p, "log1p", domain_log1p,
       transcendental_bound("log1p")},
      {unary_log2, "log2", domain_positive,
       transcendental_bound("log2")},
      {unary_sin, "sin_stable", domain_trig_stable,
       stable_trig_bound("sin_stable")},
      {unary_sin, "sin_conditioned", domain_trig_conditioned,
       conditioned_trig_bound("sin_conditioned")},
      {unary_cos, "cos_stable", domain_trig_stable,
       stable_trig_bound("cos_stable")},
      {unary_cos, "cos_conditioned", domain_trig_conditioned,
       conditioned_trig_bound("cos_conditioned")},
      {unary_tan, "tan_stable", domain_trig_stable,
       stable_trig_bound("tan_stable")},
      {unary_tan, "tan_conditioned", domain_trig_conditioned,
       conditioned_trig_bound("tan_conditioned")},
      {unary_asin, "asin", domain_unit,
       transcendental_bound("asin")},
      {unary_acos, "acos", domain_unit,
       transcendental_bound("acos")},
      {unary_atan, "atan", domain_all_moderate,
       transcendental_bound("atan")},
      {unary_sinh, "sinh", domain_hyperbolic_moderate,
       transcendental_bound("sinh")},
      {unary_cosh, "cosh", domain_hyperbolic_moderate,
       transcendental_bound("cosh")},
      {unary_tanh, "tanh", domain_hyperbolic_moderate,
       transcendental_bound("tanh")},
      {unary_cbrt, "cbrt", domain_all_moderate,
       algebraic_unary_bound("cbrt")},
      {unary_trunc, "trunc", domain_all_moderate,
       exact_add_smoke_bound()},
      {unary_round, "round", domain_all_moderate,
       exact_add_smoke_bound()}};
  *count = sizeof(entries) / sizeof(entries[0]);
  return entries;
}

} // namespace qd_oracle

#endif
