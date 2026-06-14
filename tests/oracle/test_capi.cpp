#include "mpfr_oracle.h"
#include "qd_rng.h"
#include "tap.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <qd/c_dd.h>
#include <qd/c_qd.h>
#include <qd/c_td.h>
#ifdef QD_HAVE_EDD_REAL
#include <qd/c_edd.h>
#endif
#include <qd/fpu.h>

namespace {

const int kSamples = 16;

struct Options {
  bool test_dd;
  bool test_td;
  bool test_qd;
  bool test_edd;
  bool verbose;
  bool has_seed;
  std::uint64_t seed;

  Options()
      : test_dd(false), test_td(false), test_qd(false), test_edd(false),
        verbose(false), has_seed(false), seed(0) {}
};

std::string int_text(int value) {
  std::ostringstream os;
  os << value;
  return os.str();
}

void print_usage() {
  std::cout << "oracle_test_capi [-dd] [-td] [-qd] [-edd] [-all] [-v]"
            << " [--seed=N]\n";
}

template <class T> struct CApi;

template <> struct CApi<dd_real> {
  typedef double Limb;
  typedef void (*UnaryFn)(const Limb *, Limb *);
  typedef void (*BinaryFn)(const Limb *, const Limb *, Limb *);
  static void copy(const Limb *a, Limb *b) { c_dd_copy(a, b); }
  static void copy_d(double a, Limb *b) { c_dd_copy_d(a, b); }
  static void add(const Limb *a, const Limb *b, Limb *c) { c_dd_add(a, b, c); }
  static void sub(const Limb *a, const Limb *b, Limb *c) { c_dd_sub(a, b, c); }
  static void mul(const Limb *a, const Limb *b, Limb *c) { c_dd_mul(a, b, c); }
  static void div(const Limb *a, const Limb *b, Limb *c) { c_dd_div(a, b, c); }
  static void sqrt_op(const Limb *a, Limb *b) { c_dd_sqrt(a, b); }
  static void sqr_op(const Limb *a, Limb *b) { c_dd_sqr(a, b); }
  static void abs_op(const Limb *a, Limb *b) { c_dd_abs(a, b); }
  static void neg_op(const Limb *a, Limb *b) { c_dd_neg(a, b); }
  static void exp_op(const Limb *a, Limb *b) { c_dd_exp(a, b); }
  static void log_op(const Limb *a, Limb *b) { c_dd_log(a, b); }
  static void log10_op(const Limb *a, Limb *b) { c_dd_log10(a, b); }
  static void sin_op(const Limb *a, Limb *b) { c_dd_sin(a, b); }
  static void cos_op(const Limb *a, Limb *b) { c_dd_cos(a, b); }
  static void tan_op(const Limb *a, Limb *b) { c_dd_tan(a, b); }
  static void atan2_op(const Limb *a, const Limb *b, Limb *c) {
    c_dd_atan2(a, b, c);
  }
  static void sinh_op(const Limb *a, Limb *b) { c_dd_sinh(a, b); }
  static void cosh_op(const Limb *a, Limb *b) { c_dd_cosh(a, b); }
  static void tanh_op(const Limb *a, Limb *b) { c_dd_tanh(a, b); }
  static void npwr_op(const Limb *a, int n, Limb *b) { c_dd_npwr(a, n, b); }
  static void nroot_op(const Limb *a, int n, Limb *b) { c_dd_nroot(a, n, b); }
  static void nint_op(const Limb *a, Limb *b) { c_dd_nint(a, b); }
  static void aint_op(const Limb *a, Limb *b) { c_dd_aint(a, b); }
  static void floor_op(const Limb *a, Limb *b) { c_dd_floor(a, b); }
  static void ceil_op(const Limb *a, Limb *b) { c_dd_ceil(a, b); }
  static void read(const char *s, Limb *a) { c_dd_read(s, a); }
  static void swrite(const Limb *a, int precision, char *s, int len) {
    c_dd_swrite(a, precision, s, len);
  }
  static void comp(const Limb *a, const Limb *b, int *result) {
    c_dd_comp(a, b, result);
  }
  static void pi(Limb *a) { c_dd_pi(a); }
  static void two_pi(Limb *a) { c_dd_2pi(a); }
  static Limb epsilon() { return c_dd_epsilon(); }
};

template <> struct CApi<td_real> {
  typedef double Limb;
  typedef void (*UnaryFn)(const Limb *, Limb *);
  typedef void (*BinaryFn)(const Limb *, const Limb *, Limb *);
  static void copy(const Limb *a, Limb *b) { c_td_copy(a, b); }
  static void copy_d(double a, Limb *b) { c_td_copy_d(a, b); }
  static void add(const Limb *a, const Limb *b, Limb *c) { c_td_add(a, b, c); }
  static void sub(const Limb *a, const Limb *b, Limb *c) { c_td_sub(a, b, c); }
  static void mul(const Limb *a, const Limb *b, Limb *c) { c_td_mul(a, b, c); }
  static void div(const Limb *a, const Limb *b, Limb *c) { c_td_div(a, b, c); }
  static void sqrt_op(const Limb *a, Limb *b) { c_td_sqrt(a, b); }
  static void sqr_op(const Limb *a, Limb *b) { c_td_sqr(a, b); }
  static void abs_op(const Limb *a, Limb *b) { c_td_abs(a, b); }
  static void neg_op(const Limb *a, Limb *b) { c_td_neg(a, b); }
  static void exp_op(const Limb *a, Limb *b) { c_td_exp(a, b); }
  static void log_op(const Limb *a, Limb *b) { c_td_log(a, b); }
  static void log10_op(const Limb *a, Limb *b) { c_td_log10(a, b); }
  static void sin_op(const Limb *a, Limb *b) { c_td_sin(a, b); }
  static void cos_op(const Limb *a, Limb *b) { c_td_cos(a, b); }
  static void tan_op(const Limb *a, Limb *b) { c_td_tan(a, b); }
  static void atan2_op(const Limb *a, const Limb *b, Limb *c) {
    c_td_atan2(a, b, c);
  }
  static void sinh_op(const Limb *a, Limb *b) { c_td_sinh(a, b); }
  static void cosh_op(const Limb *a, Limb *b) { c_td_cosh(a, b); }
  static void tanh_op(const Limb *a, Limb *b) { c_td_tanh(a, b); }
  static void npwr_op(const Limb *a, int n, Limb *b) { c_td_npwr(a, n, b); }
  static void read(const char *s, Limb *a) { c_td_read(s, a); }
  static void swrite(const Limb *a, int precision, char *s, int len) {
    c_td_swrite(a, precision, s, len);
  }
  static void comp(const Limb *a, const Limb *b, int *result) {
    c_td_comp(a, b, result);
  }
  static void pi(Limb *a) { c_td_pi(a); }
  static void two_pi(Limb *a) { c_td_2pi(a); }
  static Limb epsilon() { return c_td_epsilon(); }
};

template <> struct CApi<qd_real> {
  typedef double Limb;
  typedef void (*UnaryFn)(const Limb *, Limb *);
  typedef void (*BinaryFn)(const Limb *, const Limb *, Limb *);
  static void copy(const Limb *a, Limb *b) { c_qd_copy(a, b); }
  static void copy_d(double a, Limb *b) { c_qd_copy_d(a, b); }
  static void add(const Limb *a, const Limb *b, Limb *c) { c_qd_add(a, b, c); }
  static void sub(const Limb *a, const Limb *b, Limb *c) { c_qd_sub(a, b, c); }
  static void mul(const Limb *a, const Limb *b, Limb *c) { c_qd_mul(a, b, c); }
  static void div(const Limb *a, const Limb *b, Limb *c) { c_qd_div(a, b, c); }
  static void sqrt_op(const Limb *a, Limb *b) { (void) c_qd_sqrt(a, b); }
  static void sqr_op(const Limb *a, Limb *b) { c_qd_sqr(a, b); }
  static void abs_op(const Limb *a, Limb *b) { c_qd_abs(a, b); }
  static void neg_op(const Limb *a, Limb *b) { c_qd_neg(a, b); }
  static void exp_op(const Limb *a, Limb *b) { c_qd_exp(a, b); }
  static void log_op(const Limb *a, Limb *b) { c_qd_log(a, b); }
  static void log10_op(const Limb *a, Limb *b) { c_qd_log10(a, b); }
  static void sin_op(const Limb *a, Limb *b) { c_qd_sin(a, b); }
  static void cos_op(const Limb *a, Limb *b) { c_qd_cos(a, b); }
  static void tan_op(const Limb *a, Limb *b) { c_qd_tan(a, b); }
  static void atan2_op(const Limb *a, const Limb *b, Limb *c) {
    c_qd_atan2(a, b, c);
  }
  static void sinh_op(const Limb *a, Limb *b) { c_qd_sinh(a, b); }
  static void cosh_op(const Limb *a, Limb *b) { c_qd_cosh(a, b); }
  static void tanh_op(const Limb *a, Limb *b) { c_qd_tanh(a, b); }
  static void npwr_op(const Limb *a, int n, Limb *b) { c_qd_npwr(a, n, b); }
  static void nroot_op(const Limb *a, int n, Limb *b) { c_qd_nroot(a, n, b); }
  static void nint_op(const Limb *a, Limb *b) { c_qd_nint(a, b); }
  static void aint_op(const Limb *a, Limb *b) { c_qd_aint(a, b); }
  static void floor_op(const Limb *a, Limb *b) { c_qd_floor(a, b); }
  static void ceil_op(const Limb *a, Limb *b) { c_qd_ceil(a, b); }
  static void read(const char *s, Limb *a) { c_qd_read(s, a); }
  static void swrite(const Limb *a, int precision, char *s, int len) {
    c_qd_swrite(a, precision, s, len);
  }
  static void comp(const Limb *a, const Limb *b, int *result) {
    c_qd_comp(a, b, result);
  }
  static void pi(Limb *a) { c_qd_pi(a); }
  static void two_pi(Limb *a) { c_qd_2pi(a); }
  static Limb epsilon() { return c_qd_epsilon(); }
};

#ifdef QD_HAVE_EDD_REAL
template <> struct CApi<edd_real> {
  typedef edd_word Limb;
  typedef void (*UnaryFn)(const Limb *, Limb *);
  typedef void (*BinaryFn)(const Limb *, const Limb *, Limb *);
  static void copy(const Limb *a, Limb *b) { c_edd_copy(a, b); }
  static void copy_d(double a, Limb *b) { c_edd_copy_d(a, b); }
  static void add(const Limb *a, const Limb *b, Limb *c) { c_edd_add(a, b, c); }
  static void sub(const Limb *a, const Limb *b, Limb *c) { c_edd_sub(a, b, c); }
  static void mul(const Limb *a, const Limb *b, Limb *c) { c_edd_mul(a, b, c); }
  static void div(const Limb *a, const Limb *b, Limb *c) { c_edd_div(a, b, c); }
  static void sqrt_op(const Limb *a, Limb *b) { c_edd_sqrt(a, b); }
  static void sqr_op(const Limb *a, Limb *b) { c_edd_sqr(a, b); }
  static void abs_op(const Limb *a, Limb *b) { c_edd_abs(a, b); }
  static void neg_op(const Limb *a, Limb *b) { c_edd_neg(a, b); }
  static void exp_op(const Limb *a, Limb *b) { c_edd_exp(a, b); }
  static void log_op(const Limb *a, Limb *b) { c_edd_log(a, b); }
  static void log10_op(const Limb *a, Limb *b) { c_edd_log10(a, b); }
  static void sin_op(const Limb *a, Limb *b) { c_edd_sin(a, b); }
  static void cos_op(const Limb *a, Limb *b) { c_edd_cos(a, b); }
  static void tan_op(const Limb *a, Limb *b) { c_edd_tan(a, b); }
  static void atan2_op(const Limb *a, const Limb *b, Limb *c) {
    c_edd_atan2(a, b, c);
  }
  static void sinh_op(const Limb *a, Limb *b) { c_edd_sinh(a, b); }
  static void cosh_op(const Limb *a, Limb *b) { c_edd_cosh(a, b); }
  static void tanh_op(const Limb *a, Limb *b) { c_edd_tanh(a, b); }
  static void read(const char *s, Limb *a) { c_edd_read(s, a); }
  static void swrite(const Limb *a, int precision, char *s, int len) {
    c_edd_swrite(a, precision, s, len);
  }
  static void comp(const Limb *a, const Limb *b, int *result) {
    c_edd_comp(a, b, result);
  }
  static void pi(Limb *a) { c_edd_pi(a); }
  static void two_pi(Limb *a) { c_edd_2pi(a); }
  static Limb epsilon() { return c_edd_epsilon(); }
};
#endif

template <class T>
void to_array(const T &value, typename CApi<T>::Limb *out) {
  for (int i = 0; i < qd_oracle::TypeTraits<T>::limbs; ++i) {
    out[i] = qd_oracle::TypeTraits<T>::limb(value, i);
  }
}

template <class T>
T from_array(const typename CApi<T>::Limb *in) {
  typename qd_oracle::TypeTraits<T>::limb_type limbs[qd_oracle::TypeTraits<T>::limbs];
  for (int i = 0; i < qd_oracle::TypeTraits<T>::limbs; ++i) {
    limbs[i] = in[i];
  }
  return qd_oracle::TypeTraits<T>::make(limbs);
}

template <class T>
std::vector<qd_oracle::Tap::Diagnostic>
base_diag(const char *case_name, const char *op, int iteration) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  std::ostringstream replay;
  replay << "tests/oracle/test_capi -" << traits::name()
         << " --seed=" << qd_oracle::rng::active_seed();
  diag.push_back(qd_oracle::Tap::Diagnostic(
      "seed", std::to_string(qd_oracle::rng::active_seed())));
  diag.push_back(qd_oracle::Tap::Diagnostic("replay", replay.str()));
  diag.push_back(qd_oracle::Tap::Diagnostic("case", case_name));
  diag.push_back(qd_oracle::Tap::Diagnostic("operation", op));
  diag.push_back(qd_oracle::Tap::Diagnostic("iteration", int_text(iteration)));
  return diag;
}

template <class T>
void append_value_diag(std::vector<qd_oracle::Tap::Diagnostic> *diag,
                       const T &input_a, const T &input_b, const T &expected,
                       const T &got) {
  diag->push_back(qd_oracle::Tap::Diagnostic("input_a_limbs",
                                             qd_oracle::limbs_hex(input_a)));
  diag->push_back(qd_oracle::Tap::Diagnostic("input_b_limbs",
                                             qd_oracle::limbs_hex(input_b)));
  diag->push_back(qd_oracle::Tap::Diagnostic("expected_limbs",
                                             qd_oracle::limbs_hex(expected)));
  diag->push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                             qd_oracle::limbs_hex(got)));
}

template <class T>
bool record_failure(std::vector<qd_oracle::Tap::Diagnostic> *diag,
                    const char *case_name, const char *op, int iteration,
                    const T &input_a, const T &input_b, const T &expected,
                    const T &got) {
  if (diag->empty()) {
    *diag = base_diag<T>(case_name, op, iteration);
    append_value_diag(diag, input_a, input_b, expected, got);
  }
  return false;
}

template <class T>
bool check_unary(const char *case_name, const char *op, int iteration,
                 const T &input, const T &expected,
                 typename CApi<T>::UnaryFn fn,
                 std::vector<qd_oracle::Tap::Diagnostic> *diag) {
  typename CApi<T>::Limb a[qd_oracle::TypeTraits<T>::limbs];
  typename CApi<T>::Limb out[qd_oracle::TypeTraits<T>::limbs];
  to_array(input, a);
  fn(a, out);
  T got = from_array<T>(out);
  if (!qd_oracle::bit_equal(got, expected)) {
    return record_failure(diag, case_name, op, iteration, input, T(0),
                          expected, got);
  }
  return true;
}

template <class T>
bool check_binary(const char *case_name, const char *op, int iteration,
                  const T &a_value, const T &b_value, const T &expected,
                  typename CApi<T>::BinaryFn fn,
                  std::vector<qd_oracle::Tap::Diagnostic> *diag) {
  typename CApi<T>::Limb a[qd_oracle::TypeTraits<T>::limbs];
  typename CApi<T>::Limb b[qd_oracle::TypeTraits<T>::limbs];
  typename CApi<T>::Limb out[qd_oracle::TypeTraits<T>::limbs];
  to_array(a_value, a);
  to_array(b_value, b);
  fn(a, b, out);
  T got = from_array<T>(out);
  if (!qd_oracle::bit_equal(got, expected)) {
    return record_failure(diag, case_name, op, iteration, a_value, b_value,
                          expected, got);
  }
  return true;
}

template <class T>
bool finish_group(qd_oracle::Tap &tap, const char *case_name, bool pass,
                  const std::vector<qd_oracle::Tap::Diagnostic> &diag,
                  bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;
  std::string name = std::string(traits::name()) + " " + case_name;
  tap.ok(pass, name, diag);
  if (verbose || pass) {
    std::cout << "# " << traits::name() << " " << case_name
              << " samples=" << kSamples << "\n";
  }
  return pass;
}

template <class T>
bool run_copy_io(qd_oracle::Tap &tap, bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  typename CApi<T>::Limb out[traits::limbs];
  typename CApi<T>::Limb in[traits::limbs];

  CApi<T>::copy_d(1.25, out);
  T got = from_array<T>(out);
  if (!qd_oracle::bit_equal(got, T(1.25))) {
    pass = record_failure(&diag, "copy/read/swrite shims", "copy_d", 0,
                          T(1.25), T(0), T(1.25), got);
  }

  T source = qd_oracle::rng::uniform_type<T>(-4, 4);
  to_array(source, in);
  CApi<T>::copy(in, out);
  got = from_array<T>(out);
  if (!qd_oracle::bit_equal(got, source)) {
    pass = record_failure(&diag, "copy/read/swrite shims", "copy", 0, source,
                          T(0), source, got);
  }

  CApi<T>::read("1.25e-3", out);
  T expected("1.25e-3");
  got = from_array<T>(out);
  if (!qd_oracle::bit_equal(got, expected)) {
    pass = record_failure(&diag, "copy/read/swrite shims", "read", 0,
                          expected, T(0), expected, got);
  }

  char c_text[256];
  char cpp_text[256];
  CApi<T>::swrite(in, T::_ndigits, c_text, static_cast<int>(sizeof(c_text)));
  source.write(cpp_text, static_cast<int>(sizeof(cpp_text)), T::_ndigits);
  if (std::strcmp(c_text, cpp_text) != 0 && diag.empty()) {
    diag = base_diag<T>("copy/read/swrite shims", "swrite", 0);
    diag.push_back(qd_oracle::Tap::Diagnostic("input_limbs",
                                              qd_oracle::limbs_hex(source)));
    diag.push_back(qd_oracle::Tap::Diagnostic("c_text", c_text));
    diag.push_back(qd_oracle::Tap::Diagnostic("cpp_text", cpp_text));
    pass = false;
  }

  return finish_group<T>(tap, "copy/read/swrite shims", pass, diag, verbose);
}

template <class T>
bool run_arithmetic(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  for (int i = 0; i < kSamples; ++i) {
    T a = qd_oracle::rng::uniform_type<T>(-4, 4);
    T b = qd_oracle::rng::positive_type<T>(-4, 4) + T(0.125);
    pass &= check_binary<T>("arithmetic shims", "add", i, a, b, a + b,
                            CApi<T>::add, &diag);
    pass &= check_binary<T>("arithmetic shims", "sub", i, a, b, a - b,
                            CApi<T>::sub, &diag);
    pass &= check_binary<T>("arithmetic shims", "mul", i, a, b, a * b,
                            CApi<T>::mul, &diag);
    pass &= check_binary<T>("arithmetic shims", "div", i, a, b, a / b,
                            CApi<T>::div, &diag);
  }
  return finish_group<T>(tap, "arithmetic shims", pass, diag, verbose);
}

template <class T>
bool run_unary_algebraic(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  for (int i = 0; i < kSamples; ++i) {
    T a = qd_oracle::rng::uniform_type<T>(-4, 4);
    T p = abs(a) + T(0.125);
    pass &= check_unary<T>("algebraic unary shims", "sqrt", i, p, sqrt(p),
                           CApi<T>::sqrt_op, &diag);
    pass &= check_unary<T>("algebraic unary shims", "sqr", i, a, sqr(a),
                           CApi<T>::sqr_op, &diag);
    pass &= check_unary<T>("algebraic unary shims", "abs", i, a, abs(a),
                           CApi<T>::abs_op, &diag);
    pass &= check_unary<T>("algebraic unary shims", "neg", i, a, -a,
                           CApi<T>::neg_op, &diag);
  }
  return finish_group<T>(tap, "algebraic unary shims", pass, diag, verbose);
}

template <class T>
bool run_transcendental(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  for (int i = 0; i < kSamples; ++i) {
    T x = qd_oracle::rng::uniform_type<T>(-4, -1);
    T p = qd_oracle::rng::positive_type<T>(-4, 4) + T(0.5);
    pass &= check_unary<T>("transcendental shims", "exp", i, x, exp(x),
                           CApi<T>::exp_op, &diag);
    pass &= check_unary<T>("transcendental shims", "log", i, p, log(p),
                           CApi<T>::log_op, &diag);
    pass &= check_unary<T>("transcendental shims", "log10", i, p, log10(p),
                           CApi<T>::log10_op, &diag);
  }
  return finish_group<T>(tap, "transcendental shims", pass, diag, verbose);
}

template <class T>
bool run_trig(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  for (int i = 0; i < kSamples; ++i) {
    T x = qd_oracle::rng::uniform_type<T>(-4, -1);
    T y = qd_oracle::rng::uniform_type<T>(-3, 3);
    T z = qd_oracle::rng::positive_type<T>(-3, 3) + T(0.25);
    pass &= check_unary<T>("trig shims", "sin", i, x, sin(x),
                           CApi<T>::sin_op, &diag);
    pass &= check_unary<T>("trig shims", "cos", i, x, cos(x),
                           CApi<T>::cos_op, &diag);
    pass &= check_unary<T>("trig shims", "tan", i, x, tan(x),
                           CApi<T>::tan_op, &diag);
    pass &= check_binary<T>("trig shims", "atan2", i, y, z, atan2(y, z),
                            CApi<T>::atan2_op, &diag);
  }
  return finish_group<T>(tap, "trig shims", pass, diag, verbose);
}

template <class T>
bool run_hyperbolic(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  for (int i = 0; i < kSamples; ++i) {
    T x = qd_oracle::rng::uniform_type<T>(-4, -1);
    pass &= check_unary<T>("hyperbolic shims", "sinh", i, x, sinh(x),
                           CApi<T>::sinh_op, &diag);
    pass &= check_unary<T>("hyperbolic shims", "cosh", i, x, cosh(x),
                           CApi<T>::cosh_op, &diag);
    pass &= check_unary<T>("hyperbolic shims", "tanh", i, x, tanh(x),
                           CApi<T>::tanh_op, &diag);
  }
  return finish_group<T>(tap, "hyperbolic shims", pass, diag, verbose);
}

template <class T>
bool run_comp_constants(qd_oracle::Tap &tap, bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  typename CApi<T>::Limb a[traits::limbs];
  typename CApi<T>::Limb b[traits::limbs];
  typename CApi<T>::Limb out[traits::limbs];

  CApi<T>::pi(out);
  T got = from_array<T>(out);
  if (!qd_oracle::bit_equal(got, T::_pi)) {
    pass = record_failure(&diag, "comparison/constants shims", "pi", 0,
                          T::_pi, T(0), T::_pi, got);
  }

  CApi<T>::two_pi(out);
  got = from_array<T>(out);
  if (!qd_oracle::bit_equal(got, T::_2pi)) {
    pass = record_failure(&diag, "comparison/constants shims", "2pi", 0,
                          T::_2pi, T(0), T::_2pi, got);
  }

  T left(1.25);
  T right(2.0);
  to_array(left, a);
  to_array(right, b);
  int result = 0;
  CApi<T>::comp(a, b, &result);
  if (!(result < 0)) {
    if (diag.empty()) {
      diag = base_diag<T>("comparison/constants shims", "comp", 0);
      diag.push_back(qd_oracle::Tap::Diagnostic("input_a_limbs",
                                                qd_oracle::limbs_hex(left)));
      diag.push_back(qd_oracle::Tap::Diagnostic("input_b_limbs",
                                                qd_oracle::limbs_hex(right)));
      diag.push_back(qd_oracle::Tap::Diagnostic("result", int_text(result)));
    }
    pass = false;
  }

  typename CApi<T>::Limb eps = CApi<T>::epsilon();
  if (eps != static_cast<typename CApi<T>::Limb>(T::_eps)) {
    if (diag.empty()) {
      diag = base_diag<T>("comparison/constants shims", "epsilon", 0);
      diag.push_back(qd_oracle::Tap::Diagnostic(
          "result", qd_oracle::limbs_hex(T(static_cast<double>(eps)))));
    }
    pass = false;
  }

  return finish_group<T>(tap, "comparison/constants shims", pass, diag,
                         verbose);
}

template <class T>
bool run_power_rounding(qd_oracle::Tap &, bool) {
  return true;
}

template <>
bool run_power_rounding<dd_real>(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  for (int i = 0; i < kSamples; ++i) {
    dd_real a = qd_oracle::rng::uniform_type<dd_real>(-3, 3);
    dd_real p = abs(a) + dd_real(1.0);
    pass &= check_unary<dd_real>("power/rounding shims", "nint", i, a,
                                 nint(a), CApi<dd_real>::nint_op, &diag);
    pass &= check_unary<dd_real>("power/rounding shims", "aint", i, a,
                                 aint(a), CApi<dd_real>::aint_op, &diag);
    pass &= check_unary<dd_real>("power/rounding shims", "floor", i, a,
                                 floor(a), CApi<dd_real>::floor_op, &diag);
    pass &= check_unary<dd_real>("power/rounding shims", "ceil", i, a,
                                 ceil(a), CApi<dd_real>::ceil_op, &diag);

    double in[2];
    double out[2];
    to_array(p, in);
    CApi<dd_real>::npwr_op(in, 3, out);
    dd_real got = from_array<dd_real>(out);
    if (!qd_oracle::bit_equal(got, pow(p, 3))) {
      pass = record_failure(&diag, "power/rounding shims", "npwr", i, p,
                            dd_real(3), pow(p, 3), got);
    }
    CApi<dd_real>::nroot_op(in, 3, out);
    got = from_array<dd_real>(out);
    if (!qd_oracle::bit_equal(got, nroot(p, 3))) {
      pass = record_failure(&diag, "power/rounding shims", "nroot", i, p,
                            dd_real(3), nroot(p, 3), got);
    }
  }
  return finish_group<dd_real>(tap, "power/rounding shims", pass, diag,
                               verbose);
}

template <>
bool run_power_rounding<td_real>(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  for (int i = 0; i < kSamples; ++i) {
    td_real p = abs(qd_oracle::rng::uniform_type<td_real>(-3, 3)) + td_real(1.0);
    double in[3];
    double out[3];
    to_array(p, in);
    CApi<td_real>::npwr_op(in, 3, out);
    td_real got = from_array<td_real>(out);
    if (!qd_oracle::bit_equal(got, pow(p, 3))) {
      pass = record_failure(&diag, "power shims", "npwr", i, p, td_real(3),
                            pow(p, 3), got);
    }
  }
  return finish_group<td_real>(tap, "power shims", pass, diag, verbose);
}

template <>
bool run_power_rounding<qd_real>(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  for (int i = 0; i < kSamples; ++i) {
    qd_real a = qd_oracle::rng::uniform_type<qd_real>(-3, 3);
    qd_real p = abs(a) + qd_real(1.0);
    pass &= check_unary<qd_real>("power/rounding shims", "nint", i, a,
                                 nint(a), CApi<qd_real>::nint_op, &diag);
    pass &= check_unary<qd_real>("power/rounding shims", "aint", i, a,
                                 aint(a), CApi<qd_real>::aint_op, &diag);
    pass &= check_unary<qd_real>("power/rounding shims", "floor", i, a,
                                 floor(a), CApi<qd_real>::floor_op, &diag);
    pass &= check_unary<qd_real>("power/rounding shims", "ceil", i, a,
                                 ceil(a), CApi<qd_real>::ceil_op, &diag);

    double in[4];
    double out[4];
    to_array(p, in);
    CApi<qd_real>::npwr_op(in, 3, out);
    qd_real got = from_array<qd_real>(out);
    if (!qd_oracle::bit_equal(got, pow(p, 3))) {
      pass = record_failure(&diag, "power/rounding shims", "npwr", i, p,
                            qd_real(3), pow(p, 3), got);
    }
    CApi<qd_real>::nroot_op(in, 3, out);
    got = from_array<qd_real>(out);
    if (!qd_oracle::bit_equal(got, nroot(p, 3))) {
      pass = record_failure(&diag, "power/rounding shims", "nroot", i, p,
                            qd_real(3), nroot(p, 3), got);
    }
  }
  return finish_group<qd_real>(tap, "power/rounding shims", pass, diag,
                               verbose);
}

template <class T>
bool run_type(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  pass &= run_copy_io<T>(tap, verbose);
  pass &= run_arithmetic<T>(tap, verbose);
  pass &= run_unary_algebraic<T>(tap, verbose);
  pass &= run_transcendental<T>(tap, verbose);
  pass &= run_trig<T>(tap, verbose);
  pass &= run_hyperbolic<T>(tap, verbose);
  pass &= run_comp_constants<T>(tap, verbose);
  pass &= run_power_rounding<T>(tap, verbose);
  return pass;
}

int type_plan_dd_td_qd() { return 8; }
int type_plan_edd() { return 7; }

int selected_plan(const Options &options) {
  int count = 0;
  if (options.test_dd) count += type_plan_dd_td_qd();
  if (options.test_td) count += type_plan_dd_td_qd();
  if (options.test_qd) count += type_plan_dd_td_qd();
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) count += type_plan_edd();
#endif
  return count;
}

int selected_count(const Options &options) {
  int count = 0;
  if (options.test_dd) ++count;
  if (options.test_td) ++count;
  if (options.test_qd) ++count;
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) ++count;
#endif
  return count;
}

void select_all(Options *options) {
  options->test_dd = true;
  options->test_td = true;
  options->test_qd = true;
#ifdef QD_HAVE_EDD_REAL
  options->test_edd = true;
#endif
}

bool parse_args(int argc, char **argv, Options *options) {
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-h") == 0 ||
        std::strcmp(argv[i], "-help") == 0) {
      print_usage();
      std::exit(0);
    } else if (std::strcmp(argv[i], "-dd") == 0) {
      options->test_dd = true;
    } else if (std::strcmp(argv[i], "-td") == 0) {
      options->test_td = true;
    } else if (std::strcmp(argv[i], "-qd") == 0) {
      options->test_qd = true;
    } else if (std::strcmp(argv[i], "-edd") == 0) {
#ifdef QD_HAVE_EDD_REAL
      options->test_edd = true;
#else
      std::cerr << "edd_real is not enabled in this build\n";
      return false;
#endif
    } else if (std::strcmp(argv[i], "-all") == 0) {
      select_all(options);
    } else if (std::strcmp(argv[i], "-v") == 0 ||
               std::strcmp(argv[i], "-verbose") == 0) {
      options->verbose = true;
    } else {
      std::uint64_t seed = 0;
      if (qd_oracle::rng::parse_seed_arg(argv[i], &seed)) {
        options->has_seed = true;
        options->seed = seed;
      } else {
        std::cerr << "Unknown flag `" << argv[i] << "'.\n";
        return false;
      }
    }
  }
  if (selected_count(*options) == 0) select_all(options);
  return true;
}

} // namespace

int main(int argc, char **argv) {
  Options options;
  if (!parse_args(argc, argv, &options)) {
    print_usage();
    return 2;
  }
  qd_oracle::rng::configure(options.has_seed, options.seed);
  qd_oracle::Tap tap(selected_plan(options));
  std::cout << "# seed: " << qd_oracle::rng::active_seed() << "\n";

  bool pass = true;
  unsigned int old_cw = 0;
  bool fpu_fixed = true;
  fpu_fix_start(&old_cw);
  if (options.test_dd) pass &= run_type<dd_real>(tap, options.verbose);
  if (options.test_td) pass &= run_type<td_real>(tap, options.verbose);
  if (options.test_qd) pass &= run_type<qd_real>(tap, options.verbose);
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) {
    fpu_fix_end(&old_cw);
    fpu_fixed = false;
    pass &= run_type<edd_real>(tap, options.verbose);
  }
#endif
  if (fpu_fixed) fpu_fix_end(&old_cw);
  return pass ? tap.exit_status() : 1;
}
