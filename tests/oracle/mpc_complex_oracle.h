#ifndef QD_ORACLE_MPC_COMPLEX_ORACLE_H
#define QD_ORACLE_MPC_COMPLEX_ORACLE_H

#include "fn_registry.h"
#include "mpfr_oracle.h"
#include "qd_rng.h"
#include "tap.h"

#include <mpc.h>

#include <algorithm>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include <qd/complex.h>

namespace qd_oracle {
namespace complex_oracle {

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

void print_usage(const char *program);
bool parse_options(int argc, char **argv, const char *program, Options *options);
int selected_count(const Options &options);
std::string double_text(double value);
std::string int_text(int value);
const char *rounding_name();

inline ErrorBound complex_arithmetic_bound(const char *name) {
  ErrorBound bound = {
      name,
      256.0,
      "complex arithmetic composes real expansion operations; 256 eps covers "
      "the bounded random MPC smoke domain"};
  return bound;
}

inline ErrorBound complex_transcendental_bound(const char *name) {
  ErrorBound bound = {
      name,
      4096.0,
      "complex elementary functions compose real transcendental kernels and "
      "MPC references over a stable bounded domain"};
  return bound;
}

inline ErrorBound complex_special_bound(const char *name) {
  ErrorBound bound = {
      name,
      1024.0,
      "complex special helpers compare against MPFR/MPC references on a "
      "deterministic bounded grid"};
  return bound;
}

template <class Real>
mpfr_prec_t oracle_prec() {
  const mpfr_prec_t min_prec = 1024;
  return std::max(min_prec, ref_prec<Real>() + 256);
}

template <class Real>
void init_mpc(mpc_t z) {
  mpc_init2(z, oracle_prec<Real>());
}

template <class Real>
void to_mpc(mpc_t out, const qd3_complex<Real> &z) {
  mpfr_t re;
  mpfr_t im;
  mpfr_inits2(oracle_prec<Real>(), re, im, (mpfr_ptr) 0);
  to_mpfr(re, z.real());
  to_mpfr(im, z.imag());
  mpc_set_fr_fr(out, re, im, MPC_RNDNN);
  mpfr_clears(re, im, (mpfr_ptr) 0);
}

template <class Real>
std::string complex_limbs_hex(const qd3_complex<Real> &z) {
  return std::string("re=") + limbs_hex(z.real()) + ", im=" +
         limbs_hex(z.imag());
}

template <class Real>
std::string complex_value_string(const qd3_complex<Real> &z) {
  return std::string("re=") + value_to_mpfr_string(z.real()) + ", im=" +
         value_to_mpfr_string(z.imag());
}

inline std::string mpfr_ptr_to_string(mpfr_srcptr value) {
  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(value));
  mpfr_set(tmp, value, MPFR_RNDN);
  std::string result = mpfr_to_string(tmp);
  mpfr_clear(tmp);
  return result;
}

template <class Real>
double complex_relerr_in_eps(const qd3_complex<Real> &got, mpc_t ref) {
  const double re = relerr_in_eps(got.real(), mpc_realref(ref));
  const double im = relerr_in_eps(got.imag(), mpc_imagref(ref));
  return std::max(re, im);
}

template <class Real>
std::string complex_abs_error_string(const qd3_complex<Real> &got, mpc_t ref) {
  return std::string("re=") + abs_error_to_string(got.real(), mpc_realref(ref)) +
         ", im=" + abs_error_to_string(got.imag(), mpc_imagref(ref));
}

template <class Real>
std::vector<Tap::Diagnostic> base_diag(const char *program,
                                       const char *case_name,
                                       int iteration) {
  typedef TypeTraits<Real> traits;
  std::vector<Tap::Diagnostic> diag;
  std::ostringstream replay;
  replay << "tests/oracle/" << program << " -" << traits::name()
         << " --seed=" << qd_oracle::rng::active_seed();
  diag.push_back(Tap::Diagnostic("seed", std::to_string(qd_oracle::rng::active_seed())));
  diag.push_back(Tap::Diagnostic("replay", replay.str()));
  diag.push_back(Tap::Diagnostic("case", case_name));
  diag.push_back(Tap::Diagnostic("iteration", int_text(iteration)));
  diag.push_back(Tap::Diagnostic("mpc_rounding", rounding_name()));
  return diag;
}

template <class Real>
void append_complex_diag(std::vector<Tap::Diagnostic> *diag,
                         const qd3_complex<Real> &input,
                         const qd3_complex<Real> &got,
                         mpc_t ref, double relerr,
                         const ErrorBound &bound) {
  diag->push_back(Tap::Diagnostic("input_limbs", complex_limbs_hex(input)));
  diag->push_back(Tap::Diagnostic("input_value", complex_value_string(input)));
  diag->push_back(Tap::Diagnostic("mpc_reference",
                                  std::string("re=") + mpfr_ptr_to_string(mpc_realref(ref)) +
                                      ", im=" + mpfr_ptr_to_string(mpc_imagref(ref))));
  diag->push_back(Tap::Diagnostic("got_value", complex_value_string(got)));
  diag->push_back(Tap::Diagnostic("got_limbs", complex_limbs_hex(got)));
  diag->push_back(Tap::Diagnostic("abs_error_mpfr", complex_abs_error_string(got, ref)));
  diag->push_back(Tap::Diagnostic("relerr_eps", double_text(relerr)));
  diag->push_back(Tap::Diagnostic("allowed_eps_multiplier", double_text(bound.eps_multiplier)));
  diag->push_back(Tap::Diagnostic("bound_justification", bound.justification));
}

template <class Real>
void append_complex_binary_diag(std::vector<Tap::Diagnostic> *diag,
                                const qd3_complex<Real> &a,
                                const qd3_complex<Real> &b,
                                const qd3_complex<Real> &got,
                                mpc_t ref, double relerr,
                                const ErrorBound &bound) {
  diag->push_back(Tap::Diagnostic("input_a_limbs", complex_limbs_hex(a)));
  diag->push_back(Tap::Diagnostic("input_a_value", complex_value_string(a)));
  diag->push_back(Tap::Diagnostic("input_b_limbs", complex_limbs_hex(b)));
  diag->push_back(Tap::Diagnostic("input_b_value", complex_value_string(b)));
  diag->push_back(Tap::Diagnostic("mpc_reference",
                                  std::string("re=") + mpfr_ptr_to_string(mpc_realref(ref)) +
                                      ", im=" + mpfr_ptr_to_string(mpc_imagref(ref))));
  diag->push_back(Tap::Diagnostic("got_value", complex_value_string(got)));
  diag->push_back(Tap::Diagnostic("got_limbs", complex_limbs_hex(got)));
  diag->push_back(Tap::Diagnostic("abs_error_mpfr", complex_abs_error_string(got, ref)));
  diag->push_back(Tap::Diagnostic("relerr_eps", double_text(relerr)));
  diag->push_back(Tap::Diagnostic("allowed_eps_multiplier", double_text(bound.eps_multiplier)));
  diag->push_back(Tap::Diagnostic("bound_justification", bound.justification));
}

template <class Real>
qd3_complex<Real> random_complex(int emin, int emax) {
  return qd3_complex<Real>(rng::uniform_type<Real>(emin, emax),
                            rng::uniform_type<Real>(emin, emax));
}

} // namespace complex_oracle
} // namespace qd_oracle

#endif
