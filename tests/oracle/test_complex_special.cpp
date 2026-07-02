#include "mpc_complex_oracle.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

const int kCasesPerType = 7;

template <class Real>
std::vector<qd_oracle::Tap::Diagnostic>
real_diag(const char *case_name, const qd3_complex<Real> &input,
          const Real &got, mpfr_t ref, double relerr,
          const qd_oracle::ErrorBound &bound) {
  std::vector<qd_oracle::Tap::Diagnostic> diag =
      qd_oracle::complex_oracle::base_diag<Real>("test_complex_special",
                                                 case_name, 0);
  diag.push_back(qd_oracle::Tap::Diagnostic(
      "input_limbs", qd_oracle::complex_oracle::complex_limbs_hex(input)));
  diag.push_back(qd_oracle::Tap::Diagnostic(
      "input_value", qd_oracle::complex_oracle::complex_value_string(input)));
  diag.push_back(qd_oracle::Tap::Diagnostic("mpfr_reference",
                                             qd_oracle::mpfr_to_string(ref)));
  diag.push_back(qd_oracle::Tap::Diagnostic("got_value",
                                             qd_oracle::value_to_mpfr_string(got)));
  diag.push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                             qd_oracle::limbs_hex(got)));
  diag.push_back(qd_oracle::Tap::Diagnostic("abs_error_mpfr",
                                             qd_oracle::abs_error_to_string(got, ref)));
  diag.push_back(qd_oracle::Tap::Diagnostic("relerr_eps",
                                             qd_oracle::complex_oracle::double_text(relerr)));
  diag.push_back(qd_oracle::Tap::Diagnostic(
      "allowed_eps_multiplier",
      qd_oracle::complex_oracle::double_text(bound.eps_multiplier)));
  diag.push_back(qd_oracle::Tap::Diagnostic("bound_justification",
                                             bound.justification));
  return diag;
}

template <class Real>
bool run_abs(qd_oracle::Tap &tap, bool verbose) {
  typedef qd3_complex<Real> Complex;
  using std::abs;
  qd_oracle::ErrorBound bound =
      qd_oracle::complex_oracle::complex_special_bound("complex abs");
  Complex z(Real(1.25), Real(-0.5));
  Real got = abs(z);
  mpc_t z_mpc;
  mpfr_t ref;
  qd_oracle::complex_oracle::init_mpc<Real>(z_mpc);
  mpfr_init2(ref, qd_oracle::complex_oracle::oracle_prec<Real>());
  qd_oracle::complex_oracle::to_mpc(z_mpc, z);
  mpc_abs(ref, z_mpc, MPFR_RNDN);
  double relerr = qd_oracle::relerr_in_eps(got, ref);
  bool pass = std::isfinite(relerr) && relerr <= bound.eps_multiplier;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) diag = real_diag("complex abs oracle", z, got, ref, relerr, bound);
  tap.ok(pass, std::string(qd_oracle::TypeTraits<Real>::name()) + " complex abs oracle",
         diag, verbose);
  mpfr_clear(ref);
  mpc_clear(z_mpc);
  return pass;
}

template <class Real>
bool run_norm(qd_oracle::Tap &tap, bool verbose) {
  typedef qd3_complex<Real> Complex;
  qd_oracle::ErrorBound bound =
      qd_oracle::complex_oracle::complex_special_bound("complex norm");
  Complex z(Real(1.25), Real(-0.5));
  Real got = norm(z);
  mpfr_t re;
  mpfr_t im;
  mpfr_t ref;
  mpfr_inits2(qd_oracle::complex_oracle::oracle_prec<Real>(), re, im, ref,
              (mpfr_ptr) 0);
  qd_oracle::to_mpfr(re, z.real());
  qd_oracle::to_mpfr(im, z.imag());
  mpfr_sqr(re, re, MPFR_RNDN);
  mpfr_sqr(im, im, MPFR_RNDN);
  mpfr_add(ref, re, im, MPFR_RNDN);
  double relerr = qd_oracle::relerr_in_eps(got, ref);
  bool pass = std::isfinite(relerr) && relerr <= bound.eps_multiplier;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) diag = real_diag("complex norm oracle", z, got, ref, relerr, bound);
  tap.ok(pass, std::string(qd_oracle::TypeTraits<Real>::name()) + " complex norm oracle",
         diag, verbose);
  mpfr_clears(re, im, ref, (mpfr_ptr) 0);
  return pass;
}

template <class Real>
bool run_arg(qd_oracle::Tap &tap, bool verbose) {
  typedef qd3_complex<Real> Complex;
  qd_oracle::ErrorBound bound =
      qd_oracle::complex_oracle::complex_special_bound("complex arg");
  Complex z(Real(-1.25), Real(0.5));
  Real got = arg(z);
  mpfr_t re;
  mpfr_t im;
  mpfr_t ref;
  mpfr_inits2(qd_oracle::complex_oracle::oracle_prec<Real>(), re, im, ref,
              (mpfr_ptr) 0);
  qd_oracle::to_mpfr(re, z.real());
  qd_oracle::to_mpfr(im, z.imag());
  mpfr_atan2(ref, im, re, MPFR_RNDN);
  double relerr = qd_oracle::relerr_in_eps(got, ref);
  bool pass = std::isfinite(relerr) && relerr <= bound.eps_multiplier;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) diag = real_diag("complex arg oracle", z, got, ref, relerr, bound);
  tap.ok(pass, std::string(qd_oracle::TypeTraits<Real>::name()) + " complex arg oracle",
         diag, verbose);
  mpfr_clears(re, im, ref, (mpfr_ptr) 0);
  return pass;
}

template <class Real>
bool run_polar(qd_oracle::Tap &tap, bool verbose) {
  typedef qd3_complex<Real> Complex;
  using std::polar;
  qd_oracle::ErrorBound bound =
      qd_oracle::complex_oracle::complex_special_bound("complex polar");
  Real r(1.75);
  Real theta(0.375);
  Complex input(r, theta);
  Complex got = polar(r, theta);
  mpc_t ref;
  mpfr_t r_mp;
  mpfr_t theta_mp;
  mpfr_t tmp;
  qd_oracle::complex_oracle::init_mpc<Real>(ref);
  mpfr_inits2(qd_oracle::complex_oracle::oracle_prec<Real>(), r_mp, theta_mp,
              tmp, (mpfr_ptr) 0);
  qd_oracle::to_mpfr(r_mp, r);
  qd_oracle::to_mpfr(theta_mp, theta);
  mpfr_cos(tmp, theta_mp, MPFR_RNDN);
  mpfr_mul(mpc_realref(ref), r_mp, tmp, MPFR_RNDN);
  mpfr_sin(tmp, theta_mp, MPFR_RNDN);
  mpfr_mul(mpc_imagref(ref), r_mp, tmp, MPFR_RNDN);
  double relerr = qd_oracle::complex_oracle::complex_relerr_in_eps(got, ref);
  bool pass = std::isfinite(relerr) && relerr <= bound.eps_multiplier;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = qd_oracle::complex_oracle::base_diag<Real>("test_complex_special",
                                                      "complex polar oracle", 0);
    qd_oracle::complex_oracle::append_complex_diag(&diag, input, got, ref,
                                                   relerr, bound);
  }
  tap.ok(pass, std::string(qd_oracle::TypeTraits<Real>::name()) + " complex polar oracle",
         diag, verbose);
  mpfr_clears(r_mp, theta_mp, tmp, (mpfr_ptr) 0);
  mpc_clear(ref);
  return pass;
}

template <class Real>
bool run_ceil(qd_oracle::Tap &tap, bool verbose) {
  typedef qd3_complex<Real> Complex;
  using std::ceil;
  qd_oracle::ErrorBound bound =
      qd_oracle::complex_oracle::complex_special_bound("complex ceil");
  Complex z(Real(1.25), Real(-0.5));
  Complex got = ceil(z);
  mpc_t ref;
  mpfr_t re;
  mpfr_t im;
  qd_oracle::complex_oracle::init_mpc<Real>(ref);
  mpfr_inits2(qd_oracle::complex_oracle::oracle_prec<Real>(), re, im,
              (mpfr_ptr) 0);
  qd_oracle::to_mpfr(re, z.real());
  qd_oracle::to_mpfr(im, z.imag());
  mpfr_ceil(mpc_realref(ref), re);
  mpfr_ceil(mpc_imagref(ref), im);
  double relerr = qd_oracle::complex_oracle::complex_relerr_in_eps(got, ref);
  bool pass = std::isfinite(relerr) && relerr <= bound.eps_multiplier;
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = qd_oracle::complex_oracle::base_diag<Real>("test_complex_special",
                                                      "complex ceil oracle", 0);
    qd_oracle::complex_oracle::append_complex_diag(&diag, z, got, ref,
                                                   relerr, bound);
  }
  tap.ok(pass, std::string(qd_oracle::TypeTraits<Real>::name()) + " complex ceil oracle",
         diag, verbose);
  mpfr_clears(re, im, (mpfr_ptr) 0);
  mpc_clear(ref);
  return pass;
}

template <class Real>
bool run_conj(qd_oracle::Tap &tap, bool verbose) {
  typedef qd3_complex<Real> Complex;
  Complex z(Real(1.25), Real(-0.5));
  Complex got = conj(z);
  bool pass = got == Complex(Real(1.25), Real(0.5));
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = qd_oracle::complex_oracle::base_diag<Real>("test_complex_special",
                                                      "complex conj exact", 0);
    diag.push_back(qd_oracle::Tap::Diagnostic(
        "input_limbs", qd_oracle::complex_oracle::complex_limbs_hex(z)));
    diag.push_back(qd_oracle::Tap::Diagnostic(
        "got_limbs", qd_oracle::complex_oracle::complex_limbs_hex(got)));
  }
  tap.ok(pass, std::string(qd_oracle::TypeTraits<Real>::name()) + " complex conj exact",
         diag, verbose);
  return pass;
}


template <class Real>
bool run_proj(qd_oracle::Tap &tap, bool verbose) {
  typedef qd3_complex<Real> Complex;
  using std::proj;
  Complex finite(Real(1.25), Real(-0.5));
  Complex got_finite = proj(finite);
  Complex infinite(Real::_inf, Real(-2.0));
  Complex got_infinite = proj(infinite);
  bool pass = got_finite == finite && got_infinite.real().isinf() &&
              got_infinite.imag() == Real(0.0) &&
              std::signbit(static_cast<long double>(
                  qd_oracle::TypeTraits<Real>::limb(got_infinite.imag(), 0)));
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = qd_oracle::complex_oracle::base_diag<Real>("test_complex_special",
                                                      "complex proj", 0);
    diag.push_back(qd_oracle::Tap::Diagnostic(
        "finite_input_limbs", qd_oracle::complex_oracle::complex_limbs_hex(finite)));
    diag.push_back(qd_oracle::Tap::Diagnostic(
        "finite_got_limbs", qd_oracle::complex_oracle::complex_limbs_hex(got_finite)));
    diag.push_back(qd_oracle::Tap::Diagnostic(
        "infinite_got_limbs", qd_oracle::complex_oracle::complex_limbs_hex(got_infinite)));
  }
  tap.ok(pass, std::string(qd_oracle::TypeTraits<Real>::name()) + " complex proj",
         diag, verbose);
  return pass;
}

template <class Real>
bool run_type(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  pass &= run_abs<Real>(tap, verbose);
  pass &= run_norm<Real>(tap, verbose);
  pass &= run_arg<Real>(tap, verbose);
  pass &= run_polar<Real>(tap, verbose);
  pass &= run_ceil<Real>(tap, verbose);
  pass &= run_conj<Real>(tap, verbose);
  pass &= run_proj<Real>(tap, verbose);
  return pass;
}

} // namespace

int main(int argc, char **argv) {
  qd_oracle::complex_oracle::Options options;
  if (!qd_oracle::complex_oracle::parse_options(argc, argv, "test_complex_special", &options)) {
    return 2;
  }
  qd_oracle::rng::configure(options.has_seed, options.seed);
  qd_oracle::Tap tap(qd_oracle::complex_oracle::selected_count(options) * kCasesPerType);
  std::cout << "# seed: " << qd_oracle::rng::active_seed() << "\n";

  bool pass = true;
  if (options.test_dd) pass &= run_type<dd_real>(tap, options.verbose);
  if (options.test_td) pass &= run_type<td_real>(tap, options.verbose);
  if (options.test_qd) pass &= run_type<qd_real>(tap, options.verbose);
#ifdef QD_HAVE_EDD_REAL
  if (options.test_edd) pass &= run_type<edd_real>(tap, options.verbose);
#endif
  return pass ? tap.exit_status() : 1;
}
