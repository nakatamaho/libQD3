#include "mpc_complex_oracle.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

const int kSamples = 20;
const int kCasesPerType = 17;

enum UnaryCase {
  case_sqrt,
  case_exp,
  case_log,
  case_log10,
  case_log2,
  case_sin,
  case_cos,
  case_tan,
  case_asin,
  case_acos,
  case_atan,
  case_sinh,
  case_cosh,
  case_tanh,
  case_asinh,
  case_acosh,
  case_atanh
};

const char *case_name(UnaryCase op) {
  switch (op) {
  case case_sqrt: return "complex sqrt oracle";
  case case_exp: return "complex exp oracle";
  case case_log: return "complex log oracle";
  case case_log10: return "complex log10 oracle";
  case case_log2: return "complex log2 oracle";
  case case_sin: return "complex sin oracle";
  case case_cos: return "complex cos oracle";
  case case_tan: return "complex tan oracle";
  case case_asin: return "complex asin oracle";
  case case_acos: return "complex acos oracle";
  case case_atan: return "complex atan oracle";
  case case_sinh: return "complex sinh oracle";
  case case_cosh: return "complex cosh oracle";
  case case_tanh: return "complex tanh oracle";
  case case_asinh: return "complex asinh oracle";
  case case_acosh: return "complex acosh oracle";
  case case_atanh: return "complex atanh oracle";
  }
  return "complex unknown oracle";
}

template <class Real>
qd3_complex<Real> apply_lib(UnaryCase op, const qd3_complex<Real> &z) {
  using std::acos;
  using std::acosh;
  using std::asin;
  using std::asinh;
  using std::atan;
  using std::atanh;
  using std::cos;
  using std::cosh;
  using std::exp;
  using std::log;
  using std::log10;
  using std::log2;
  using std::sin;
  using std::sinh;
  using std::sqrt;
  using std::tan;
  using std::tanh;

  switch (op) {
  case case_sqrt: return sqrt(z);
  case case_exp: return exp(z);
  case case_log: return log(z);
  case case_log10: return log10(z);
  case case_log2: return log2(z);
  case case_sin: return sin(z);
  case case_cos: return cos(z);
  case case_tan: return tan(z);
  case case_asin: return asin(z);
  case case_acos: return acos(z);
  case case_atan: return atan(z);
  case case_sinh: return sinh(z);
  case case_cosh: return cosh(z);
  case case_tanh: return tanh(z);
  case case_asinh: return asinh(z);
  case case_acosh: return acosh(z);
  case case_atanh: return atanh(z);
  }
  return z;
}

template <class Real>
void apply_mpc(UnaryCase op, mpc_t out, mpc_t z) {
  switch (op) {
  case case_sqrt:
    mpc_sqrt(out, z, MPC_RNDNN);
    break;
  case case_exp:
    mpc_exp(out, z, MPC_RNDNN);
    break;
  case case_log:
    mpc_log(out, z, MPC_RNDNN);
    break;
  case case_log10:
    mpc_log10(out, z, MPC_RNDNN);
    break;
  case case_log2: {
    mpc_t denom;
    mpfr_t l2;
    mpfr_t zero;
    qd_oracle::complex_oracle::init_mpc<Real>(denom);
    mpfr_inits2(qd_oracle::complex_oracle::oracle_prec<Real>(), l2, zero,
                (mpfr_ptr) 0);
    mpfr_log_ui(l2, 2, MPFR_RNDN);
    mpfr_set_zero(zero, 0);
    mpc_log(out, z, MPC_RNDNN);
    mpc_set_fr_fr(denom, l2, zero, MPC_RNDNN);
    mpc_div(out, out, denom, MPC_RNDNN);
    mpfr_clears(l2, zero, (mpfr_ptr) 0);
    mpc_clear(denom);
    break;
  }
  case case_sin:
    mpc_sin(out, z, MPC_RNDNN);
    break;
  case case_cos:
    mpc_cos(out, z, MPC_RNDNN);
    break;
  case case_tan:
    mpc_tan(out, z, MPC_RNDNN);
    break;
  case case_asin:
    mpc_asin(out, z, MPC_RNDNN);
    break;
  case case_acos:
    mpc_acos(out, z, MPC_RNDNN);
    break;
  case case_atan:
    mpc_atan(out, z, MPC_RNDNN);
    break;
  case case_sinh:
    mpc_sinh(out, z, MPC_RNDNN);
    break;
  case case_cosh:
    mpc_cosh(out, z, MPC_RNDNN);
    break;
  case case_tanh:
    mpc_tanh(out, z, MPC_RNDNN);
    break;
  case case_asinh:
    mpc_asinh(out, z, MPC_RNDNN);
    break;
  case case_acosh:
    mpc_acosh(out, z, MPC_RNDNN);
    break;
  case case_atanh:
    mpc_atanh(out, z, MPC_RNDNN);
    break;
  }
}

template <class Real>
bool run_case(qd_oracle::Tap &tap, UnaryCase op, bool verbose) {
  typedef qd3_complex<Real> Complex;
  using namespace qd_oracle::complex_oracle;

  qd_oracle::ErrorBound bound = complex_transcendental_bound(case_name(op));
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  Complex worst_input;
  Complex worst_got;
  mpc_t input_mpc;
  mpc_t ref;
  mpc_t worst_ref;
  init_mpc<Real>(input_mpc);
  init_mpc<Real>(ref);
  init_mpc<Real>(worst_ref);

  for (int i = 0; i < kSamples; ++i) {
    Complex z = random_complex<Real>(-2, 1);
    if (op == case_log || op == case_log10 || op == case_log2) {
      z += Complex(Real(1.5), Real(0.25));
    }
    Complex got = apply_lib(op, z);
    to_mpc(input_mpc, z);
    apply_mpc<Real>(op, ref, input_mpc);
    double relerr = complex_relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr;
      worst_iteration = i;
      worst_input = z;
      worst_got = got;
      mpc_set(worst_ref, ref, MPC_RNDNN);
    }
    if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) {
      pass = false;
    }
  }

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = base_diag<Real>("test_complex_unary", case_name(op), worst_iteration);
    append_complex_diag(&diag, worst_input, worst_got, worst_ref, worst, bound);
  }
  std::string name = std::string(qd_oracle::TypeTraits<Real>::name()) + " " +
                     case_name(op);
  tap.ok(pass, name, diag, verbose);

  mpc_clear(input_mpc);
  mpc_clear(ref);
  mpc_clear(worst_ref);
  return pass;
}

template <class Real>
bool run_type(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  pass &= run_case<Real>(tap, case_sqrt, verbose);
  pass &= run_case<Real>(tap, case_exp, verbose);
  pass &= run_case<Real>(tap, case_log, verbose);
  pass &= run_case<Real>(tap, case_log10, verbose);
  pass &= run_case<Real>(tap, case_log2, verbose);
  pass &= run_case<Real>(tap, case_sin, verbose);
  pass &= run_case<Real>(tap, case_cos, verbose);
  pass &= run_case<Real>(tap, case_tan, verbose);
  pass &= run_case<Real>(tap, case_asin, verbose);
  pass &= run_case<Real>(tap, case_acos, verbose);
  pass &= run_case<Real>(tap, case_atan, verbose);
  pass &= run_case<Real>(tap, case_sinh, verbose);
  pass &= run_case<Real>(tap, case_cosh, verbose);
  pass &= run_case<Real>(tap, case_tanh, verbose);
  pass &= run_case<Real>(tap, case_asinh, verbose);
  pass &= run_case<Real>(tap, case_acosh, verbose);
  pass &= run_case<Real>(tap, case_atanh, verbose);
  return pass;
}

} // namespace

int main(int argc, char **argv) {
  qd_oracle::complex_oracle::Options options;
  if (!qd_oracle::complex_oracle::parse_options(argc, argv, "test_complex_unary", &options)) {
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
