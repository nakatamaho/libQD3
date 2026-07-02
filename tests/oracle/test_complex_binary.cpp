#include "mpc_complex_oracle.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

const int kSamples = 20;
const int kCasesPerType = 5;

enum BinaryCase {
  case_pow_int,
  case_pow_real,
  case_real_pow,
  case_pow_complex,
  case_ldexp
};

const char *case_name(BinaryCase op) {
  switch (op) {
  case case_pow_int: return "complex pow int oracle";
  case case_pow_real: return "complex pow real oracle";
  case case_real_pow: return "real pow complex oracle";
  case case_pow_complex: return "complex pow complex oracle";
  case case_ldexp: return "complex ldexp oracle";
  }
  return "complex unknown oracle";
}

template <class Real>
void mpc_pow_int_ref(mpc_t out, mpc_t z, int n) {
  mpc_t base;
  mpc_t result;
  mpc_t one;
  qd_oracle::complex_oracle::init_mpc<Real>(base);
  qd_oracle::complex_oracle::init_mpc<Real>(result);
  qd_oracle::complex_oracle::init_mpc<Real>(one);

  mpc_set(base, z, MPC_RNDNN);
  mpc_set_ui(result, 1, MPC_RNDNN);
  mpc_set_ui(one, 1, MPC_RNDNN);

  unsigned int exp = n < 0 ? static_cast<unsigned int>(-(n + 1)) + 1U
                           : static_cast<unsigned int>(n);
  while (exp != 0U) {
    if ((exp & 1U) != 0U) {
      mpc_mul(result, result, base, MPC_RNDNN);
    }
    exp >>= 1U;
    if (exp != 0U) {
      mpc_mul(base, base, base, MPC_RNDNN);
    }
  }

  if (n < 0) {
    mpc_div(out, one, result, MPC_RNDNN);
  } else {
    mpc_set(out, result, MPC_RNDNN);
  }

  mpc_clear(base);
  mpc_clear(result);
  mpc_clear(one);
}

template <class Real>
bool run_case(qd_oracle::Tap &tap, BinaryCase op, bool verbose) {
  typedef qd3_complex<Real> Complex;
  using namespace qd_oracle::complex_oracle;
  using std::ldexp;
  using std::pow;

  qd_oracle::ErrorBound bound =
      (op == case_ldexp) ? qd_oracle::exact_ldexp_bound()
                         : complex_transcendental_bound(case_name(op));
  bool pass = true;
  double worst = -1.0;
  int worst_iteration = -1;
  Complex worst_a;
  Complex worst_b;
  Complex worst_got;
  mpc_t a_mpc;
  mpc_t b_mpc;
  mpc_t ref;
  mpc_t worst_ref;
  init_mpc<Real>(a_mpc);
  init_mpc<Real>(b_mpc);
  init_mpc<Real>(ref);
  init_mpc<Real>(worst_ref);

  for (int i = 0; i < kSamples; ++i) {
    Complex a = random_complex<Real>(-2, 1) + Complex(Real(1.25), Real(0.125));
    Complex b = random_complex<Real>(-3, 0);
    Complex got;
    to_mpc(a_mpc, a);
    to_mpc(b_mpc, b);

    if (op == case_pow_int) {
      int n = (i % 7) - 3;
      if (n == 0) n = 3;
      got = pow(a, n);
      mpc_pow_int_ref<Real>(ref, a_mpc, n);
      b = Complex(Real(n));
    } else if (op == case_pow_real) {
      Real e = qd_oracle::rng::uniform_type<Real>(-3, 0);
      got = pow(a, e);
      b = Complex(e);
      to_mpc(b_mpc, b);
      mpc_pow(ref, a_mpc, b_mpc, MPC_RNDNN);
    } else if (op == case_real_pow) {
      Real base = Real(1.5) + qd_oracle::rng::positive_type<Real>(-3, 0);
      got = pow(base, b);
      a = Complex(base);
      to_mpc(a_mpc, a);
      mpc_pow(ref, a_mpc, b_mpc, MPC_RNDNN);
    } else if (op == case_pow_complex) {
      got = pow(a, b);
      mpc_pow(ref, a_mpc, b_mpc, MPC_RNDNN);
    } else {
      int n = (i % 9) - 4;
      got = ldexp(a, n);
      mpc_set(ref, a_mpc, MPC_RNDNN);
      mpfr_mul_2si(mpc_realref(ref), mpc_realref(ref), n, MPFR_RNDN);
      mpfr_mul_2si(mpc_imagref(ref), mpc_imagref(ref), n, MPFR_RNDN);
      b = Complex(Real(n));
    }

    double relerr = complex_relerr_in_eps(got, ref);
    if (relerr > worst || !std::isfinite(relerr)) {
      worst = relerr;
      worst_iteration = i;
      worst_a = a;
      worst_b = b;
      worst_got = got;
      mpc_set(worst_ref, ref, MPC_RNDNN);
    }
    if (!std::isfinite(relerr) || relerr > bound.eps_multiplier) {
      pass = false;
    }
  }

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag = base_diag<Real>("test_complex_binary", case_name(op), worst_iteration);
    append_complex_binary_diag(&diag, worst_a, worst_b, worst_got, worst_ref,
                               worst, bound);
  }
  std::string name = std::string(qd_oracle::TypeTraits<Real>::name()) + " " +
                     case_name(op);
  tap.ok(pass, name, diag, verbose);

  mpc_clear(a_mpc);
  mpc_clear(b_mpc);
  mpc_clear(ref);
  mpc_clear(worst_ref);
  return pass;
}

template <class Real>
bool run_type(qd_oracle::Tap &tap, bool verbose) {
  bool pass = true;
  pass &= run_case<Real>(tap, case_pow_int, verbose);
  pass &= run_case<Real>(tap, case_pow_real, verbose);
  pass &= run_case<Real>(tap, case_real_pow, verbose);
  pass &= run_case<Real>(tap, case_pow_complex, verbose);
  pass &= run_case<Real>(tap, case_ldexp, verbose);
  return pass;
}

} // namespace

int main(int argc, char **argv) {
  qd_oracle::complex_oracle::Options options;
  if (!qd_oracle::complex_oracle::parse_options(argc, argv, "test_complex_binary", &options)) {
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
