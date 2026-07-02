#include "mpc_complex_oracle.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

const int kSamples = 24;
const int kCasesPerType = 4;

enum ArithCase {
  case_add,
  case_sub,
  case_mul,
  case_div
};

const char *case_name(ArithCase op) {
  switch (op) {
  case case_add: return "complex add oracle";
  case case_sub: return "complex sub oracle";
  case case_mul: return "complex mul oracle";
  case case_div: return "complex div oracle";
  }
  return "complex unknown oracle";
}

template <class Real>
qd3_complex<Real> apply_lib(ArithCase op, const qd3_complex<Real> &a,
                             const qd3_complex<Real> &b) {
  switch (op) {
  case case_add: return a + b;
  case case_sub: return a - b;
  case case_mul: return a * b;
  case case_div: return a / b;
  }
  return a;
}

void apply_mpc(ArithCase op, mpc_t out, mpc_t a, mpc_t b) {
  switch (op) {
  case case_add: mpc_add(out, a, b, MPC_RNDNN); break;
  case case_sub: mpc_sub(out, a, b, MPC_RNDNN); break;
  case case_mul: mpc_mul(out, a, b, MPC_RNDNN); break;
  case case_div: mpc_div(out, a, b, MPC_RNDNN); break;
  }
}

template <class Real>
bool run_case(qd_oracle::Tap &tap, ArithCase op, bool verbose) {
  typedef qd3_complex<Real> Complex;
  using namespace qd_oracle::complex_oracle;

  qd_oracle::ErrorBound bound = complex_arithmetic_bound(case_name(op));
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
    Complex a = random_complex<Real>(-2, 2);
    Complex b = random_complex<Real>(-2, 2);
    if (op == case_div) {
      b += Complex(Real(1.25), Real(-0.75));
    }
    Complex got = apply_lib(op, a, b);
    to_mpc(a_mpc, a);
    to_mpc(b_mpc, b);
    apply_mpc(op, ref, a_mpc, b_mpc);
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
    diag = base_diag<Real>("test_complex_arith", case_name(op), worst_iteration);
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
  pass &= run_case<Real>(tap, case_add, verbose);
  pass &= run_case<Real>(tap, case_sub, verbose);
  pass &= run_case<Real>(tap, case_mul, verbose);
  pass &= run_case<Real>(tap, case_div, verbose);
  return pass;
}

} // namespace

int main(int argc, char **argv) {
  qd_oracle::complex_oracle::Options options;
  if (!qd_oracle::complex_oracle::parse_options(argc, argv, "test_complex_arith", &options)) {
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
