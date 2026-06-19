#include "mpfr_oracle.h"
#include "qd_rng.h"
#include "tap.h"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <qd/fpu.h>

namespace {

const int kSamples = 1024;

struct Options {
  bool test_dd;
  bool test_td;
  bool test_qd;
  bool verbose;
  bool has_seed;
  std::uint64_t seed;

  Options()
      : test_dd(false), test_td(false), test_qd(false), verbose(false),
        has_seed(false), seed(0) {}
};

void print_usage() {
  std::cout << "oracle_test_bf_commutativity [-dd] [-td] [-qd] [-all] [-v]"
            << " [--seed=N]\n";
}

template <class T>
bool run_case(qd_oracle::Tap &tap, bool verbose) {
  typedef qd_oracle::TypeTraits<T> traits;
  bool pass = true;
  int failed_iteration = -1;
  T failed_a;
  T failed_b;
  T failed_ab;
  T failed_ba;

#if defined(QD_BF_MUL)
  for (int i = 0; i < kSamples; ++i) {
    const T a = qd_oracle::rng::uniform_type<T>(-140, 140);
    const T b = qd_oracle::rng::uniform_type<T>(-140, 140);
    const T ab = a * b;
    const T ba = b * a;
    if (!qd_oracle::bit_equal(ab, ba)) {
      pass = false;
      failed_iteration = i;
      failed_a = a;
      failed_b = b;
      failed_ab = ab;
      failed_ba = ba;
      break;
    }
  }
#endif

  std::vector<qd_oracle::Tap::Diagnostic> diag;
  if (!pass || verbose) {
    diag.push_back(qd_oracle::Tap::Diagnostic(
        "seed", std::to_string(qd_oracle::rng::active_seed())));
    diag.push_back(qd_oracle::Tap::Diagnostic(
        "replay", std::string("tests/oracle/test_bf_commutativity -") +
                      traits::name() + " --seed=" +
                      std::to_string(qd_oracle::rng::active_seed())));
#if defined(QD_BF_MUL)
    diag.push_back(qd_oracle::Tap::Diagnostic(
        "iteration", std::to_string(failed_iteration)));
    diag.push_back(qd_oracle::Tap::Diagnostic("input_a_limbs",
                                              qd_oracle::limbs_hex(failed_a)));
    diag.push_back(qd_oracle::Tap::Diagnostic("input_b_limbs",
                                              qd_oracle::limbs_hex(failed_b)));
    diag.push_back(qd_oracle::Tap::Diagnostic("a_times_b_limbs",
                                              qd_oracle::limbs_hex(failed_ab)));
    diag.push_back(qd_oracle::Tap::Diagnostic("b_times_a_limbs",
                                              qd_oracle::limbs_hex(failed_ba)));
#else
    diag.push_back(qd_oracle::Tap::Diagnostic("skipped", "QD_BF_MUL not defined"));
#endif
  }

#if !defined(QD_BF_MUL)
  const std::string suffix = " skipped; QD_BF_MUL not defined";
#else
  const std::string suffix = " exact BF mul commutativity";
#endif
  tap.ok(pass, std::string(traits::name()) + suffix, diag, verbose);
  return pass;
}

Options parse_options(int argc, char **argv) {
  Options opts;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-dd") == 0) {
      opts.test_dd = true;
    } else if (std::strcmp(argv[i], "-td") == 0) {
      opts.test_td = true;
    } else if (std::strcmp(argv[i], "-qd") == 0) {
      opts.test_qd = true;
    } else if (std::strcmp(argv[i], "-all") == 0) {
      opts.test_dd = true;
      opts.test_td = true;
      opts.test_qd = true;
    } else if (std::strcmp(argv[i], "-v") == 0) {
      opts.verbose = true;
    } else if (qd_oracle::rng::parse_seed_arg(argv[i], &opts.seed)) {
      opts.has_seed = true;
    } else if (std::strcmp(argv[i], "-h") == 0 ||
               std::strcmp(argv[i], "--help") == 0) {
      print_usage();
      std::exit(0);
    } else {
      print_usage();
      std::exit(2);
    }
  }
  if (!opts.test_dd && !opts.test_td && !opts.test_qd) {
    opts.test_dd = true;
    opts.test_td = true;
    opts.test_qd = true;
  }
  return opts;
}

} // namespace

int main(int argc, char **argv) {
  fpu_fix_start(0);
  Options opts = parse_options(argc, argv);
  qd_oracle::rng::configure(opts.has_seed, opts.seed);

  int planned = 0;
  if (opts.test_dd) ++planned;
  if (opts.test_td) ++planned;
  if (opts.test_qd) ++planned;

  qd_oracle::Tap tap(planned);
  if (opts.test_dd) run_case<dd_real>(tap, opts.verbose);
  if (opts.test_td) run_case<td_real>(tap, opts.verbose);
  if (opts.test_qd) run_case<qd_real>(tap, opts.verbose);

  fpu_fix_end(0);
  return tap.exit_status();
}
