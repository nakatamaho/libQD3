#include "mpc_complex_oracle.h"
#include "qd_rng.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>

namespace qd_oracle {
namespace complex_oracle {

const char *rounding_name() {
  return "MPC_RNDNN";
}

std::string double_text(double value) {
  std::ostringstream os;
  os << value;
  return os.str();
}

std::string int_text(int value) {
  std::ostringstream os;
  os << value;
  return os.str();
}

void print_usage(const char *program) {
  std::cout << program << " [-dd] [-td] [-qd] [-edd] [-all] [-v] [--seed=N]\n";
}

bool parse_options(int argc, char **argv, const char *program, Options *options) {
  for (int i = 1; i < argc; ++i) {
    const char *arg = argv[i];
    std::uint64_t seed = 0;
    if (std::strcmp(arg, "-dd") == 0) {
      options->test_dd = true;
    } else if (std::strcmp(arg, "-td") == 0) {
      options->test_td = true;
    } else if (std::strcmp(arg, "-qd") == 0) {
      options->test_qd = true;
    } else if (std::strcmp(arg, "-edd") == 0) {
#ifdef QD_HAVE_EDD_REAL
      options->test_edd = true;
#else
      std::cerr << "edd_real is not available in this build\n";
      return false;
#endif
    } else if (std::strcmp(arg, "-all") == 0) {
      options->test_dd = true;
      options->test_td = true;
      options->test_qd = true;
#ifdef QD_HAVE_EDD_REAL
      options->test_edd = true;
#endif
    } else if (std::strcmp(arg, "-v") == 0 ||
               std::strcmp(arg, "-verbose") == 0) {
      options->verbose = true;
    } else if (qd_oracle::rng::parse_seed_arg(arg, &seed)) {
      options->has_seed = true;
      options->seed = seed;
    } else if (std::strcmp(arg, "-h") == 0 || std::strcmp(arg, "--help") == 0) {
      print_usage(program);
      std::exit(0);
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      print_usage(program);
      return false;
    }
  }

  if (!options->test_dd && !options->test_td && !options->test_qd
#ifdef QD_HAVE_EDD_REAL
      && !options->test_edd
#endif
  ) {
    options->test_dd = true;
    options->test_td = true;
    options->test_qd = true;
#ifdef QD_HAVE_EDD_REAL
    options->test_edd = true;
#endif
  }

  return true;
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

} // namespace complex_oracle
} // namespace qd_oracle
