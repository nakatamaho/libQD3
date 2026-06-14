#include "qd_rng.h"

#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace qd_oracle {
namespace rng {

namespace {

std::mt19937_64 global_engine;
std::uint64_t global_seed = default_seed;
bool configured = false;

bool parse_u64(const char *text, std::uint64_t *seed) {
  if (text == 0 || *text == '\0') {
    return false;
  }

  errno = 0;
  char *end = 0;
  unsigned long long value = std::strtoull(text, &end, 0);
  if (errno != 0 || end == text || *end != '\0') {
    return false;
  }

  *seed = static_cast<std::uint64_t>(value);
  return true;
}

} // namespace

bool parse_seed_arg(const char *arg, std::uint64_t *seed) {
  static const char prefix[] = "--seed=";
  if (arg == 0 || std::strncmp(arg, prefix, sizeof(prefix) - 1) != 0) {
    return false;
  }
  return parse_u64(arg + sizeof(prefix) - 1, seed);
}

void configure(bool has_cli_seed, std::uint64_t cli_seed) {
  const char *env_seed = std::getenv("QD_TEST_SEED");
  std::uint64_t selected = default_seed;

  if (env_seed != 0 && *env_seed != '\0') {
    if (!parse_u64(env_seed, &selected)) {
      std::cerr << "Invalid QD_TEST_SEED: " << env_seed << std::endl;
      std::exit(2);
    }
  } else if (has_cli_seed) {
    selected = cli_seed;
  }

  global_seed = selected;
  global_engine.seed(global_seed);
  configured = true;
}

std::mt19937_64 &engine() {
  if (!configured) {
    configure(false, default_seed);
  }
  return global_engine;
}

std::uint64_t active_seed() {
  if (!configured) {
    configure(false, default_seed);
  }
  return global_seed;
}

bool coin() {
  return (engine()() & 1U) != 0;
}

int uniform_int(int lo, int hi) {
  std::uniform_int_distribution<int> dist(lo, hi);
  return dist(engine());
}

std::uint64_t random_word52() {
  return engine()() & ((1ULL << 52) - 1ULL);
}

} // namespace rng
} // namespace qd_oracle
