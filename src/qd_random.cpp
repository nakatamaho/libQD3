/*
 * src/qd_random.cpp
 *
 * xoshiro256** generator (Blackman and Vigna) seeded through splitmix64,
 * shared by all libQD3 random functions.  See include/qd/qd_random.h.
 */
#include <mutex>

#include <qd/qd_random.h>

namespace {

const uint64_t kDefaultSeed = 0x9E3779B97F4A7C15ULL;

inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

inline uint64_t splitmix64(uint64_t &x) {
  uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

struct Generator {
  uint64_t s[4];
  explicit Generator(uint64_t seed) { reseed(seed); }
  void reseed(uint64_t seed) {
    for (int i = 0; i < 4; ++i) s[i] = splitmix64(seed);
  }
  uint64_t next() {
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
  }
};

std::mutex &generator_mutex() {
  static std::mutex m;
  return m;
}

Generator &generator() {
  static Generator g(kDefaultSeed);
  return g;
}

} // namespace

extern "C" {

void qd_srand(uint64_t seed) {
  std::lock_guard<std::mutex> lock(generator_mutex());
  generator().reseed(seed);
}

uint64_t qd_rand_u64(void) {
  std::lock_guard<std::mutex> lock(generator_mutex());
  return generator().next();
}

} // extern "C"
