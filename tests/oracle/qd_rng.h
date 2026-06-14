#ifndef QD_ORACLE_QD_RNG_H
#define QD_ORACLE_QD_RNG_H

#include "mpfr_oracle.h"

#include <cstdint>
#include <limits>
#include <random>

namespace qd_oracle {
namespace rng {

static const std::uint64_t default_seed = 0x9E3779B97F4A7C15ULL;

bool parse_seed_arg(const char *arg, std::uint64_t *seed);
void configure(bool has_cli_seed, std::uint64_t cli_seed);
std::mt19937_64 &engine();
std::uint64_t active_seed();
bool coin();
int uniform_int(int lo, int hi);
std::uint64_t random_word52();

template <class T>
void uniform_mpfr(mpfr_t out, int emin, int emax) {
  mpfr_set_prec(out, ref_prec<T>());
  mpfr_set_ui(out, 1, MPFR_RNDN);

  mpfr_t chunk;
  mpfr_init2(chunk, ref_prec<T>());
  const int bits = std::numeric_limits<T>::digits + 8;
  for (int offset = 0; offset < bits; offset += 52) {
    mpfr_set_ui(chunk, random_word52(), MPFR_RNDN);
    mpfr_mul_2si(chunk, chunk, -(offset + 52), MPFR_RNDN);
    mpfr_add(out, out, chunk, MPFR_RNDN);
  }

  mpfr_mul_2si(out, out, uniform_int(emin, emax), MPFR_RNDN);
  if (coin()) {
    mpfr_neg(out, out, MPFR_RNDN);
  }
  mpfr_clear(chunk);
}

template <class T>
T uniform_type(int emin, int emax) {
  mpfr_t value;
  mpfr_init2(value, ref_prec<T>());
  uniform_mpfr<T>(value, emin, emax);
  T result = from_mpfr_limbs<T>(value);
  mpfr_clear(value);
  return result;
}

template <class T>
T positive_type(int emin, int emax) {
  T value = uniform_type<T>(emin, emax);
  return value < T(0) ? -value : value;
}

template <class T>
T unit_type() {
  return uniform_type<T>(-12, -1);
}

} // namespace rng
} // namespace qd_oracle

#endif
