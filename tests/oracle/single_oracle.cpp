/* MPFR oracle coverage for ds_real, ts_real, and qs_real. */

#include "mpfr_oracle.h"
#include "qd_rng.h"
#include "tap.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <qd/ds_real.h>
#include <qd/fpu.h>
#include <qd/qs_real.h>
#include <qd/ts_real.h>

namespace {

struct Options {
  bool ds;
  bool ts;
  bool qs;
  bool verbose;
  bool has_seed;
  std::uint64_t seed;

  Options()
      : ds(false), ts(false), qs(false), verbose(false), has_seed(false),
        seed(0) {}
};

void usage() {
  std::cout << "oracle_test_single [-ds] [-ts] [-qs] [-all] [-v]"
            << " [--seed=N]\n";
}

template <class T>
double bound_for_arithmetic() {
  return qd_oracle::TypeTraits<T>::limbs == 2 ? 64.0
       : qd_oracle::TypeTraits<T>::limbs == 3 ? 256.0 : 1024.0;
}

template <class T>
double bound_for_transcendental() {
  return qd_oracle::TypeTraits<T>::limbs == 2 ? 2048.0
       : qd_oracle::TypeTraits<T>::limbs == 3 ? 8192.0 : 32768.0;
}

template <class T>
std::vector<qd_oracle::Tap::Diagnostic> failure_diag(
    const char *group, int iteration, const T &a, const T &b, const T &got,
    mpfr_t reference, double observed, double allowed) {
  std::vector<qd_oracle::Tap::Diagnostic> diag;
  std::ostringstream replay;
  replay << "tests/oracle/single_oracle -" << qd_oracle::TypeTraits<T>::name()
         << " --seed=" << qd_oracle::rng::active_seed();
  diag.push_back(qd_oracle::Tap::Diagnostic(
      "seed", std::to_string(qd_oracle::rng::active_seed())));
  diag.push_back(qd_oracle::Tap::Diagnostic("replay", replay.str()));
  diag.push_back(qd_oracle::Tap::Diagnostic("group", group));
  diag.push_back(qd_oracle::Tap::Diagnostic("iteration",
                                             std::to_string(iteration)));
  diag.push_back(qd_oracle::Tap::Diagnostic("input_a_limbs",
                                             qd_oracle::limbs_hex(a)));
  diag.push_back(qd_oracle::Tap::Diagnostic("input_a_value",
                                             qd_oracle::value_to_mpfr_string(a)));
  diag.push_back(qd_oracle::Tap::Diagnostic("input_b_limbs",
                                             qd_oracle::limbs_hex(b)));
  diag.push_back(qd_oracle::Tap::Diagnostic("input_b_value",
                                             qd_oracle::value_to_mpfr_string(b)));
  diag.push_back(qd_oracle::Tap::Diagnostic("mpfr_reference",
                                             qd_oracle::mpfr_to_string(reference)));
  diag.push_back(qd_oracle::Tap::Diagnostic("got_value",
                                             qd_oracle::value_to_mpfr_string(got)));
  diag.push_back(qd_oracle::Tap::Diagnostic("got_limbs",
                                             qd_oracle::limbs_hex(got)));
  diag.push_back(qd_oracle::Tap::Diagnostic(
      "abs_error_mpfr", qd_oracle::abs_error_to_string(got, reference)));
  diag.push_back(qd_oracle::Tap::Diagnostic("relerr_eps",
                                             std::to_string(observed)));
  diag.push_back(qd_oracle::Tap::Diagnostic("allowed_eps_multiplier",
                                             std::to_string(allowed)));
  return diag;
}

template <class T>
bool value_matches(const T &got, mpfr_t reference, double allowed) {
  if (mpfr_nan_p(reference)) return got.isnan();
  if (mpfr_inf_p(reference)) {
    return got.isinf() &&
        (std::signbit(static_cast<long double>(got.x[0])) ==
         (mpfr_sgn(reference) < 0));
  }
  return got.isfinite() && qd_oracle::relerr_in_eps(got, reference) <= allowed;
}

template <class T>
bool compare_sample(const char *group, int iteration, const T &a, const T &b,
                    const T &got, mpfr_t reference, double allowed,
                    std::vector<qd_oracle::Tap::Diagnostic> *failure) {
  const double observed = got.isfinite()
      ? qd_oracle::relerr_in_eps(got, reference) : HUGE_VAL;
  if (value_matches(got, reference, allowed)) return true;
  *failure = failure_diag(group, iteration, a, b, got, reference, observed,
                          allowed);
  return false;
}

template <class T>
bool arithmetic_group(std::vector<qd_oracle::Tap::Diagnostic> *failure) {
  const double allowed = bound_for_arithmetic<T>();
  mpfr_t a_mp;
  mpfr_t b_mp;
  mpfr_t reference;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, reference,
              (mpfr_ptr) 0);

  for (int i = 0; i < 64; ++i) {
    const T a = qd_oracle::rng::uniform_type<T>(-18, 18);
    T b = qd_oracle::rng::uniform_type<T>(-18, 18);
    if (b.is_zero()) b = T("1.125");
    qd_oracle::to_mpfr(a_mp, a);
    qd_oracle::to_mpfr(b_mp, b);

    mpfr_add(reference, a_mp, b_mp, MPFR_RNDN);
    if (!compare_sample("add", i, a, b, a + b, reference, allowed, failure))
      goto done;
    mpfr_sub(reference, a_mp, b_mp, MPFR_RNDN);
    if (!compare_sample("sub", i, a, b, a - b, reference, allowed, failure))
      goto done;
    mpfr_mul(reference, a_mp, b_mp, MPFR_RNDN);
    if (!compare_sample("mul", i, a, b, a * b, reference, allowed, failure))
      goto done;
    mpfr_div(reference, a_mp, b_mp, MPFR_RNDN);
    if (!compare_sample("div", i, a, b, a / b, reference, allowed, failure))
      goto done;
    mpfr_sqr(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("sqr", i, a, b, sqr(a), reference, allowed, failure))
      goto done;
  }

done:
  const bool pass = failure->empty();
  mpfr_clears(a_mp, b_mp, reference, (mpfr_ptr) 0);
  return pass;
}

template <class T>
bool algebraic_group(std::vector<qd_oracle::Tap::Diagnostic> *failure) {
  const double allowed = bound_for_arithmetic<T>() * 32.0;
  mpfr_t a_mp;
  mpfr_t b_mp;
  mpfr_t reference;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, reference,
              (mpfr_ptr) 0);

  for (int i = 0; i < 48; ++i) {
    const T positive = qd_oracle::rng::positive_type<T>(-8, 8);
    const T a = qd_oracle::rng::uniform_type<T>(-8, 8);
    T b = qd_oracle::rng::positive_type<T>(-6, 6);
    if (b.is_zero()) b = T("0.75");
    qd_oracle::to_mpfr(a_mp, a);
    qd_oracle::to_mpfr(b_mp, b);

    qd_oracle::to_mpfr(a_mp, positive);
    mpfr_sqrt(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("sqrt", i, positive, b, sqrt(positive), reference,
                       allowed, failure))
      goto done;
    mpfr_cbrt(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("nroot3", i, positive, b, nroot(positive, 3),
                       reference, allowed, failure))
      goto done;

    qd_oracle::to_mpfr(a_mp, a);
    mpfr_pow_si(reference, a_mp, 3, MPFR_RNDN);
    if (!compare_sample("pow_int", i, a, b, pow(a, 3), reference, allowed,
                       failure))
      goto done;
    mpfr_atan2(reference, a_mp, b_mp, MPFR_RNDN);
    if (!compare_sample("atan2", i, a, b, atan2(a, b), reference, allowed,
                       failure))
      goto done;
    mpfr_hypot(reference, a_mp, b_mp, MPFR_RNDN);
    if (!compare_sample("hypot", i, a, b, hypot(a, b), reference, allowed,
                       failure))
      goto done;
    mpfr_fmod(reference, a_mp, b_mp, MPFR_RNDN);
    if (!compare_sample("fmod", i, a, b, fmod(a, b), reference, allowed,
                       failure))
      goto done;
  }

done:
  const bool pass = failure->empty();
  mpfr_clears(a_mp, b_mp, reference, (mpfr_ptr) 0);
  return pass;
}

template <class T>
bool transcendental_group(std::vector<qd_oracle::Tap::Diagnostic> *failure) {
  const double allowed = bound_for_transcendental<T>();
  mpfr_t a_mp;
  mpfr_t b_mp;
  mpfr_t reference;
  mpfr_inits2(qd_oracle::ref_prec<T>(), a_mp, b_mp, reference,
              (mpfr_ptr) 0);

  for (int i = 0; i < 40; ++i) {
    const T a = qd_oracle::rng::uniform_type<T>(-2, 2);
    const T positive = qd_oracle::rng::positive_type<T>(-2, 2);
    const T unit = qd_oracle::rng::uniform_type<T>(-2, -1);
    qd_oracle::to_mpfr(a_mp, a);
    qd_oracle::to_mpfr(b_mp, positive);

    mpfr_exp(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("exp", i, a, positive, exp(a), reference, allowed,
                       failure))
      goto done;
    mpfr_expm1(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("expm1", i, a, positive, expm1(a), reference, allowed,
                       failure))
      goto done;
    mpfr_log(reference, b_mp, MPFR_RNDN);
    if (!compare_sample("log", i, positive, a, log(positive), reference,
                       allowed, failure))
      goto done;
    mpfr_log1p(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("log1p", i, a, positive, log1p(a), reference, allowed,
                       failure))
      goto done;
    mpfr_sin(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("sin", i, a, positive, sin(a), reference, allowed,
                       failure))
      goto done;
    mpfr_cos(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("cos", i, a, positive, cos(a), reference, allowed,
                       failure))
      goto done;
    mpfr_tanh(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("tanh", i, a, positive, tanh(a), reference, allowed,
                       failure))
      goto done;
    qd_oracle::to_mpfr(a_mp, unit);
    mpfr_asin(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("asin", i, unit, positive, asin(unit), reference,
                       allowed, failure))
      goto done;
    mpfr_acos(reference, a_mp, MPFR_RNDN);
    if (!compare_sample("acos", i, unit, positive, acos(unit), reference,
                       allowed, failure))
      goto done;
  }

done:
  const bool pass = failure->empty();
  mpfr_clears(a_mp, b_mp, reference, (mpfr_ptr) 0);
  return pass;
}

template <class T>
bool special_and_io_group(std::vector<qd_oracle::Tap::Diagnostic> *failure) {
  const double allowed = bound_for_arithmetic<T>() * 16.0;
  mpfr_t reference;
  mpfr_init2(reference, qd_oracle::ref_prec<T>());
  const char *inputs[] = {
      "1.23456789012345678901234567890123456789",
      "-9.87654321098765432109876543210987654321e-7",
      "3.1415926535897932384626433832795028841971"};
  for (int i = 0; i < 3; ++i) {
    T value(inputs[i]);
    mpfr_set_str(reference, inputs[i], 10, MPFR_RNDN);
    T zero;
    if (!compare_sample("decimal", i, value, zero, value, reference, allowed,
                       failure)) {
      mpfr_clear(reference);
      return false;
    }
    char text[256];
    value.write(text, sizeof(text), T::_ndigits);
    T parsed(text);
    if (!compare_sample("write_read", i, value, zero, parsed, reference,
                       allowed * 4.0, failure)) {
      mpfr_clear(reference);
      return false;
    }
  }

  const T nan = T::_nan;
  const T inf = T::_inf;
  const bool pass = nan.isnan() && inf.isinf() &&
      (T(1) / T(0)).isinf() && (inf - inf).isnan();
  if (!pass) {
    std::ostringstream replay;
    replay << "tests/oracle/single_oracle -"
           << qd_oracle::TypeTraits<T>::name() << " --seed="
           << qd_oracle::rng::active_seed();
    failure->push_back(qd_oracle::Tap::Diagnostic("replay", replay.str()));
  }
  mpfr_clear(reference);
  return pass && failure->empty();
}

template <class T>
void run_type(qd_oracle::Tap &tap, bool verbose) {
  std::vector<qd_oracle::Tap::Diagnostic> failure;
  bool pass = arithmetic_group<T>(&failure);
  tap.ok(pass, std::string(qd_oracle::TypeTraits<T>::name()) + " arithmetic",
         failure);
  if (!pass) return;

  failure.clear();
  pass = algebraic_group<T>(&failure);
  tap.ok(pass, std::string(qd_oracle::TypeTraits<T>::name()) + " algebraic",
         failure);
  if (!pass) return;

  failure.clear();
  pass = transcendental_group<T>(&failure);
  tap.ok(pass, std::string(qd_oracle::TypeTraits<T>::name()) +
             " transcendental",
         failure);
  if (!pass) return;

  failure.clear();
  pass = special_and_io_group<T>(&failure);
  tap.ok(pass, std::string(qd_oracle::TypeTraits<T>::name()) + " special_io",
         failure);
  if (verbose) {
    std::cout << "# " << qd_oracle::TypeTraits<T>::name()
              << " oracle groups passed=" << (pass ? "yes" : "no")
              << " seed=" << qd_oracle::rng::active_seed() << std::endl;
  }
}

} // namespace

int main(int argc, char **argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "-ds") options.ds = true;
    else if (arg == "-ts") options.ts = true;
    else if (arg == "-qs") options.qs = true;
    else if (arg == "-all") options.ds = options.ts = options.qs = true;
    else if (arg == "-v") options.verbose = true;
    else if (qd_oracle::rng::parse_seed_arg(argv[i], &options.seed))
      options.has_seed = true;
    else {
      usage();
      return 2;
    }
  }
  if (!options.ds && !options.ts && !options.qs) {
    options.ds = options.ts = options.qs = true;
  }
  qd_oracle::rng::configure(options.has_seed, options.seed);

  const int selected = (options.ds ? 1 : 0) + (options.ts ? 1 : 0) +
                       (options.qs ? 1 : 0);
  qd_oracle::Tap tap(selected * 4);
  unsigned int old_cw = 0;
  fpu_fix_start(&old_cw);
  const bool old_ds = ds_suppress_error_messages;
  const bool old_ts = ts_suppress_error_messages;
  const bool old_qs = qs_suppress_error_messages;
  ds_suppress_error_messages = true;
  ts_suppress_error_messages = true;
  qs_suppress_error_messages = true;
  if (options.ds) run_type<ds_real>(tap, options.verbose);
  if (options.ts) run_type<ts_real>(tap, options.verbose);
  if (options.qs) run_type<qs_real>(tap, options.verbose);
  ds_suppress_error_messages = old_ds;
  ts_suppress_error_messages = old_ts;
  qs_suppress_error_messages = old_qs;
  fpu_fix_end(&old_cw);
  return tap.exit_status();
}
