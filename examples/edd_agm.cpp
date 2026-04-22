/*
 * examples/edd_agm.cpp
 *
 * Compute pi, log(2), log(10), log10(2), and log10(e)
 * with AGM-based formulas using edd_real.
 */
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include <qd/edd_real.h>

namespace {

edd_real agm(edd_real a, edd_real b, int max_iter, bool verbose,
    const char *label = NULL) {
  for (int i = 1; i <= max_iter; ++i) {
    edd_real a_next = (a + b) * (edd_word) 0.5;
    edd_real b_next = sqrt(a * b);

    if (verbose && label != NULL) {
      std::cout << label << " iteration " << i
                << ": a=" << a_next << " b=" << b_next << '\n';
    }

    a = a_next;
    b = b_next;

    if (abs(a - b) <= edd_real::_eps * abs(a)) {
      break;
    }
  }

  return (a + b) * (edd_word) 0.5;
}

edd_real compute_pi_agm(int max_iter, bool verbose) {
  edd_real a((edd_word) 1.0);
  edd_real b = sqrt(edd_real((edd_word) 0.5));
  edd_real t((edd_word) 0.25);
  edd_real p((edd_word) 1.0);

  if (verbose) {
    std::cout << "iteration 0: a=" << a << " b=" << b
              << " t=" << t << " p=" << p << '\n';
  }

  for (int i = 1; i <= max_iter; ++i) {
    edd_real a_next = (a + b) * (edd_word) 0.5;
    edd_real b_next = sqrt(a * b);
    edd_real delta = a - a_next;

    t -= p * sqr(delta);
    a = a_next;
    b = b_next;
    p *= (edd_word) 2.0;

    if (verbose) {
      std::cout << "iteration " << i << ": a=" << a << " b=" << b
                << " t=" << t << " p=" << p << '\n';
    }

    if (abs(a - b) <= edd_real::_eps * abs(a)) {
      break;
    }
  }

  return sqr(a + b) / (t * (edd_word) 4.0);
}

edd_real compute_log2_agm(int max_iter, int scale_shift, bool verbose) {
  edd_real one((edd_word) 1.0);
  edd_real b = ldexp(one, 1 - scale_shift);
  int agm_iters = (max_iter < 12) ? 12 : max_iter;
  edd_real g = agm(one, b, agm_iters, verbose, "log2 AGM");
  edd_real denom = g * (edd_word) (2 * (scale_shift + 1));
  return edd_real::_pi / denom;
}

edd_real compute_log_agm(const edd_real &x, const edd_real &log2_agm,
    int max_iter, int scale_shift, bool verbose) {
  edd_real one((edd_word) 1.0);
  edd_real scaled = ldexp(x, scale_shift);
  edd_real b = (edd_word) 4.0 / scaled;
  int agm_iters = (max_iter < 12) ? 12 : max_iter;
  edd_real g = agm(one, b, agm_iters, verbose, "log AGM");
  edd_real approx = edd_real::_pi / (g * (edd_word) 2.0);
  return approx - log2_agm * (edd_word) scale_shift;
}

int parse_int_arg(const char *s, int default_value) {
  if (s == NULL) {
    return default_value;
  }

  char *end = NULL;
  long value = std::strtol(s, &end, 10);
  if (end == s || *end != '\0' || value <= 0) {
    return default_value;
  }

  return static_cast<int>(value);
}

} // namespace

int main(int argc, char **argv) {
  int max_iter = 8;
  const int log2_shift = 80;
  bool verbose = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-v" || arg == "--verbose") {
      verbose = true;
    } else {
      max_iter = parse_int_arg(argv[i], max_iter);
    }
  }

  std::cout << std::setprecision(edd_real::_ndigits);

  edd_real pi_agm = compute_pi_agm(max_iter, verbose);
  edd_real log2_agm = compute_log2_agm(max_iter, log2_shift, verbose);
  edd_real log10_agm = compute_log_agm(edd_real((edd_word) 10.0),
      log2_agm, max_iter, log2_shift, verbose);
  edd_real log10_2_agm = log2_agm / log10_agm;
  edd_real log10_e_agm = edd_real((edd_word) 1.0) / log10_agm;
  edd_real err = abs(pi_agm - edd_real::_pi);
  edd_real log2_err = abs(log2_agm - edd_real::_log2);
  edd_real log10_err = abs(log10_agm - edd_real::_log10);
  edd_real log10_2_err = abs(log10_2_agm - (edd_real::_log2 / edd_real::_log10));
  edd_real log10_e_err = abs(log10_e_agm - ((edd_word) 1.0 / edd_real::_log10));

  std::cout << "pi (AGM)   = " << pi_agm << '\n';
  std::cout << "pi const   = " << edd_real::_pi << '\n';
  std::cout << "abs error  = " << err << '\n';
  std::cout << "error / eps = "
            << static_cast<long double>(to_float64x(err) / edd_real::_eps) << '\n';
  std::cout << '\n';
  std::cout << "log2 (AGM) = " << log2_agm << '\n';
  std::cout << "log2 const = " << edd_real::_log2 << '\n';
  std::cout << "abs error  = " << log2_err << '\n';
  std::cout << "error / eps = "
            << static_cast<long double>(to_float64x(log2_err) / edd_real::_eps) << '\n';
  std::cout << '\n';
  std::cout << "log(10) AGM = " << log10_agm << '\n';
  std::cout << "log(10) const = " << edd_real::_log10 << '\n';
  std::cout << "abs error    = " << log10_err << '\n';
  std::cout << "error / eps  = "
            << static_cast<long double>(to_float64x(log10_err) / edd_real::_eps) << '\n';
  std::cout << '\n';
  std::cout << "log10(2) AGM = " << log10_2_agm << '\n';
  std::cout << "log10(2) ref = " << (edd_real::_log2 / edd_real::_log10) << '\n';
  std::cout << "abs error    = " << log10_2_err << '\n';
  std::cout << "error / eps  = "
            << static_cast<long double>(to_float64x(log10_2_err) / edd_real::_eps) << '\n';
  std::cout << '\n';
  std::cout << "log10(e) AGM = " << log10_e_agm << '\n';
  std::cout << "log10(e) ref = " << ((edd_word) 1.0 / edd_real::_log10) << '\n';
  std::cout << "abs error    = " << log10_e_err << '\n';
  std::cout << "error / eps  = "
            << static_cast<long double>(to_float64x(log10_e_err) / edd_real::_eps) << '\n';

  return 0;
}
