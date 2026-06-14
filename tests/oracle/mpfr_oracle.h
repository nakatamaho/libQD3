#ifndef QD_ORACLE_MPFR_ORACLE_H
#define QD_ORACLE_MPFR_ORACLE_H

#include <mpfr.h>

#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

#include <qd/dd_real.h>
#include <qd/qd_real.h>
#include <qd/td_real.h>
#ifdef QD_HAVE_EDD_REAL
#include <qd/edd_real.h>
#endif

namespace qd_oracle {

void require_exact(int ternary, const char *operation);
std::string mpfr_to_string(mpfr_t value);

template <class T> struct TypeTraits;

template <>
struct TypeTraits<dd_real> {
  typedef double limb_type;
  static const int limbs = 2;
  static const char *name() { return "dd"; }
  static limb_type limb(const dd_real &x, int i) { return x.x[i]; }
  static dd_real make(const limb_type *x) { return dd_real(x); }
  static int set_limb(mpfr_t out, limb_type x) {
    return mpfr_set_d(out, x, MPFR_RNDN);
  }
  static int set_eps(mpfr_t out) { return mpfr_set_d(out, dd_real::_eps, MPFR_RNDN); }
  static int set_min_normalized(mpfr_t out) {
    return mpfr_set_d(out, dd_real::_min_normalized, MPFR_RNDN);
  }
};

template <>
struct TypeTraits<td_real> {
  typedef double limb_type;
  static const int limbs = 3;
  static const char *name() { return "td"; }
  static limb_type limb(const td_real &x, int i) { return x.x[i]; }
  static td_real make(const limb_type *x) { return td_real(x); }
  static int set_limb(mpfr_t out, limb_type x) {
    return mpfr_set_d(out, x, MPFR_RNDN);
  }
  static int set_eps(mpfr_t out) { return mpfr_set_d(out, td_real::_eps, MPFR_RNDN); }
  static int set_min_normalized(mpfr_t out) {
    return mpfr_set_d(out, td_real::_min_normalized, MPFR_RNDN);
  }
};

template <>
struct TypeTraits<qd_real> {
  typedef double limb_type;
  static const int limbs = 4;
  static const char *name() { return "qd"; }
  static limb_type limb(const qd_real &x, int i) { return x.x[i]; }
  static qd_real make(const limb_type *x) { return qd_real(x); }
  static int set_limb(mpfr_t out, limb_type x) {
    return mpfr_set_d(out, x, MPFR_RNDN);
  }
  static int set_eps(mpfr_t out) { return mpfr_set_d(out, qd_real::_eps, MPFR_RNDN); }
  static int set_min_normalized(mpfr_t out) {
    return mpfr_set_d(out, qd_real::_min_normalized, MPFR_RNDN);
  }
};

#ifdef QD_HAVE_EDD_REAL
template <>
struct TypeTraits<edd_real> {
  typedef edd_word limb_type;
  static const int limbs = 2;
  static const char *name() { return "edd"; }
  static limb_type limb(const edd_real &x, int i) { return x.x[i]; }
  static edd_real make(const limb_type *x) { return edd_real(x); }
  static int set_limb(mpfr_t out, limb_type x) {
    return mpfr_set_ld(out, static_cast<long double>(x), MPFR_RNDN);
  }
  static int set_eps(mpfr_t out) {
    return mpfr_set_ld(out, static_cast<long double>(edd_real::_eps), MPFR_RNDN);
  }
  static int set_min_normalized(mpfr_t out) {
    return mpfr_set_ld(out, static_cast<long double>(edd_real::_min_normalized),
                       MPFR_RNDN);
  }
};
#endif

template <class T>
constexpr mpfr_prec_t ref_prec() {
  return static_cast<mpfr_prec_t>(std::numeric_limits<T>::digits + 160);
}

template <class T>
constexpr mpfr_prec_t acc_prec() {
  return static_cast<mpfr_prec_t>(std::numeric_limits<T>::digits + 160);
}

template <class T>
void to_mpfr(mpfr_t out, const T &x) {
  mpfr_set_prec(out, acc_prec<T>());
  mpfr_set_zero(out, 0);

  mpfr_t term;
  mpfr_init2(term, acc_prec<T>());

  for (int i = 0; i < TypeTraits<T>::limbs; ++i) {
    require_exact(TypeTraits<T>::set_limb(term, TypeTraits<T>::limb(x, i)),
                  "set limb");
    require_exact(mpfr_add(out, out, term, MPFR_RNDN), "add limb");
  }

  mpfr_clear(term);
}

template <class T>
T from_mpfr_limbs(mpfr_t in) {
  typedef TypeTraits<T> traits;
  typename traits::limb_type limbs[traits::limbs];

  mpfr_t residual;
  mpfr_t term;
  mpfr_init2(residual, acc_prec<T>());
  mpfr_init2(term, acc_prec<T>());
  mpfr_set(residual, in, MPFR_RNDN);

  for (int i = 0; i < traits::limbs; ++i) {
#ifdef QD_HAVE_EDD_REAL
    if (std::strcmp(traits::name(), "edd") == 0) {
      limbs[i] = static_cast<typename traits::limb_type>(
          mpfr_get_ld(residual, MPFR_RNDN));
    } else
#endif
    {
      limbs[i] = static_cast<typename traits::limb_type>(
          mpfr_get_d(residual, MPFR_RNDN));
    }
    require_exact(traits::set_limb(term, limbs[i]), "set extracted limb");
    mpfr_sub(residual, residual, term, MPFR_RNDN);
  }

  mpfr_clear(term);
  mpfr_clear(residual);
  return traits::make(limbs);
}

template <class T>
T round_nearest_expansion(mpfr_t exact) {
  return from_mpfr_limbs<T>(exact);
}

template <class Limb>
bool same_limb_value(Limb a, Limb b) {
  if (a != b) {
    return false;
  }
  if (a == static_cast<Limb>(0)) {
    return std::signbit(static_cast<long double>(a)) ==
           std::signbit(static_cast<long double>(b));
  }
  return true;
}

template <class T>
bool bit_equal(const T &a, const T &b) {
  for (int i = 0; i < TypeTraits<T>::limbs; ++i) {
    typename TypeTraits<T>::limb_type ax = TypeTraits<T>::limb(a, i);
    typename TypeTraits<T>::limb_type bx = TypeTraits<T>::limb(b, i);
    if (!same_limb_value(ax, bx)) {
      return false;
    }
  }
  return true;
}

template <class T>
double relerr_in_eps(const T &got, mpfr_t ref) {
  mpfr_t got_mp;
  mpfr_t diff;
  mpfr_t denom;
  mpfr_t eps;
  mpfr_t scaled;
  mpfr_inits2(ref_prec<T>(), got_mp, diff, denom, eps, scaled, (mpfr_ptr) 0);

  to_mpfr(got_mp, got);
  mpfr_sub(diff, got_mp, ref, MPFR_RNDN);
  mpfr_abs(diff, diff, MPFR_RNDN);

  mpfr_abs(denom, ref, MPFR_RNDN);
  require_exact(TypeTraits<T>::set_min_normalized(eps), "set min normalized");
  if (mpfr_zero_p(ref) || mpfr_cmp(denom, eps) < 0) {
    mpfr_set(denom, eps, MPFR_RNDN);
  }

  mpfr_div(scaled, diff, denom, MPFR_RNDN);
  require_exact(TypeTraits<T>::set_eps(eps), "set eps");
  mpfr_div(scaled, scaled, eps, MPFR_RNDN);

  double result = mpfr_get_d(scaled, MPFR_RNDN);
  mpfr_clears(got_mp, diff, denom, eps, scaled, (mpfr_ptr) 0);
  return result;
}

template <class T>
std::string value_to_mpfr_string(const T &x) {
  mpfr_t value;
  mpfr_init2(value, acc_prec<T>());
  to_mpfr(value, x);
  std::string result = mpfr_to_string(value);
  mpfr_clear(value);
  return result;
}

template <class T>
std::string abs_error_to_string(const T &got, mpfr_t ref) {
  mpfr_t got_mp;
  mpfr_t diff;
  mpfr_inits2(ref_prec<T>(), got_mp, diff, (mpfr_ptr) 0);
  to_mpfr(got_mp, got);
  mpfr_sub(diff, got_mp, ref, MPFR_RNDN);
  mpfr_abs(diff, diff, MPFR_RNDN);
  std::string result = mpfr_to_string(diff);
  mpfr_clears(got_mp, diff, (mpfr_ptr) 0);
  return result;
}

template <class T>
double ulp_error_estimate(const T &got, mpfr_t ref) {
  return relerr_in_eps(got, ref);
}

template <class Limb>
std::string hexfloat_fixed_precision(const Limb &value, int fraction_digits) {
  std::ostringstream text;
  text << std::showpoint << std::hexfloat << std::setprecision(fraction_digits)
       << value;
  std::string result = text.str();

  const std::string::size_type p_pos = result.find('p');
  if (p_pos == std::string::npos) {
    return result;
  }
  const std::string::size_type dot_pos = result.rfind('.', p_pos);
  if (dot_pos == std::string::npos) {
    return result;
  }

  const int digits = static_cast<int>(p_pos - dot_pos - 1);
  if (digits < fraction_digits) {
    result.insert(dot_pos + 1 + digits, fraction_digits - digits, '0');
  }
  return result;
}

template <class T>
std::string limbs_hex(const T &x) {
  std::ostringstream os;
  constexpr int limb_precision =
      (std::numeric_limits<typename TypeTraits<T>::limb_type>::digits - 1 + 3) / 4;
  os << "[";
  for (int i = 0; i < TypeTraits<T>::limbs; ++i) {
    if (i != 0) {
      os << ", ";
    }
    os << hexfloat_fixed_precision(TypeTraits<T>::limb(x, i), limb_precision);
  }
  os << "]";
  return os.str();
}

} // namespace qd_oracle

#endif
