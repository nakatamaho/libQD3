#ifndef QD_DETAIL_COMPLEX_IMPL_H
#define QD_DETAIL_COMPLEX_IMPL_H

#include <algorithm>
#include <cmath>
#include <complex>

template <class Real>
class qd3_complex {
public:
  typedef Real value_type;

  qd3_complex() : re_(0.0), im_(0.0) {}
  qd3_complex(const Real &re) : re_(re), im_(0.0) {}
  qd3_complex(const Real &re, const Real &im) : re_(re), im_(im) {}
  qd3_complex(double re) : re_(re), im_(0.0) {}
  qd3_complex(double re, double im) : re_(re), im_(im) {}
  qd3_complex(const std::complex<double> &z) : re_(z.real()), im_(z.imag()) {}

  Real &real() { return re_; }
  const Real &real() const { return re_; }
  Real &imag() { return im_; }
  const Real &imag() const { return im_; }

  void real(const Real &re) { re_ = re; }
  void imag(const Real &im) { im_ = im; }

  qd3_complex &operator+=(const qd3_complex &rhs) {
    re_ += rhs.re_;
    im_ += rhs.im_;
    return *this;
  }

  qd3_complex &operator+=(const Real &rhs) {
    re_ += rhs;
    return *this;
  }

  qd3_complex &operator-=(const qd3_complex &rhs) {
    re_ -= rhs.re_;
    im_ -= rhs.im_;
    return *this;
  }

  qd3_complex &operator-=(const Real &rhs) {
    re_ -= rhs;
    return *this;
  }

  qd3_complex &operator*=(const qd3_complex &rhs) {
    *this = *this * rhs;
    return *this;
  }

  qd3_complex &operator*=(const Real &rhs) {
    re_ *= rhs;
    im_ *= rhs;
    return *this;
  }

  qd3_complex &operator/=(const qd3_complex &rhs) {
    *this = *this / rhs;
    return *this;
  }

  qd3_complex &operator/=(const Real &rhs) {
    re_ /= rhs;
    im_ /= rhs;
    return *this;
  }

  qd3_complex operator+() const { return *this; }
  qd3_complex operator-() const { return qd3_complex(-re_, -im_); }

  operator std::complex<double>() const {
    return std::complex<double>(to_double(re_), to_double(im_));
  }

private:
  Real re_;
  Real im_;
};

template <class Real>
inline Real real(const qd3_complex<Real> &z) {
  return z.real();
}

template <class Real>
inline Real imag(const qd3_complex<Real> &z) {
  return z.imag();
}

template <class Real>
inline bool operator==(const qd3_complex<Real> &a, const qd3_complex<Real> &b) {
  return a.real() == b.real() && a.imag() == b.imag();
}

template <class Real>
inline bool operator==(const qd3_complex<Real> &a, const Real &b) {
  return a.real() == b && a.imag() == Real(0.0);
}

template <class Real>
inline bool operator==(const Real &a, const qd3_complex<Real> &b) {
  return b == a;
}

template <class Real>
inline bool operator!=(const qd3_complex<Real> &a, const qd3_complex<Real> &b) {
  return !(a == b);
}

template <class Real>
inline bool operator!=(const qd3_complex<Real> &a, const Real &b) {
  return !(a == b);
}

template <class Real>
inline bool operator!=(const Real &a, const qd3_complex<Real> &b) {
  return !(a == b);
}

template <class Real>
inline bool operator==(const qd3_complex<Real> &a, double b) {
  return a == Real(b);
}

template <class Real>
inline bool operator==(double a, const qd3_complex<Real> &b) {
  return Real(a) == b;
}

template <class Real>
inline bool operator!=(const qd3_complex<Real> &a, double b) {
  return !(a == b);
}

template <class Real>
inline bool operator!=(double a, const qd3_complex<Real> &b) {
  return !(a == b);
}

template <class Real>
inline qd3_complex<Real> operator+(qd3_complex<Real> a, const qd3_complex<Real> &b) {
  a += b;
  return a;
}

template <class Real>
inline qd3_complex<Real> operator+(qd3_complex<Real> a, const Real &b) {
  a += b;
  return a;
}

template <class Real>
inline qd3_complex<Real> operator+(const Real &a, qd3_complex<Real> b) {
  b += a;
  return b;
}

template <class Real>
inline qd3_complex<Real> operator+(qd3_complex<Real> a, double b) {
  a += Real(b);
  return a;
}

template <class Real>
inline qd3_complex<Real> operator+(double a, qd3_complex<Real> b) {
  b += Real(a);
  return b;
}

template <class Real>
inline qd3_complex<Real> operator-(qd3_complex<Real> a, const qd3_complex<Real> &b) {
  a -= b;
  return a;
}

template <class Real>
inline qd3_complex<Real> operator-(qd3_complex<Real> a, const Real &b) {
  a -= b;
  return a;
}

template <class Real>
inline qd3_complex<Real> operator-(const Real &a, const qd3_complex<Real> &b) {
  return qd3_complex<Real>(a - b.real(), -b.imag());
}

template <class Real>
inline qd3_complex<Real> operator-(qd3_complex<Real> a, double b) {
  a -= Real(b);
  return a;
}

template <class Real>
inline qd3_complex<Real> operator-(double a, const qd3_complex<Real> &b) {
  return qd3_complex<Real>(Real(a) - b.real(), -b.imag());
}

template <class Real>
inline qd3_complex<Real> operator*(const qd3_complex<Real> &a,
                               const qd3_complex<Real> &b) {
  /* TODO: replace this with a scaled complex product before claiming
     overflow-safe multiplication for large-magnitude inputs. */
  return qd3_complex<Real>(a.real() * b.real() - a.imag() * b.imag(),
                       a.real() * b.imag() + a.imag() * b.real());
}

template <class Real>
inline qd3_complex<Real> operator*(qd3_complex<Real> a, const Real &b) {
  a *= b;
  return a;
}

template <class Real>
inline qd3_complex<Real> operator*(const Real &a, qd3_complex<Real> b) {
  b *= a;
  return b;
}

template <class Real>
inline qd3_complex<Real> operator*(qd3_complex<Real> a, double b) {
  a *= Real(b);
  return a;
}

template <class Real>
inline qd3_complex<Real> operator*(double a, qd3_complex<Real> b) {
  b *= Real(a);
  return b;
}

template <class Real>
inline qd3_complex<Real> operator/(const qd3_complex<Real> &a,
                               const qd3_complex<Real> &b) {
  using std::abs;
  const Real abr = abs(b.real());
  const Real abi = abs(b.imag());

  if (abr >= abi) {
    const Real r = b.imag() / b.real();
    const Real den = b.real() + b.imag() * r;
    return qd3_complex<Real>((a.real() + a.imag() * r) / den,
                         (a.imag() - a.real() * r) / den);
  }

  const Real r = b.real() / b.imag();
  const Real den = b.real() * r + b.imag();
  return qd3_complex<Real>((a.real() * r + a.imag()) / den,
                       (a.imag() * r - a.real()) / den);
}

template <class Real>
inline qd3_complex<Real> operator/(qd3_complex<Real> a, const Real &b) {
  a /= b;
  return a;
}

template <class Real>
inline qd3_complex<Real> operator/(const Real &a, const qd3_complex<Real> &b) {
  return qd3_complex<Real>(a) / b;
}

template <class Real>
inline qd3_complex<Real> operator/(qd3_complex<Real> a, double b) {
  a /= Real(b);
  return a;
}

template <class Real>
inline qd3_complex<Real> operator/(double a, const qd3_complex<Real> &b) {
  return qd3_complex<Real>(Real(a)) / b;
}

template <class Real>
inline qd3_complex<Real> operator+(qd3_complex<Real> a, const std::complex<double> &b) {
  return a + qd3_complex<Real>(b);
}

template <class Real>
inline qd3_complex<Real> operator+(const std::complex<double> &a, qd3_complex<Real> b) {
  return qd3_complex<Real>(a) + b;
}

template <class Real>
inline qd3_complex<Real> operator-(qd3_complex<Real> a, const std::complex<double> &b) {
  return a - qd3_complex<Real>(b);
}

template <class Real>
inline qd3_complex<Real> operator-(const std::complex<double> &a, const qd3_complex<Real> &b) {
  return qd3_complex<Real>(a) - b;
}

template <class Real>
inline qd3_complex<Real> operator*(qd3_complex<Real> a, const std::complex<double> &b) {
  return a * qd3_complex<Real>(b);
}

template <class Real>
inline qd3_complex<Real> operator*(const std::complex<double> &a, qd3_complex<Real> b) {
  return qd3_complex<Real>(a) * b;
}

template <class Real>
inline qd3_complex<Real> operator/(qd3_complex<Real> a, const std::complex<double> &b) {
  return a / qd3_complex<Real>(b);
}

template <class Real>
inline qd3_complex<Real> operator/(const std::complex<double> &a, const qd3_complex<Real> &b) {
  return qd3_complex<Real>(a) / b;
}

template <class Real>
inline qd3_complex<Real> conj(const qd3_complex<Real> &z) {
  return qd3_complex<Real>(z.real(), -z.imag());
}

template <class Real>
inline Real norm(const qd3_complex<Real> &z) {
  return z.real() * z.real() + z.imag() * z.imag();
}

template <class Real>
inline Real abs(const qd3_complex<Real> &z) {
  using std::abs;
  using std::sqrt;

  Real a = abs(z.real());
  Real b = abs(z.imag());
  if (a < b) {
    std::swap(a, b);
  }
  if (a == Real(0.0)) {
    return Real(0.0);
  }

  const Real t = b / a;
  return a * sqrt(Real(1.0) + t * t);
}

template <class Real>
inline Real arg(const qd3_complex<Real> &z) {
  using ::atan2;
  return atan2(z.imag(), z.real());
}

template <class Real>
inline qd3_complex<Real> polar(const Real &r, const Real &theta) {
  using std::cos;
  using std::sin;
  return qd3_complex<Real>(r * cos(theta), r * sin(theta));
}

template <class Real>
inline qd3_complex<Real> sqrt(const qd3_complex<Real> &z) {
  using std::sqrt;

  if (z.real() == Real(0.0) && z.imag() == Real(0.0)) {
    return qd3_complex<Real>();
  }

  const Real r = abs(z);
  if (z.real() >= Real(0.0)) {
    const Real t = sqrt((r + z.real()) / Real(2.0));
    return qd3_complex<Real>(t, z.imag() / (Real(2.0) * t));
  }

  Real t = sqrt((r - z.real()) / Real(2.0));
  if (std::signbit(to_double(z.imag()))) {
    t = -t;
  }
  return qd3_complex<Real>(z.imag() / (Real(2.0) * t), t);
}

template <class Real>
inline qd3_complex<Real> exp(const qd3_complex<Real> &z) {
  using std::cos;
  using std::exp;
  using std::sin;

  const Real er = exp(z.real());
  return qd3_complex<Real>(er * cos(z.imag()), er * sin(z.imag()));
}

template <class Real>
inline qd3_complex<Real> log(const qd3_complex<Real> &z) {
  using std::log;
  return qd3_complex<Real>(log(abs(z)), arg(z));
}

template <class Real>
inline qd3_complex<Real> log10(const qd3_complex<Real> &z) {
  using std::log;
  return log(z) / log(Real(10.0));
}

template <class Real>
inline qd3_complex<Real> log2(const qd3_complex<Real> &z) {
  using std::log;
  return log(z) / log(Real(2.0));
}

template <class Real>
inline qd3_complex<Real> sin(const qd3_complex<Real> &z) {
  using std::cos;
  using std::cosh;
  using std::sin;
  using std::sinh;

  return qd3_complex<Real>(sin(z.real()) * cosh(z.imag()),
                       cos(z.real()) * sinh(z.imag()));
}

template <class Real>
inline qd3_complex<Real> cos(const qd3_complex<Real> &z) {
  using std::cos;
  using std::cosh;
  using std::sin;
  using std::sinh;

  return qd3_complex<Real>(cos(z.real()) * cosh(z.imag()),
                       -sin(z.real()) * sinh(z.imag()));
}

template <class Real>
inline qd3_complex<Real> sinh(const qd3_complex<Real> &z) {
  using std::cos;
  using std::cosh;
  using std::sin;
  using std::sinh;

  return qd3_complex<Real>(sinh(z.real()) * cos(z.imag()),
                       cosh(z.real()) * sin(z.imag()));
}

template <class Real>
inline qd3_complex<Real> cosh(const qd3_complex<Real> &z) {
  using std::cos;
  using std::cosh;
  using std::sin;
  using std::sinh;

  return qd3_complex<Real>(cosh(z.real()) * cos(z.imag()),
                       sinh(z.real()) * sin(z.imag()));
}


template <class Real>
inline qd3_complex<Real> tan(const qd3_complex<Real> &z) {
  return sin(z) / cos(z);
}

template <class Real>
inline qd3_complex<Real> tanh(const qd3_complex<Real> &z) {
  return sinh(z) / cosh(z);
}

template <class Real>
inline qd3_complex<Real> asin(const qd3_complex<Real> &z) {
  const qd3_complex<Real> i(Real(0.0), Real(1.0));
  return -i * log(i * z + sqrt(qd3_complex<Real>(Real(1.0)) - z * z));
}

template <class Real>
inline qd3_complex<Real> acos(const qd3_complex<Real> &z) {
  return qd3_complex<Real>(Real::_pi2) - asin(z);
}

template <class Real>
inline qd3_complex<Real> atan(const qd3_complex<Real> &z) {
  const qd3_complex<Real> i(Real(0.0), Real(1.0));
  const qd3_complex<Real> one(Real(1.0));
  return (i * Real(0.5)) * (log(one - i * z) - log(one + i * z));
}

template <class Real>
inline qd3_complex<Real> asinh(const qd3_complex<Real> &z) {
  return log(z + sqrt(z * z + Real(1.0)));
}

template <class Real>
inline qd3_complex<Real> acosh(const qd3_complex<Real> &z) {
  return log(z + sqrt(z + Real(1.0)) * sqrt(z - Real(1.0)));
}

template <class Real>
inline qd3_complex<Real> atanh(const qd3_complex<Real> &z) {
  const qd3_complex<Real> one(Real(1.0));
  return Real(0.5) * (log(one + z) - log(one - z));
}

template <class Real>
inline qd3_complex<Real> pow(const qd3_complex<Real> &z, int n) {
  if (n == 0) {
    return qd3_complex<Real>(Real(1.0));
  }

  unsigned int exp = n < 0 ? static_cast<unsigned int>(-(n + 1)) + 1U
                           : static_cast<unsigned int>(n);
  qd3_complex<Real> base = z;
  qd3_complex<Real> result(Real(1.0));

  while (exp != 0U) {
    if ((exp & 1U) != 0U) {
      result *= base;
    }
    exp >>= 1U;
    if (exp != 0U) {
      base *= base;
    }
  }

  if (n < 0) {
    result = qd3_complex<Real>(Real(1.0)) / result;
  }
  return result;
}

template <class Real>
inline qd3_complex<Real> pow(const qd3_complex<Real> &z, const Real &a) {
  return exp(qd3_complex<Real>(a) * log(z));
}

template <class Real>
inline qd3_complex<Real> pow(const qd3_complex<Real> &z, double a) {
  return pow(z, Real(a));
}

template <class Real>
inline qd3_complex<Real> pow(const Real &a, const qd3_complex<Real> &z) {
  return exp(z * log(qd3_complex<Real>(a)));
}

template <class Real>
inline qd3_complex<Real> pow(double a, const qd3_complex<Real> &z) {
  return pow(Real(a), z);
}

template <class Real>
inline qd3_complex<Real> pow(const qd3_complex<Real> &a, const qd3_complex<Real> &b) {
  return exp(b * log(a));
}

template <class Real>
inline qd3_complex<Real> ldexp(const qd3_complex<Real> &z, int n) {
  using std::ldexp;
  return qd3_complex<Real>(ldexp(z.real(), n), ldexp(z.imag(), n));
}

template <class Real>
inline qd3_complex<Real> ceil(const qd3_complex<Real> &z) {
  using std::ceil;
  /* Component-wise extension; C++ does not define ceil(std::complex<T>). */
  return qd3_complex<Real>(ceil(z.real()), ceil(z.imag()));
}


template <class Real>
inline qd3_complex<Real> proj(const qd3_complex<Real> &z) {
  if (z.real().isinf() || z.imag().isinf()) {
    return qd3_complex<Real>(Real::_inf,
                             Real(std::copysign(0.0, to_double(z.imag()))));
  }
  return z;
}


#endif /* QD_DETAIL_COMPLEX_IMPL_H */
