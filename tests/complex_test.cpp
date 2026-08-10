#include <qd/complex.h>
#include <qd/dd_complex.h>
#include <qd/td_complex.h>
#include <qd/qd_complex.h>
#include <qd/edd_complex.h>
#include <qd/fpu.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <string>

namespace {

template <class Real>
double real_eps() {
  return static_cast<double>(std::numeric_limits<Real>::epsilon());
}

template <class Real>
bool close_real(const Real &got, const Real &expected, double factor) {
  const double err = std::abs(to_double(got - expected));
  const double scale = std::max(1.0, std::abs(to_double(expected)));
  return err <= factor * real_eps<Real>() * scale;
}

template <class Complex>
bool close_complex(const Complex &got, const Complex &expected, double factor) {
  typedef typename Complex::value_type Real;
  return close_real(got.real(), expected.real(), factor) &&
         close_real(got.imag(), expected.imag(), factor);
}


inline bool real_signbit(const dd_real &x) {
  return std::signbit(x.x[0]);
}

inline bool real_signbit(const td_real &x) {
  return std::signbit(x[0]);
}

inline bool real_signbit(const qd_real &x) {
  return std::signbit(x[0]);
}

#ifdef QD_HAVE_EDD_REAL
inline bool real_signbit(const edd_real &x) {
  return std::signbit(static_cast<long double>(x[0]));
}
#endif

template <class Complex>
bool finite_complex(const Complex &z) {
  return std::isfinite(to_double(z.real())) &&
         std::isfinite(to_double(z.imag()));
}


template <class Real>
bool run_real_math_type(const char *name) {
  bool pass = true;

  using std::acosh;
  using std::asinh;
  using std::atanh;
  using std::cbrt;
  using std::exp;
  using std::exp2;
  using std::expm1;
  using std::fmod;
  using std::hypot;
  using std::log1p;
  using std::log2;
  using std::pow;
  using std::proj;
  using std::round;
  using std::sinh;
  using std::tan;
  using std::tanh;
  using std::tanh;
  using std::trunc;

  const Real small(0.03125);
  pass &= close_real(log2(Real(8.0)), Real(3.0), 4096.0);
  pass &= close_real(exp2(Real(3.0)), Real(8.0), 4096.0);
  pass &= close_real(log1p(expm1(small)), small, 4096.0);
  pass &= close_real(exp(log1p(small)), Real(1.0) + small, 4096.0);
  pass &= close_real(hypot(Real(3.0), Real(4.0)), Real(5.0), 4096.0);
  pass &= close_real(cbrt(Real(8.0)), Real(2.0), 4096.0);
  pass &= close_real(cbrt(Real(-8.0)), Real(-2.0), 4096.0);
  pass &= trunc(Real(-1.75)) == Real(-1.0);
  pass &= trunc(Real(1.75)) == Real(1.0);
  pass &= round(Real(-1.5)) == Real(-2.0);
  pass &= round(Real(1.5)) == Real(2.0);
  pass &= close_real(fmod(Real(5.5), Real(2.0)), Real(1.5), 4096.0);
  pass &= close_real(pow(Real(2.0), Real(3.0)), Real(8.0), 4096.0);
  pass &= close_real(sinh(asinh(Real(0.25))), Real(0.25), 4096.0);
  pass &= close_real(tanh(atanh(Real(0.25))), Real(0.25), 4096.0);
  pass &= close_real(cosh(acosh(Real(1.25))), Real(1.25), 4096.0);

  if (!pass) {
    std::cerr << name << " real math compatibility test failed\n";
  }
  return pass;
}

template <class Complex>
bool run_type(const char *name) {
  typedef typename Complex::value_type Real;
  bool pass = true;
  pass &= run_real_math_type<Real>(name);

  Complex z0;
  pass &= z0 == Real(0.0);

  Complex z1(Real(1.25));
  pass &= z1 == Real(1.25);

  Complex z2(Real(1.25), Real(-0.5));
  pass &= z2.real() == Real(1.25);
  pass &= z2.imag() == Real(-0.5);

  Complex z3(2.0, -3.0);
  pass &= z3.real() == Real(2.0);
  pass &= z3.imag() == Real(-3.0);

  std::complex<double> stdz(0.75, -0.125);
  Complex z4(stdz);
  std::complex<double> roundtrip = z4;
  pass &= roundtrip.real() == stdz.real();
  pass &= roundtrip.imag() == stdz.imag();

  z2.real() = Real(3.0);
  z2.imag() = Real(4.0);
  pass &= z2.real() == Real(3.0);
  pass &= z2.imag() == Real(4.0);
  z2.real(Real(5.0));
  z2.imag(Real(-6.0));
  const Complex &cz2 = z2;
  pass &= cz2.real() == Real(5.0);
  pass &= cz2.imag() == Real(-6.0);

  const Complex a(Real(1.0), Real(2.0));
  const Complex b(Real(3.0), Real(-4.0));
  pass &= (a + b) == Complex(Real(4.0), Real(-2.0));
  pass &= (a - b) == Complex(Real(-2.0), Real(6.0));
  pass &= (a * b) == Complex(Real(11.0), Real(2.0));
  pass &= close_complex(a / b,
                        Complex(-Real(1.0) / Real(5.0),
                                Real(2.0) / Real(5.0)),
                        64.0);

  Complex c = a;
  c += b;
  pass &= c == Complex(Real(4.0), Real(-2.0));
  c -= b;
  pass &= c == a;
  c *= b;
  pass &= c == Complex(Real(11.0), Real(2.0));
  c /= b;
  pass &= close_complex(c, a, 128.0);

  pass &= (a + Real(2.0)) == Complex(Real(3.0), Real(2.0));
  pass &= (Real(2.0) + a) == Complex(Real(3.0), Real(2.0));
  pass &= (a - Real(2.0)) == Complex(Real(-1.0), Real(2.0));
  pass &= (Real(2.0) - a) == Complex(Real(1.0), Real(-2.0));
  pass &= (a * Real(2.0)) == Complex(Real(2.0), Real(4.0));
  pass &= (Real(2.0) * a) == Complex(Real(2.0), Real(4.0));
  pass &= (a / Real(2.0)) == Complex(Real(0.5), Real(1.0));
  pass &= close_complex(Real(2.0) / a,
                        Complex(Real(2.0) / Real(5.0),
                                -Real(4.0) / Real(5.0)),
                        64.0);
  pass &= (a + 2.0) == Complex(Real(3.0), Real(2.0));
  pass &= (Complex(Real(2.0)) == Real(2.0));
  pass &= (Real(2.0) == Complex(Real(2.0)));
  pass &= (a != Real(1.0));

  using std::abs;
  using std::arg;
  using std::ceil;
  using std::acos;
  using std::acosh;
  using std::asin;
  using std::asinh;
  using std::atan;
  using std::atanh;
  using std::conj;
  using std::cos;
  using std::cosh;
  using std::exp;
  using std::imag;
  using std::ldexp;
  using std::log;
  using std::log10;
  using std::log2;
  using std::norm;
  using std::polar;
  using std::pow;
  using std::real;
  using std::sin;
  using std::sinh;
  using std::sqrt;

  Complex z(Real(1.25), Real(-0.5));
  Complex r_sqrt = sqrt(z);
  Complex r_exp = exp(z);
  Complex r_log = log(z);
  Complex r_log10 = log10(z);
  Complex r_log2 = log2(z);
  Complex r_sin = sin(z);
  Complex r_cos = cos(z);
  Complex r_tan = tan(z);
  Complex r_asin = asin(z);
  Complex r_acos = acos(z);
  Complex r_atan = atan(z);
  Complex r_sinh = sinh(z);
  Complex r_cosh = cosh(z);
  Complex r_tanh = tanh(z);
  Complex r_asinh = asinh(z);
  Complex r_acosh = acosh(z);
  Complex r_atanh = atanh(Complex(Real(0.25), Real(-0.125)));
  Complex r_pow = pow(z, z);
  Complex r_pow_i = pow(z, 3);
  Complex r_pow_r = pow(z, Real(2.0));
  Complex r_rpow = pow(Real(2.0), z);
  Complex r_ldexp = ldexp(z, 3);
  Complex r_ceil = ceil(z);
  Complex r_proj = proj(z);
  Complex r_proj_inf = proj(Complex(Real::_inf, Real(-2.0)));

  pass &= finite_complex(r_sqrt);
  pass &= finite_complex(r_exp);
  pass &= finite_complex(r_log);
  pass &= finite_complex(r_log10);
  pass &= finite_complex(r_log2);
  pass &= finite_complex(r_sin);
  pass &= finite_complex(r_cos);
  pass &= finite_complex(r_tan);
  pass &= finite_complex(r_asin);
  pass &= finite_complex(r_acos);
  pass &= finite_complex(r_atan);
  pass &= finite_complex(r_sinh);
  pass &= finite_complex(r_cosh);
  pass &= finite_complex(r_tanh);
  pass &= finite_complex(r_asinh);
  pass &= finite_complex(r_acosh);
  pass &= finite_complex(r_atanh);
  pass &= finite_complex(r_pow);
  pass &= finite_complex(r_pow_i);
  pass &= finite_complex(r_pow_r);
  pass &= finite_complex(r_rpow);

  bool identity_pass = true;
  identity_pass &= real(z) == z.real();
  identity_pass &= imag(z) == z.imag();
  identity_pass &= conj(conj(z)) == z;
  identity_pass &= real(conj(z)) == real(z);
  identity_pass &= imag(conj(z)) == -imag(z);
  identity_pass &= close_real(abs(z), sqrt(norm(z)), 128.0);
  identity_pass &= close_real(arg(z), atan2(imag(z), real(z)), 128.0);
  identity_pass &= close_complex(polar(Real(2.0), Real(0.25)),
                        Complex(Real(2.0) * cos(Real(0.25)),
                                Real(2.0) * sin(Real(0.25))),
                        128.0);
  identity_pass &= close_complex(r_sqrt * r_sqrt, z, 1024.0);
  identity_pass &= close_complex(exp(log(z)), z, 2048.0);
  identity_pass &= close_complex(r_sin * r_sin + r_cos * r_cos,
                        Complex(Real(1.0)), 4096.0);
  identity_pass &= close_complex(r_cosh * r_cosh - r_sinh * r_sinh,
                        Complex(Real(1.0)), 4096.0);
  identity_pass &= close_complex(r_tan, r_sin / r_cos, 4096.0);
  identity_pass &= close_complex(r_tanh, r_sinh / r_cosh, 4096.0);
  identity_pass &= close_complex(sin(r_asin), z, 8192.0);
  identity_pass &= close_complex(cos(r_acos), z, 8192.0);
  identity_pass &= close_complex(tan(r_atan), z, 8192.0);
  identity_pass &= close_complex(sinh(r_asinh), z, 8192.0);
  identity_pass &= close_complex(cosh(r_acosh), z, 8192.0);
  identity_pass &= close_complex(tanh(r_atanh),
                        Complex(Real(0.25), Real(-0.125)), 8192.0);
  identity_pass &= close_complex(r_pow_i, z * z * z, 1024.0);
  identity_pass &= close_complex(r_pow_r, z * z, 4096.0);
  identity_pass &= r_ldexp == Complex(Real(10.0), Real(-4.0));
  identity_pass &= r_ceil == Complex(Real(2.0), Real(0.0));
  identity_pass &= r_proj == z;
  identity_pass &= r_proj_inf.real().isinf();
  identity_pass &= r_proj_inf.imag() == Real(0.0);
  identity_pass &= real_signbit(r_proj_inf.imag());
  pass &= identity_pass;

  if (!pass) {
    std::cerr << name << " complex test failed\n";
  }
  return pass;
}

} // namespace

int main() {
  bool pass = true;
  {
    unsigned int old_cw = 0;
    fpu_fix_start(&old_cw);
    pass &= run_type<dd_complex>("dd_complex");
    pass &= run_type<td_complex>("td_complex");
    pass &= run_type<qd_complex>("qd_complex");
    fpu_fix_end(&old_cw);
  }
#ifdef QD_HAVE_EDD_REAL
  {
    unsigned int old_edd_cw = 0;
    fpu_fix_start_80bit(&old_edd_cw);
    pass &= run_type<edd_complex>("edd_complex");
    fpu_fix_end(&old_edd_cw);
  }
#endif

  if (pass) {
    std::cout << "complex_test passed\n";
    return 0;
  }
  return 1;
}
