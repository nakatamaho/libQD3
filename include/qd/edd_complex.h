#ifndef QD_EDD_COMPLEX_H
#define QD_EDD_COMPLEX_H

#include <qd/qd_config.h>

#ifdef QD_HAVE_EDD_REAL
#include <qd/edd_real.h>
#include <qd/detail/complex_impl.h>

using edd_complex = qd3_complex<edd_real>;

inline edd_complex polar(const edd_real &r, const edd_real &theta) {
  return ::polar<edd_real>(r, theta);
}
#endif

#endif /* QD_EDD_COMPLEX_H */
