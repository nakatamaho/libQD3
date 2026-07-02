#ifndef QD_QD_COMPLEX_H
#define QD_QD_COMPLEX_H

#include <qd/qd_real.h>
#include <qd/detail/complex_impl.h>

using qd_complex = qd3_complex<qd_real>;

inline qd_complex polar(const qd_real &r, const qd_real &theta) {
  return ::polar<qd_real>(r, theta);
}

#endif /* QD_QD_COMPLEX_H */
