#ifndef QD_DD_COMPLEX_H
#define QD_DD_COMPLEX_H

#include <qd/dd_real.h>
#include <qd/detail/complex_impl.h>

using dd_complex = qd3_complex<dd_real>;

inline dd_complex polar(const dd_real &r, const dd_real &theta) {
  return ::polar<dd_real>(r, theta);
}

#endif /* QD_DD_COMPLEX_H */
