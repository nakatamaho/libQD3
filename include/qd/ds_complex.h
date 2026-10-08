#ifndef QD_DS_COMPLEX_H
#define QD_DS_COMPLEX_H

#include <qd/ds_real.h>
#include <qd/detail/complex_impl.h>

using ds_complex = qd3_complex<ds_real>;

inline ds_complex polar(const ds_real &r, const ds_real &theta) {
  return ::polar<ds_real>(r, theta);
}

#endif /* QD_DS_COMPLEX_H */
