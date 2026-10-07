#ifndef QD_QS_COMPLEX_H
#define QD_QS_COMPLEX_H

#include <qd/qs_real.h>
#include <qd/detail/complex_impl.h>

using qs_complex = qd3_complex<qs_real>;

inline qs_complex polar(const qs_real &r, const qs_real &theta) {
  return ::polar<qs_real>(r, theta);
}

#endif /* QD_QS_COMPLEX_H */
