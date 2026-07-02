#ifndef QD_TD_COMPLEX_H
#define QD_TD_COMPLEX_H

#include <qd/td_real.h>
#include <qd/detail/complex_impl.h>

using td_complex = qd3_complex<td_real>;

inline td_complex polar(const td_real &r, const td_real &theta) {
  return ::polar<td_real>(r, theta);
}

#endif /* QD_TD_COMPLEX_H */
