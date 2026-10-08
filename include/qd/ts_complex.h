#ifndef QD_TS_COMPLEX_H
#define QD_TS_COMPLEX_H

#include <qd/ts_real.h>
#include <qd/detail/complex_impl.h>

using ts_complex = qd3_complex<ts_real>;

inline ts_complex polar(const ts_real &r, const ts_real &theta) {
  return ::polar<ts_real>(r, theta);
}

#endif /* QD_TS_COMPLEX_H */
