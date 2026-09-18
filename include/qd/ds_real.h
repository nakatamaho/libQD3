/* Public double-single real type. */
#ifndef _QD_DS_REAL_H
#define _QD_DS_REAL_H

#include <qd/single_real.h>

using ds_real = single_real<2>;

inline ds_real dsrand() { return ds_real::rand(); }

#endif /* _QD_DS_REAL_H */
