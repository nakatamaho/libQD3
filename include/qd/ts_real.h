/* Public triple-single real type. */
#ifndef _QD_TS_REAL_H
#define _QD_TS_REAL_H

#include <qd/single_real.h>

using ts_real = single_real<3>;

inline ts_real tsrand() { return ts_real::rand(); }

#endif /* _QD_TS_REAL_H */
