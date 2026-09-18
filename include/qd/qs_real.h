/* Public quad-single real type. */
#ifndef _QD_QS_REAL_H
#define _QD_QS_REAL_H

#include <qd/single_real.h>

using qs_real = single_real<4>;

inline qs_real qsrand() { return qs_real::rand(); }

#endif /* _QD_QS_REAL_H */
