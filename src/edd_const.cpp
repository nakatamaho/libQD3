/*
 * src/edd_const.cpp
 *
 * Constants for edd_real. The elementary constants are derived explicitly
 * from the existing qd_real reference values and then truncated into a
 * normalized two-limb binary80 expansion.
 */
#include "config.h"
#include <limits>
#include <qd/edd_real.h>

#ifndef QD_INLINE
#include <qd/edd_inline.h>
#endif

const edd_real edd_real::_2pi = to_edd_real(qd_real::_2pi);
const edd_real edd_real::_pi = to_edd_real(qd_real::_pi);
const edd_real edd_real::_3pi4 = to_edd_real(qd_real::_3pi4);
const edd_real edd_real::_pi2 = to_edd_real(qd_real::_pi2);
const edd_real edd_real::_pi4 = to_edd_real(qd_real::_pi4);
const edd_real edd_real::_e = to_edd_real(qd_real::_e);
const edd_real edd_real::_log2 = to_edd_real(qd_real::_log2);
const edd_real edd_real::_log10 = to_edd_real(qd_real::_log10);
const edd_real edd_real::_nan = edd_real(edd::d_nan(), (edd_word) 0.0);
const edd_real edd_real::_inf = edd_real(edd::d_inf(), (edd_word) 0.0);

const edd_word edd_real::_eps = edd::ldexpx((edd_word) 1.0, -126);
const edd_word edd_real::_min_normalized =
    edd::ldexpx((edd_word) 1.0, -16318);
const edd_real edd_real::_max = edd_real(
    std::numeric_limits<edd_word>::max(),
    edd::ldexpx(std::numeric_limits<edd_word>::max(), -QD_EDD_FLT64X_MANT_DIG));
const edd_real edd_real::_safe_max = edd_real(
    edd::ldexpx(std::numeric_limits<edd_word>::max(), -1),
    edd::ldexpx(std::numeric_limits<edd_word>::max(),
        -(QD_EDD_FLT64X_MANT_DIG + 1)));
const int edd_real::_ndigits = 38;
