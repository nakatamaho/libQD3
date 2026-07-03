/*
 * src/edd_const.cpp
 *
 * Constants for edd_real as normalized two-limb binary80 expansions.
 */
#include "config.h"
#include <limits>
#include <qd/edd_real.h>

#ifndef QD_INLINE
#include <qd/edd_inline.h>
#endif

const edd_real edd_real::_2pi = edd_real(
    (edd_word) 0xc.90fdaa22168c235p-1L,
    (edd_word) -0xe.ce675d1fc8f8cbbp-67L);
const edd_real edd_real::_pi = edd_real(
    (edd_word) 0xc.90fdaa22168c235p-2L,
    (edd_word) -0xe.ce675d1fc8f8cbbp-68L);
const edd_real edd_real::_3pi4 = edd_real(
    (edd_word) 0x9.6cbe3f9990e91a8p-2L,
    (edd_word) -0xd.8d66c2ebeb5d4c6p-67L);
const edd_real edd_real::_pi2 = edd_real(
    (edd_word) 0xc.90fdaa22168c235p-3L,
    (edd_word) -0xe.ce675d1fc8f8cbbp-69L);
const edd_real edd_real::_pi4 = edd_real(
    (edd_word) 0xc.90fdaa22168c235p-4L,
    (edd_word) -0xe.ce675d1fc8f8cbbp-70L);
const edd_real edd_real::_e = edd_real(
    (edd_word) 0xa.df85458a2bb4a9bp-2L,
    (edd_word) -0xa.04753bfb185861cp-67L);
const edd_real edd_real::_log2 = edd_real(
    (edd_word) 0xb.17217f7d1cf79acp-4L,
    (edd_word) -0xd.871319ff0342543p-70L);
const edd_real edd_real::_log10 = edd_real(
    (edd_word) 0x9.35d8dddaaa8ac17p-2L,
    (edd_word) -0xa.d494ea3e967aeb9p-69L);
const edd_real edd_real::_nan = edd_real(edd::d_nan(), (edd_word) 0.0);
const edd_real edd_real::_inf = edd_real(edd::d_inf(), (edd_word) 0.0);

const edd_word edd_real::_eps = edd::ldexpx((edd_word) 1.0, -126);
const edd_word edd_real::_min_normalized =
    edd::ldexpx((edd_word) 1.0, -16318);
const edd_real edd_real::_max = edd_real(
    std::numeric_limits<edd_word>::max(),
    edd::ldexpx(std::numeric_limits<edd_word>::max(), -QD_EDD_WORD_MANT_DIG));
const edd_real edd_real::_safe_max = edd_real(
    edd::ldexpx(std::numeric_limits<edd_word>::max(), -1),
    edd::ldexpx(std::numeric_limits<edd_word>::max(),
        -(QD_EDD_WORD_MANT_DIG + 1)));
const int edd_real::_ndigits = 38;
