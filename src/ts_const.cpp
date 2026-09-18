#include "config.h"

#include <qd/ts_real.h>

template <> const single_real<3> single_real<3>::_nan(
    std::numeric_limits<float>::quiet_NaN());
template <> const single_real<3> single_real<3>::_inf(
    std::numeric_limits<float>::infinity());
template <> const single_real<3> single_real<3>::_2pi(
    "6.2831853071795864769252867665590057683943387987502");
template <> const single_real<3> single_real<3>::_pi(
    "3.1415926535897932384626433832795028841971693993751");
template <> const single_real<3> single_real<3>::_3pi4(
    "2.3561944901923449288469825374596271631478770495313");
template <> const single_real<3> single_real<3>::_pi2(
    "1.5707963267948966192313216916397514420985846996876");
template <> const single_real<3> single_real<3>::_pi4(
    "0.78539816339744830961566084581987572104929234984378");
template <> const single_real<3> single_real<3>::_e(
    "2.7182818284590452353602874713526624977572470937000");
template <> const single_real<3> single_real<3>::_log2(
    "0.69314718055994530941723212145817656807550013436026");
template <> const single_real<3> single_real<3>::_log10(
    "2.3025850929940456840179914546843642076011014886288");
template <> const single_real<3> single_real<3>::_max(
    std::numeric_limits<float>::max());
template <> const single_real<3> single_real<3>::_safe_max(
    std::numeric_limits<float>::max());

template <> const float single_real<3>::_eps = 0x1p-70f;
template <> const float single_real<3>::_min_normalized = 0x1p-78f;
template <> const int single_real<3>::_ndigits = 21;

bool ts_suppress_error_messages = false;
