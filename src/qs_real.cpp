/* Float quad-single template instantiation unit. */
#include <qd/qs_real.h>
#include <qd/qd_real.h>

#include <string>

void qd_single_detail::decimal_to_qd(const char *digits, int count,
                                     int scale10, double out[4]) {
  std::string text(digits, static_cast<std::size_t>(count));
  text += 'e';
  text += std::to_string(scale10);
  qd_real r;
  if (qd_real::read(text.c_str(), r) != 0) r = qd_real::_nan;
  for (int i = 0; i < 4; ++i) out[i] = r[i];
}

template struct single_real<4>;
