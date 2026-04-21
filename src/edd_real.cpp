/*
 * src/edd_real.cpp
 *
 * Non-inline implementation for edd_real.
 */
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <iostream>
#include <string>

#include "config.h"
#include <qd/edd_real.h>
#include "util.h"

#ifndef QD_INLINE
#include <qd/edd_inline.h>
#endif

using std::cerr;
using std::endl;
using std::ios_base;
using std::istream;
using std::ostream;
using std::string;

bool edd_suppress_error_messages = false;

namespace {

inline void round_string_edd(char *s, int precision, int *offset) {
  int i;

  if (precision > 0 && s[precision] >= '5') {
    s[precision - 1]++;
    i = precision - 1;
    while (i > 0 && s[i] > '9') {
      s[i] -= 10;
      s[--i]++;
    }
  }

  if (s[0] > '9') {
    for (i = precision; i >= 1; i--) {
      s[i + 1] = s[i];
    }
    s[0] = '1';
    s[1] = '0';
    (*offset)++;
    precision++;
  }

  s[precision] = 0;
}

inline edd_real sqrt_karp(const edd_real &a) {
  edd_word x = (edd_word) 1.0 / edd::sqrtx(a[0]);
  edd_word ax = a[0] * x;
  return edd_real::add(ax, (a - edd_real::sqr(ax))[0] * (x * (edd_word) 0.5));
}

inline bool edd_sqrt_needs_rescale(edd_word a_hi) {
  return edd::fabsx(a_hi) > edd::sqrt_rescale_thresh();
}

inline bool is_space_char(char ch) {
  return std::isspace(static_cast<unsigned char>(ch)) != 0;
}

inline std::string lower_ascii(const char *s) {
  std::string result(s);
  while (!result.empty() && is_space_char(result[result.size() - 1])) {
    result.erase(result.size() - 1);
  }
  for (std::string::size_type i = 0; i < result.size(); ++i) {
    result[i] = static_cast<char>(
        std::tolower(static_cast<unsigned char>(result[i])));
  }
  return result;
}

} // namespace

void edd_real::error(const char *msg) {
  if (edd_suppress_error_messages)
    return;
  if (msg)
    cerr << "ERROR " << msg << endl;
}

edd_real::edd_real(const char *s) {
  if (edd_real::read(s, *this)) {
    edd_real::error("(edd_real::edd_real): INPUT ERROR.");
    *this = edd_real::_nan;
  }
}

edd_real &edd_real::operator=(const char *s) {
  if (edd_real::read(s, *this)) {
    edd_real::error("(edd_real::operator=): INPUT ERROR.");
    *this = edd_real::_nan;
  }
  return *this;
}

edd_real sqrt(const edd_real &a) {
  if (a.is_zero())
    return (edd_word) 0.0;

  if (a.is_negative()) {
    edd_real::error("(edd_real::sqrt): Negative argument.");
    return edd_real::_nan;
  }

  if (edd_sqrt_needs_rescale(a[0])) {
    return mul_pwr2(sqrt_karp(mul_pwr2(a, (edd_word) 0.25)), (edd_word) 2.0);
  }

  return sqrt_karp(a);
}

edd_real edd_real::sqrt(edd_word a) {
  return ::sqrt(edd_real(a));
}

edd_real npwr(const edd_real &a, int n) {
  if (n == 0) {
    if (a.is_zero()) {
      edd_real::error("(edd_real::npwr): Invalid argument.");
      return edd_real::_nan;
    }
    return (edd_word) 1.0;
  }

  edd_real r = a;
  edd_real s((edd_word) 1.0);
  int N = std::abs(n);

  if (N > 1) {
    while (N > 0) {
      if (N % 2 == 1)
        s *= r;
      N /= 2;
      if (N > 0)
        r = sqr(r);
    }
  } else {
    s = r;
  }

  if (n < 0)
    return ((edd_word) 1.0) / s;

  return s;
}

edd_real pow(const edd_real &a, int n) {
  return npwr(a, n);
}

ostream &operator<<(ostream &os, const edd_real &edd_value) {
  bool showpos = (os.flags() & ios_base::showpos) != 0;
  bool uppercase = (os.flags() & ios_base::uppercase) != 0;
  return os << edd_value.to_string(os.precision(), os.width(), os.flags(),
      showpos, uppercase, os.fill());
}

istream &operator>>(istream &s, edd_real &a) {
  std::string str;
  s >> str;
  if (!s)
    return s;
  a = edd_real(str.c_str());
  return s;
}

int edd_real::read(const char *s, edd_real &a) {
  const char *p = s;
  char ch;
  int sign = 0;
  int point = -1;
  int nd = 0;
  int e = 0;
  bool done = false;
  edd_real r((edd_word) 0.0);

  while (is_space_char(*p))
    p++;

  const std::string lower = lower_ascii(p);
  if (lower == "nan" || lower == "+nan" || lower == "-nan") {
    a = edd_real::_nan;
    return 0;
  }
  if (lower == "inf" || lower == "+inf") {
    a = edd_real::_inf;
    return 0;
  }
  if (lower == "-inf") {
    a = -edd_real::_inf;
    return 0;
  }

  while (!done && (ch = *p) != '\0') {
    if (ch >= '0' && ch <= '9') {
      int d = ch - '0';
      r *= (edd_word) 10.0;
      r += (edd_word) d;
      nd++;
    } else {
      switch (ch) {
      case '.':
        if (point >= 0)
          return -1;
        point = nd;
        break;
      case '-':
      case '+':
        if (sign != 0 || nd > 0)
          return -1;
        sign = (ch == '-') ? -1 : 1;
        break;
      case 'E':
      case 'e':
        if (std::sscanf(p + 1, "%d", &e) != 1)
          return -1;
        done = true;
        break;
      case ' ':
      case '\t':
      case '\n':
      case '\r':
      case '\f':
      case '\v':
        done = true;
        break;
      default:
        return -1;
      }
    }
    p++;
  }

  if (point >= 0)
    e -= (nd - point);

  if (e != 0)
    r *= (edd_real((edd_word) 10.0) ^ e);

  a = (sign < 0) ? -r : r;
  return 0;
}

void edd_real::to_digits(char *s, int &expn, int precision) const {
  const int D = precision + 1;
  edd_real r = abs(*this);
  int e;
  int i;
  int d;

  if (x[0] == (edd_word) 0.0) {
    expn = 0;
    for (i = 0; i < precision; i++)
      s[i] = '0';
    s[precision] = 0;
    return;
  }

  e = static_cast<int>(edd::floorx(edd::log10x(edd::fabsx(x[0]))));

  if (e < -4000) {
    r *= edd_real((edd_word) 10.0) ^ 4000;
    r /= edd_real((edd_word) 10.0) ^ (e + 4000);
  } else if (e > 4000) {
    r = ldexp(r, -QD_EDD_FLT64X_MANT_DIG);
    r /= edd_real((edd_word) 10.0) ^ e;
    r = ldexp(r, QD_EDD_FLT64X_MANT_DIG);
  } else {
    r /= edd_real((edd_word) 10.0) ^ e;
  }

  if (r >= (edd_word) 10.0) {
    r /= (edd_word) 10.0;
    e++;
  } else if (r < (edd_word) 1.0) {
    r *= (edd_word) 10.0;
    e--;
  }

  if (r >= (edd_word) 10.0 || r < (edd_word) 1.0) {
    edd_real::error("(edd_real::to_digits): can't compute exponent.");
    return;
  }

  for (i = 0; i < D; i++) {
    d = static_cast<int>(r[0]);
    r -= (edd_word) d;
    r *= (edd_word) 10.0;
    s[i] = static_cast<char>(d + '0');
  }

  for (i = D - 1; i > 0; i--) {
    if (s[i] < '0') {
      s[i - 1]--;
      s[i] += 10;
    } else if (s[i] > '9') {
      s[i - 1]++;
      s[i] -= 10;
    }
  }

  if (s[0] <= '0') {
    edd_real::error("(edd_real::to_digits): non-positive leading digit.");
    return;
  }

  if (s[D - 1] >= '5') {
    s[D - 2]++;
    i = D - 2;
    while (i > 0 && s[i] > '9') {
      s[i] -= 10;
      s[--i]++;
    }
  }

  if (s[0] > '9') {
    e++;
    for (i = precision; i >= 2; i--) {
      s[i] = s[i - 1];
    }
    s[0] = '1';
    s[1] = '0';
  }

  s[precision] = 0;
  expn = e;
}

void edd_real::write(char *s, int len, int precision, bool showpos, bool uppercase) const {
  string str = to_string(precision, 0, ios_base::scientific, showpos, uppercase);
  std::strncpy(s, str.c_str(), len - 1);
  s[len - 1] = 0;
}

string edd_real::to_string(int precision, int width, ios_base::fmtflags fmt,
    bool showpos, bool uppercase, char fill) const {
  string s;
  bool fixed = (fmt & ios_base::fixed) != 0;
  bool sgn = true;
  int i;
  int e = 0;

  if (isinf()) {
    if (*this < (edd_word) 0.0) {
      s += '-';
    } else if (showpos) {
      s += '+';
    } else {
      sgn = false;
    }
    s += uppercase ? "INF" : "inf";
  } else if (isnan()) {
    s = uppercase ? "NAN" : "nan";
    sgn = false;
  } else {
    if (*this < (edd_word) 0.0) {
      s += '-';
    } else if (showpos) {
      s += '+';
    } else {
      sgn = false;
    }

    if (*this == (edd_word) 0.0) {
      s += '0';
      if (precision > 0) {
        s += '.';
        s.append(precision, '0');
      }
    } else {
      int off = fixed ?
          (1 + static_cast<int>(edd::floorx(edd::log10x(edd::fabsx(x[0]))))) : 1;
      int d = precision + off;
      int d_with_extra = fixed ? std::max(80, d) : d;

      if (fixed && precision == 0 && abs(*this) < (edd_word) 1.0) {
        s += (abs(*this) >= (edd_word) 0.5) ? '1' : '0';
        return s;
      }

      if (fixed && d <= 0) {
        s += '0';
        if (precision > 0) {
          s += '.';
          s.append(precision, '0');
        }
      } else {
        char *t;
        int j;

        if (fixed) {
          t = new char[d_with_extra + 1];
          to_digits(t, e, d_with_extra);
        } else {
          t = new char[d + 1];
          to_digits(t, e, d);
        }

        off = e + 1;

        if (fixed) {
          round_string_edd(t, d, &off);
          if (off > 0) {
            for (i = 0; i < off; i++)
              s += t[i];
            if (precision > 0) {
              s += '.';
              for (j = 0; j < precision; j++, i++)
                s += t[i];
            }
          } else {
            s += "0.";
            if (off < 0)
              s.append(-off, '0');
            for (i = 0; i < d; i++)
              s += t[i];
          }
        } else {
          s += t[0];
          if (precision > 0)
            s += '.';
          for (i = 1; i <= precision; i++)
            s += t[i];
        }

        delete [] t;
      }

      if (!fixed && !isinf()) {
        s += uppercase ? 'E' : 'e';
        append_expn(s, e);
      }
    }
  }

  const int len = static_cast<int>(s.length());
  if (len < width) {
    int delta = width - len;
    if (fmt & ios_base::internal) {
      if (sgn)
        s.insert(static_cast<string::size_type>(1), delta, fill);
      else
        s.insert(static_cast<string::size_type>(0), delta, fill);
    } else if (fmt & ios_base::left) {
      s.append(delta, fill);
    } else {
      s.insert(static_cast<string::size_type>(0), delta, fill);
    }
  }

  return s;
}
