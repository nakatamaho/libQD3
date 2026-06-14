#ifndef QD_ORACLE_TAP_H
#define QD_ORACLE_TAP_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace qd_oracle {

class Tap {
public:
  typedef std::pair<std::string, std::string> Diagnostic;

  explicit Tap(int planned, std::ostream &out = std::cout)
      : out_(out), planned_(planned), count_(0), failed_(0) {
    out_ << "TAP version 13\n";
    out_ << "1.." << planned_ << "\n";
  }

  void ok(bool pass, const std::string &name,
          const std::vector<Diagnostic> &diag = std::vector<Diagnostic>(),
          bool emit_diag_on_pass = false) {
    ++count_;
    if (!pass) {
      ++failed_;
    }

    const int display_digits = display_digits_for_test(name);
    out_ << (pass ? "ok " : "not ok ") << count_ << " - " << name << "\n";
    if ((!pass || emit_diag_on_pass) && !diag.empty()) {
      out_ << "  ---\n";
      for (std::vector<Diagnostic>::const_iterator it = diag.begin();
           it != diag.end(); ++it) {
        std::string value = it->second;
        if (should_shorten_value(it->first)) {
          value = trim_mantissa_digits(value, display_digits);
        }
        if (it->first == "mpfr_reference") {
          out_ << "  digit_ruler:    "
               << yaml_quote(digit_ruler(display_digits))
               << "\n";
        }
        out_ << "  " << it->first << ": ";
        if (it->first == "got_value") {
          out_ << "     ";
        }
        out_ << yaml_quote(value) << "\n";
      }
      out_ << "  ...\n";
    }
  }

  int exit_status() const {
    return failed_ == 0 && count_ == planned_ ? 0 : 1;
  }

private:
  static int display_digits_for_test(const std::string &name) {
    if (name.compare(0, 2, "dd") == 0) {
      return 33;
    }
    if (name.compare(0, 2, "td") == 0) {
      return 49;
    }
    if (name.compare(0, 2, "qd") == 0) {
      return 65;
    }
    if (name.compare(0, 3, "edd") == 0) {
      return 39;
    }
    return 80;
  }

  static std::string digit_ruler(int digits) {
    static const std::string base = "12345678901234567890";
    if (digits <= 0) {
      return "";
    }
    std::string ruler;
    ruler.reserve(digits);
    while (static_cast<int>(ruler.size()) < digits) {
      ruler += base;
    }
    ruler.resize(digits);
    return ruler;
  }

  static std::string yaml_quote(const std::string &text) {
    std::string quoted(1, static_cast<char>(39));
    for (std::string::const_iterator it = text.begin(); it != text.end(); ++it) {
      if (*it == 39) {
        quoted += static_cast<char>(39);
        quoted += static_cast<char>(39);
      } else {
        quoted += *it;
      }
    }
    quoted += static_cast<char>(39);
    return quoted;
  }

  static bool should_shorten_value(const std::string &name) {
    return name == "mpfr_reference" || name == "got_value" ||
           name == "input_value" || name == "input_a_value" ||
           name == "input_b_value" || name == "mpfr_sin" ||
           name == "mpfr_cos" || name == "mpfr_tan" ||
           name == "abs_error_mpfr" || name == "td_reduced_arg" ||
           name == "td_reduced_ref";
  }

  static std::string trim_mantissa_digits(const std::string &text,
                                         int max_significant_digits) {
    if (max_significant_digits <= 0 || text.empty()) {
      return text;
    }

    const std::string::size_type e_pos = text.find_first_of("eE");
    if (e_pos == std::string::npos) {
      return text;
    }

    const std::string::size_type dot_pos = text.find(static_cast<char>(46));
    if (dot_pos == std::string::npos || dot_pos >= e_pos) {
      return text;
    }

    std::size_t sign_pos = 0;
    if (!text.empty() && (text[0] == '-' || text[0] == '+')) {
      sign_pos = 1;
    }
    if (sign_pos >= e_pos) {
      return text;
    }

    const int current_digits = static_cast<int>(e_pos - sign_pos);
    if (current_digits <= max_significant_digits) {
      return text;
    }

    const int max_fraction_digits = max_significant_digits - 1;
    if (max_fraction_digits <= 0) {
      return text.substr(0, dot_pos) + text.substr(e_pos);
    }

    const std::string::size_type max_len =
        static_cast<std::string::size_type>(dot_pos + 1 + max_fraction_digits);
    if (max_len >= e_pos) {
      return text;
    }

    return text.substr(0, max_len) + text.substr(e_pos);
  }

  std::ostream &out_;
  int planned_;
  int count_;
  int failed_;
};

} // namespace qd_oracle

#endif
