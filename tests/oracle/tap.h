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
          const std::vector<Diagnostic> &diag = std::vector<Diagnostic>()) {
    ++count_;
    if (!pass) {
      ++failed_;
    }

    out_ << (pass ? "ok " : "not ok ") << count_ << " - " << name << "\n";
    if (!pass && !diag.empty()) {
      out_ << "  ---\n";
      for (std::vector<Diagnostic>::const_iterator it = diag.begin();
           it != diag.end(); ++it) {
        out_ << "  " << it->first << ": " << yaml_quote(it->second) << "\n";
      }
      out_ << "  ...\n";
    }
  }

  int exit_status() const {
    return failed_ == 0 && count_ == planned_ ? 0 : 1;
  }

private:
  static std::string yaml_quote(const std::string &text) {
    std::string quoted = "'";
    for (std::string::const_iterator it = text.begin(); it != text.end(); ++it) {
      if (*it == '\'') {
        quoted += "''";
      } else {
        quoted += *it;
      }
    }
    quoted += "'";
    return quoted;
  }

  std::ostream &out_;
  int planned_;
  int count_;
  int failed_;
};

} // namespace qd_oracle

#endif
