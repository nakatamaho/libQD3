#include "mpfr_oracle.h"

#include <cstdlib>
#include <iostream>

namespace qd_oracle {

void require_exact(int ternary, const char *operation) {
  if (ternary != 0) {
    std::cerr << "MPFR oracle conversion was inexact during "
              << operation << std::endl;
    std::abort();
  }
}

std::string mpfr_to_string(mpfr_t value) {
  char *buffer = 0;
  mpfr_asprintf(&buffer, "%.80Re", value);
  std::string result(buffer ? buffer : "");
  mpfr_free_str(buffer);
  return result;
}

} // namespace qd_oracle
