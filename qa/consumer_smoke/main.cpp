#include <iostream>
#include <qd/qd_real.h>

int main() {
  qd_real x = sqrt(qd_real(2.0));
  qd_real y = x * x - qd_real(2.0);
  std::cout << y << "\n";
  return 0;
}
