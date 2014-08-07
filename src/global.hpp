#ifndef GLOBAL_HPP
#define GLOBAL_HPP

#include <limits>

namespace PLBPD_NS{
  const unsigned int MAXPRMRYS = 10;
  const unsigned int NPRMRYDBLES = 7;
  const double PI = 3.1415926535897931;
  const double machineEpsilonD = std::numeric_limits<double>::epsilon();
}
#endif
