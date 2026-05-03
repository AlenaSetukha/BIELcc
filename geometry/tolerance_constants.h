#ifndef _TOLERANCE_CONSTANTS_H
#define _TOLERANCE_CONSTANTS_H

#include <complex>
#include <limits>
#include <iostream>

namespace bielcc {
    
namespace Tolerance_Constants {
    const double MACHINE_ZERO = std::numeric_limits<double>::epsilon();
    const double POINT_TOLERANCE = 1e-5;
    const double POINT_TOLERANCE_SQ = POINT_TOLERANCE * POINT_TOLERANCE;
}
}       // namespace bielcc
#endif  // _TOLERANCE_CONSTANTS_H