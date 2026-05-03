#ifndef _EM_FIELD_CONSTANTS_H
#define _EM_FIELD_CONSTANTS_H

#include <complex>
#include <limits>
#include <iostream>

namespace bielcc {

namespace ED_Constants {
    constexpr double ra = 57.295779513;
    constexpr double eps0 = 8.8541878128 * 1e-12;
    constexpr double mu0 = 4 * M_PI * 1e-7;
    constexpr int c_vacuum = 299792458; // м/с
}

namespace Calculation_Constants {
    const double MACHINE_ZERO = std::numeric_limits<double>::epsilon();
    const double INTEGRAL_ACCURACY = 1e-6;            // calculation integralsy
    const double ANALYTIC_CALC_DIST = 3.0;            // analytical calc R (relative to grid step h)
    const double SMOOTHING_DIST_SURF_INT = 0.1;        // surface smoothing radius (relative to small cell step h2
    const double SMOOTHING_DIST_SEG_INT = 0.1;         // segment smoothing radius (relative to small cell step h2
    const int START_CELL_SPLIT = 20;                  // starting cell split
    const int START_SEG_SPLIT = 40;                   // starting segment split
    const int PMAX_CELL_SPLIT = 2;                    // max cell split index (2 ^{PMAX})
    const int PMAX_SEG_SPLIT = 2;                     // max segment split index (2^{PMAX})
}
}       // namespace bielcc
#endif  // _EM_FIELD_CONSTANTS_H
