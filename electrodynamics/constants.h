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
    const double ANALYTIC_CALC_DIST = 3.0;
    const double SMOOTHING_DIST = 0.1; // Для рисования поля
    const int START_CELL_SPLIT = 8;
    const int START_SEG_SPLIT = 10;
    const int PMAX_CELL_SPLIT = 10; // 2 ^{PMAX}
    const int PMAX_SEG_SPLIT = 10;  // 2^{PMAX}
    const double INTEGRAL_ACCURACY = 1e-6;
}
}       // namespace bielcc
#endif  // _EM_FIELD_CONSTANTS_H
