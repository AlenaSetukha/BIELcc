#ifndef _DOUBLE_SEG_CELL_INTEGRAL_SPH_H_
#define _DOUBLE_SEG_CELL_INTEGRAL_SPH_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>


#include "kernel_param.h"
#include "integral_param.h"
#include "integral_universal.h"
#include "integral_universal_seg.h"
#include "element_geom.h"
#include "kernel_lib.h"
#include "integrals_analytic.h"


namespace bielcc {

constexpr static const int _PMAX_SMOOTH_SEG = 8;

/**
    * @brief Double integral over the segment of the
    * semi-analytical surface integral over F(x-y) = e^{ikr} / r integrand.
    * @param a start of the integration interval
    * @param b end of the integration interval
    * @param cell cell of the inner surface integral
    * @param ker_param integrand parameters
    * @param int_param_seg segment integration parameters (partition, dimension, etc.)
    * @param int_param_cell surface integration parameters (partition, dimension, etc.)
    * @param res result
    * @note Rectangle formula in outer integral, inner integral - highlighting
    * singularity. No dividing by 4 pi.
*/
template <size_t CellPoints, typename KType>
std::complex<double>
DoubleSegmentCellIntegral_SimplePot_H_HS(const double* a, const double* b,
                                      const double (&cell)[CellPoints][3],
                                const KernelParam<KType>& ker_param_const,
                   double analyticR_r, const IntegralParam& int_param_seg,
                                      const IntegralParam& int_param_cell)
{
    KernelParam<KType> ker_param = ker_param_const;
    double d[3], xk[3], dl, delta = 0.;
    int n = int_param_seg.n_start;
    const int PMax = int_param_seg.GetPMax();
    const double epsOut = int_param_seg.GetEpsAccur();

    std::complex<double> ff = std::complex<double>(0., 0.);
    std::complex<double> res_prev = std::complex<double>(0., 0.);
    std::complex<double> res = std::complex<double>(0., 0.);

    KType k_wave = param.k;

    double start_rs = ker_param.smoothR;

    int p_n = 0;
    for (p_n = 0; p_n < PMax; p_n++) {
        for (int j = 0; j < 3; j++) {
            d[j] = (b[j] - a[j]) / n;
        }
        dl = vec_length(d);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                xk[j] = a[j] + (i + 0.5) * d[j];
            }
            // Расчет внутреннего интеграла
            ff = IntegralSurf_SimplePot_H_HS(cell, xk, analyticR_r, int_param_cell, k_wave);
            res += ff * static_cast<std::complex<double>>(dl);
        }

        delta = std::abs(res - res_prev) * std::abs(res - res_prev);

        if (delta <  int_param_seg.eps * int_param_seg.eps && p_n != 0) {
            break;
        }

        n = n * 2; // the step is reduced by 2 times -> smoothing radius divided by 2
        res_prev = res;
        res = std::complex<double>(0., 0.);
        if (p_n < _PMAX_SMOOTH_SEG && ker_param.smoothR > ker_param.smoothRMin) {
            ker_param.smoothR /= 2.;
        }
    }

    if (p_n == PMax) {
        res = res_prev;
    }
    
    return res;
}


}        // namespace bielcc
#endif   // _DOUBLE_SEG_CELL_INTEGRAL_SPH_H_