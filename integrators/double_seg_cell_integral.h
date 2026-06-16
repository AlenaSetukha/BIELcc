#ifndef _DOUBLE_SEG_CELL_INTEGRAL_H_
#define _DOUBLE_SEG_CELL_INTEGRAL_H_

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

constexpr static const int _DIM_MAX_SEG = 3;
constexpr static const int _PMAX_SMOOTH_SEG = 8;

/**
    * @brief Double integral over the segment of the
    * cell integral. Completely numerical.
    * @param a start of the integration interval
    * @param b end of the integration interval
    * @param cell cell of the inner surface integral
    * @param f_0 integrand F(x - y)
    * @param ker_param integrand parameters
    * @param int_param_seg segment integration parameters (partition, dimension, etc.)
    * @param int_param_cell surface integration parameters (partition, dimension, etc.)
    * @param res result
    * @note Rectangle formula
*/
template <size_t CellPoints, typename KType, typename P>
void DoubleSegmentCellIntegral(const double* a, const double* b,
                             const double (&cell)[CellPoints][3],
        void (*f_0)(const double*, const double*, const KernelParam<KType>&, P*),
               KernelParam<KType>& ker_param, const IntegralParam& int_param_seg,
                                     const IntegralParam& int_param_cell, P* res)
{
    double d[3], xk[3], dl, delta = 0.;
    int n = int_param_seg.GetNStart();
    int idim = int_param_seg.GetIDim();
    double eps_acc = int_param_seg.GetEpsAccur();
    const int PMax = int_param_seg.GetPMax();

    double start_rs = ker_param.smoothR;

    P ff[_DIM_MAX_SEG]{}, res_prev[_DIM_MAX_SEG]{}; 
    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
    }

    int p_n = 0;
    for (p_n = 0; p_n < PMax; p_n++) {
        for (int j = 0; j < 3; j++) {
            d[j] = (b[j] - a[j]) / n;
        }
        dl = vec_length(d); // шаг при текущем n

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                xk[j] = a[j] + (i + 0.5) * d[j]; // середина k-го отрезка разбиения
            }
            
            // Расчет внутреннего интеграла
            IntegralUniversalPnt(xk, cell, f_0, ker_param, int_param_cell, ff);

            for (int g = 0; g < idim; g++) {
                res[g] += ff[g] * static_cast<P>(dl);
            }
        }

        delta = 0.;
        for (int g = 0; g < idim; g++) {
            delta += std::abs(res[g] - res_prev[g]) * std::abs(res[g] - res_prev[g]);
        }

        if (delta < eps_acc * eps_acc && p_n != 0) {
            break;
        }

        n = n * 2; // the step is reduced by 2 times -> smoothing radius divided by 2
        for (int g = 0; g < idim; g++) {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
            if (p_n < _PMAX_SMOOTH_SEG && ker_param.smoothR > ker_param.smoothRMin) {
                ker_param.smoothR /= 2.;
            }
        }
    }

    if (p_n == PMax) {
        for (int g = 0; g < idim; g++) {
            res[g] = res_prev[g];
        }
    }

    ker_param.smoothR = start_rs;
}

}        // namespace bielcc
#endif   // _DOUBLE_SEG_CELL_INTEGRAL_H_