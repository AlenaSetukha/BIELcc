#ifndef _INTEGRAL_UNIVERSAL_SEG_H_
#define _INTEGRAL_UNIVERSAL_SEG_H_

#include <iostream>

#include "kernel_param.h"
#include "integral_param.h"
#include "element_geom.h"


//===========================================================
//------Integral in point over a segment by F(xk, y)---------
//===========================================================
/**
    * Calculates the curvilinear integral over a segment of
    * the function F(xk,y) at the fixed point xk.
    * Rectangle formula, integration over dy.
 */

constexpr static const int _DIM_MAX = 3;


namespace bielcc {
/**
    * @brief Integral over a segment [a, b] from the function F(xk, y)
    * @param a Point1 [3]
    * @param b Point2 [3]
    * @param x Collocation point F(x, y) (fixed)
    * @param f_0 Integrand
    * @param ker_param Integrand parameters
    * @param int_param Integration parameters
    * @param res Result
    * @note Rectangle formula, integration by dy
*/
template <typename P, typename KType>
void IntegralUniversalSegPnt(const double* a, const double* b, const double* x,
        void (*f_0)(const double*, const double*, const KernelParam<KType>&, P*),
        const KernelParam<KType>& ker_param, const IntegralParam& int_param,
        P* res)
{
    int p_n;
    int n = int_param.GetNStart();
    double d[3], y[3], dl, delta = 0.;
    int idim = int_param.GetIDim();


    P ff[_DIM_MAX]{}, res_prev[_DIM_MAX]{}; 
    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
    }


    for (p_n = 0; p_n < int_param.GetPMax(); p_n++) {
        for (int j = 0; j < 3; j++) {
            d[j] = (b[j] - a[j]) / n;
        }
        dl = vec_length(d);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                y[j] = a[j] + (i + 0.5) * d[j];
            }
            f_0(x, y, ker_param, ff);

            for (int g = 0; g < idim; g++) {
                res[g] += ff[g] * static_cast<P>(dl);
            }
        }

        delta = 0.;
        for (int g = 0; g < idim; g++) {
            delta += std::abs(res[g] - res_prev[g]) * std::abs(res[g] - res_prev[g]);
        }

        if (delta <  int_param.GetEpsAccur() * int_param.GetEpsAccur() && p_n != 0) {
            break;
        }

        n = n * 2;
        for (int g = 0; g < idim; g++) {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
        }
    }

    if (p_n == int_param.GetPMax()) {
        for (int g = 0; g < idim; g++) {
            res[g] = res_prev[g];
        }
    }
}



//===========================================================
//------------Integral over a segment by F(x)------------
//===========================================================
/**
    * Calculates the curvilinear integral over a segment of
    * the function F(x). Rectangle formula, integration over dx.
*/

/**
    * @brief Integral over a segment [a, b] from the function F(x)
    * @param a Point1 [3]
    * @param b Point2 [3]
    * @param f_0 Integrand
    * @param ker_param Integrand parameters
    * @param int_param Integration parameters
    * @param res Result
    * @note Rectangle formula, integration by dx
*/
template <typename P, typename KType>
void IntegralUniversalSegPnt(const double* a, const double* b,
        void (*f_0)(const double*, const KernelParam<KType>&, P*),
                              const KernelParam<KType>& ker_param,
                           const IntegralParam& int_param, P* res)
{
    int p_n;
    int n = int_param.GetNStart();
    double d[3], y[3], dl, delta = 0.;
    int idim = int_param.GetIDim();

    P ff[_DIM_MAX]{}, res_prev[_DIM_MAX]{};
    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
    }


    for (p_n = 0; p_n < int_param.GetPMax(); p_n++) {
        for (int j = 0; j < 3; j++) {
            d[j] = (b[j] - a[j]) / n;
        }
        dl = vec_length(d);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                y[j] = a[j] + (i + 0.5) * d[j];
            }
            f_0(y, ker_param, ff);

            for (int g = 0; g < idim; g++) {
                res[g] += ff[g] * static_cast<P>(dl);
            }
        }

        delta = 0.;
        for (int g = 0; g < idim; g++) {
            delta += std::abs(res[g] - res_prev[g]) * std::abs(res[g] - res_prev[g]);
        }

        if (delta <  int_param.GetEpsAccur() * int_param.GetEpsAccur() && p_n != 0) {
            break;
        }

        n = n * 2;
        for (int g = 0; g < idim; g++) {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
        }
    }

    if (p_n == int_param.GetPMax()) {
        for (int g = 0; g < idim; g++) {
            res[g] = res_prev[g];
        }
    }
}
}       // namespace bielcc
#endif  // _INTEGRAL_UNIVERSAL_SEG_PNT_H_
