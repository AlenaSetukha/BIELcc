#ifndef _INTEGRALS_ANALYTIC_H_
#define _INTEGRALS_ANALYTIC_H_

#include <iostream>
#include <array>
#include <cstdlib>
#include <complex>



#include "kernel_param.h"
#include "integral_param.h"
#include "integral_universal.h"
#include "kernel_lib.h"
#include "element_geom.h"
#include "tolerance_constants.h"

namespace bielcc {

//===========================================================
//-----------------Integral 1 / |x - y|----------------------
//-------------Written by Ilya Mass @vnebelunka--------------
//===========================================================
/**
    * @brief Analytical surface integral F(x, y) = 1 / |x - y|
    * @param rut0 Integration cell [4][3]
    * @param x Collocation point for F(x, y) (fixed)
    * @result Integral value
    * @note Integration by dy, the point x lies inside the cell
*/
double integral1Divr(const double (&rut0)[4][3], const double* x);
/**
    * @brief Analytical surface integral F(x, y) = 1 / |x - y|
    * @param rut0 Integration cell [3][3]
    * @param x Collocation point for F(x, y) (fixed)
    * @result Integral value
    * @note Integration by dy, the point x lies inside the cell
*/
double integral1Divr(const double (&rut0)[3][3], const double *x);


/**
    * @brief Analytical curvilinear integral grad_x F(x, y) = grad_x (1 / |x - y|)
    * @param rut0 Integration cell [3][3]
    * @param x Collocation point for F(x, y) (fixed)
    * @result Integral value
    * @note Integration by dy, the point x lies inside the cell
*/
void integralGrad1DivrSeg(const double* A, const double* B,
                            const double *x, double* res);







//===========================================================
//---------------Integral d / dn(1 / |x - y|)----------------
//===========================================================
/**
    * @brief Analytical surface integral
    * F(xk, y) = d / dn(1 / |x - y|)
    * @param rut0 Integration cell [4][3]
    * @param x Collocation point for F(x, y) (fixed)
    * @result Integral value
    * @note Integration by dy. This is the oriented solid
    * angle. x does not lie on the cell (otherwise = 0)
*/
double integral_ddn_1Divr(const double (&rut0)[4][3],
                                    const double* x);






//===========================================================
//-----------Integral (x - y) / (|x - y| ^ 2)----------------
//===========================================================
/**
    * @brief Analytical surface integral
    * F(xk, y) = (x - y) / (|x - y| ^ 2)
    * @param rut0 Integration cell [4][3]
    * @param x Collocation point for F(x, y) (fixed)
    * @param res Integral vector value result
    * @note Integration by dy, x must lie on the cell
*/
void integral_xmyDivr2(const double (&rut0)[4][3],
                       const double* x, double* res);



                       


//===========================================================
//-----------Integral nu / |x - y| (curvilinear)-------------
//===========================================================
/**
    * @brief Analytical curvilinear integral 
    * F(xk, y) = nu / |x - y|
    * @param rut0 Integration cell [4][3]
    * @param x Collocation point for F(x, y) (fixed)
    * @param eps Smoothing radius (proximity of a point to the edge of a cell)
    * @param res Integral vector value result
    * @note Integration by dy, nu - outward normals to each side of the cell
*/
void integralnu1Divr(const double (&rut0)[4][3], const double* x,
                                   const double eps, double* res);








//===========================================================
//------------Semi-analytical surface integral---------------
//----------------over a cell of a function------------------
//--------------F(x-y) = e^ik|x - y| / |x - y|---------------
//===========================================================
/**
    * @brief Semi-analytical surface integral of a SimplePot_H kernel
    * F(x, y) = e^{ikr} / r using hilghlighting singularity.
    * @param cell integration cell
    * @param x collocation point for F(x, y) (fixed)
    * @param analyticR_r radius for highlighting singularity (relative to cell diam)
    * @param int_param numerical parameters for integration
    * @param k_wave wave number
    * @result Complex value
    * @note Integration by dy, highlight singularity in the integrand.
    * No dividing on 4pi.
*/
template<size_t CellPoints, typename KType>
std::complex<double> IntegralSurf_SimplePot_H_HS(const double (&cell)[CellPoints][3],
                                                 const double *x, double analyticR_r,
                                        const IntegralParam& int_param, KType k_wave)
{
    double y[3];
    get_center_mass(cell, y);
    std::complex<double> tmp[1];
    std::complex<double> ff = std::complex<double>(0., 0.);

    KernelParam<KType> param;
    param.k = k_wave;
    param.smoothR = Tolerance_Constants::MACHINE_ZERO;

    double calc_dist = analyticR_r * get_diam(cell);
    double r = dist(x, y);
    if (r < calc_dist) {
        // close, calculate analytically, dividing into 2 integrals
        ff = static_cast<std::complex<double>>(integral1Divr(cell, x));
        IntegralUniversalPnt(x, cell, bielcc::F_SimplePot_Hm1, param, int_param, tmp);
        ff += tmp[0];
    } else {
        // far, calculate the integral numerically. No smoothing
        IntegralUniversalPnt(x, cell, bielcc::F_SimplePot_H, param, int_param, tmp);
        ff = tmp[0];
    }
    return ff;
}









//===========================================================
//-----------Semi-analytical curvilinear integral------------
//---------------over a segment of a function----------------
//---------grad_x F(x-y) = grad_x (e^ik|x - y| / |x - y|)----
//----------------using higlighting singularity-------------
//===========================================================
/**
    * @brief Semi-analytical curv integral of a GradSimplePot_H kernel
    * gradx F(x, y) = grad_x (e^{ikr} / r) using hilghlighting singularity.
    * @param A point A
    * @param B point B
    * @param x collocation point for gradx F(x, y) (fixed)
    * @param analyticR_r radius for highlighting singularity (relative to cell diam)
    * @param int_param numerical parameters for integration
    * @param k_wave wave number
    * @result Complex vector [3]
    * @note Integration by dy, highlight singularity in the integrand.
    * No dividing on 4pi.
*/
template<typename KType>
void IntegralSeg_GradSimplePot_H_HS(const double* A, const double* B,
                                  const double *x, double analyticR_r,
                         const IntegralParam& int_param, KType k_wave,
                                            std::complex<double>* res)
{
    KernelParam<KType> param;
    param.k = k_wave;
    param.smoothR = Tolerance_Constants::MACHINE_ZERO;

    double y[3] = {(B[0] + A[0]) / 2.,
                   (B[1] + A[1]) / 2.,
                   (B[2] + A[2]) / 2.};
    double calc_dist = analyticR_r * dist(A, B);
    double r = dist(x, y);
    if (r < calc_dist) {
        double tmp[3]{};
        // close, calculate analytically, dividing into 2 integrals
        integralGrad1DivrSeg(A, B, x, tmp);
        IntegralUniversalSegPnt(A, B, x, bielcc::F_GradSimplePot_Hm1, param, int_param, res);
        
        for (int i = 0; i < 3; i++) {
            res[i] += std::complex<double>(tmp[i]);
        }
    } else {
        // far, calculate the integral numerically. No smoothing
        IntegralUniversalSegPnt(A, B, x, bielcc::F_GradSimplePot_H, param, int_param, res);
    }
}

}       // namespace bielcc
#endif  // _INTEGRALS_ANALYTIC_H_