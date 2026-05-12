#ifndef _R_H_
#define _R_H_

#include <iostream>
#include <complex>
#include <variant>

#include "num_param.h"
#include "kernel_param.h"
#include "integral_param.h"
#include "kernel_lib.h"
#include "element_geom.h"
#include "tolerance_constants.h"
#include "integral_universal.h"
#include "K.h"


namespace bielcc {

//===========================================================
//------------Operator R = Rot [surf integral]---------------
//---------in the near zone using smoothing kernel-----------
//===========================================================
/**
    * @brief Near-field integral operator R[cell, j_vec](x) (smoothing)
    * @details Fully numerical calculation of the near-field integral
    * operator
    *        R[cell, j_vec](x) = rot surf_int(j_vec, F) =
    *               = -surf_int( j_vec x grad_x F)
    * where F = eikr/r.
    * Integral is calculated with smoothing.
    * 
    * Calculation at fixed point x. Rectangle formula, 
    * integration over dy.
    * @param j argument vector [3]
    * @param x calc point
    * @param cell integration cell [4 / 3][3]
    * @param num_param numerical parameters for integration
    * @param k wave number for integrand
    * @param res result
    * @note Operator R[] calculated with smoothing kernel.
    * 
    * num_param.rs - the integral smoothing radius
    * relative to the second-level cell size (~5 - 10 h2)
*/
template<typename JType, typename KType, size_t CellPoints>
void R_Rot_Smooth(const JType* j, const double* x, const double (&cell)[CellPoints][3],
                    const NumParam& num_param, const KType k, std::complex<double>* res)
{
    std::complex<double> cur_res3[3]{};

    KernelParam<KType> param;
    param.smoothR = num_param.rs * get_diam(cell) / num_param.n_start;
    param.smoothRMin = get_diam(cell) / 100.;
    param.k = k;
    IntegralParam int_parGradF(3, num_param.n_start, num_param.p_max, num_param.eps);

    IntegralUniversalPnt(x, cell, F_GradSimplePot_H, param, int_parGradF, cur_res3);
    vec_prod(cur_res3, j, res);
    res[0] /= (4. * M_PI), res[1] /= (4. * M_PI), res[2] /= (4. * M_PI);
}



//===========================================================
//------------Operator R = Rot [surf integral]---------------
//----------------in the collocation points------------------
//===========================================================
/**
    * @brief Integral operator R[cell, j_vec](x) only in colloc points. No smooth.
    * 
    * @details Numerical calculation of the integral operator
    * 
    *        R[cell, j_vec](x) = rot surf_int(j_vec, F) =
    *               = -surf_int( j_vec x grad_x F)
    * 
    * where F = eikr/r.
    * If x is cell center: R[] = 0, else - surf integral with no smoothing.
    * 
    * Calculation at fixed point x. Rectangle formula, integration over dy.
    * @param j argument vector [3]
    * @param x collocation point
    * @param cell Integration cell [4 / 3][3]
    * @param num_param Numerical parameters for integration
    * @param k Wave number for integrand
    * @param res Result
    * @note Operator R[] calculated in collocation points. No smoothing.
*/
template<typename JType, typename KType, size_t CellPoints>
void R_Rot_Colloc(const JType* j, const double* x, const double (&cell)[CellPoints][3],
                   const NumParam& num_param, const KType k,  std::complex<double>* res)
{
    double y[3];
    get_center_mass(cell, y);
    if (dist(x, y) < Tolerance_Constants::MACHINE_ZERO) { 
        // i == j
        res[0] = 0., res[1] = 0., res[2] = 0.;
    } else {
        std::complex<double> cur_res3[3]{};
        KernelParam<KType> param;
        param.smoothR = Calculation_Constants::MACHINE_ZERO;
        param.smoothRMin = Calculation_Constants::MACHINE_ZERO;
        param.k = k;
        
        IntegralParam int_parGradF(3, num_param.n_start, num_param.p_max, num_param.eps);

        IntegralUniversalPnt(x, cell, F_GradSimplePot_H, param, int_parGradF, cur_res3);
        vec_prod(cur_res3, j, res);
        res[0] /= (4. * M_PI), res[1] /= (4. * M_PI), res[2] /= (4. * M_PI);
    }
}





//===========================================================
//------------Operator R = Rot [surf integral]---------------
//---------in the near zone using smoothing kernel-----------
//===========================================================
/**
    * @brief Near-field integral operator R[cell, j_vec](x) as a curl and surf integral sum.
    * @details Fully numerical calculation of the near-field integral
    * operator
    *        R[cell, j_vec](x) = j_vec x curl_int[F(x-y) * nu] dy +
    * j_vec x surf_int[n * dF(x-y) / dn] dy
    * where F = eikr/r.
    * Surf integral is calculated with smoothing, curl - with no smoothing.
    * 
    * Calculation at fixed point x. Rectangle formula, 
    * integration over dy.
    * @param j argument vector [3]
    * @param x calc point
    * @param cell integration cell [4 / 3][3]
    * @param num_param numerical parameters for integration
    * @param k wave number for integrand
    * @param res result
    * @note Operator R[] calculated with surf and curl integral sum (smooth and no smooth).
*/
template<typename JType, typename KType, size_t CellPoints>
void R_Rot_Near(const JType* j, const double* x, const double (&cell)[CellPoints][3],
            const NumParam& num_param, const KType k, std::complex<double>* res)
{
    // surf_int (n * F_DoublePot_H)
    double norm[3];
    norm_func(cell, norm);
    std::complex<double> surf_int[3]{}, res_surf[3]{};
    KernelParam<KType> param;
    param.smoothR = num_param.rs * get_diam(cell) / num_param.n_start;
    param.smoothRMin = get_diam(cell) / 100.;
    param.k = k;
    param.n[0] = norm[0], param.n[1] = norm[1], param.n[2] = norm[2];
    IntegralParam int_parFDpot(3, num_param.n_start, num_param.p_max, num_param.eps);
    IntegralUniversalPnt(x, cell, nF_DoublePot_H, param, int_parFDpot, res_surf);


    // -curl_int (nu * F_DoublePot_H)
    std::complex<double> res_curl[3]{};
    KernelParam<KType> param_seg;
    param_seg.smoothR = Calculation_Constants::MACHINE_ZERO;
    param_seg.smoothRMin = Calculation_Constants::MACHINE_ZERO;
    param_seg.k = k;
    
    IntegralParam int_parnuF(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);
    K0_Common(x, cell, nuF_SimplePot_H, param_seg, int_parnuF, res_curl);


    // result
    std::complex<double> res_tmp[3]{};
    for (int i = 0; i < 3; i++) {
        res_tmp[i] = (res_surf[i] + res_curl[i]) / (4. * M_PI);
    }

    vec_prod(j, res_tmp, res);
}




}       // namespace bielcc
#endif  // _R_H_