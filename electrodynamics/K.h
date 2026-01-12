#ifndef _K_H_
#define _K_H_

#include <iostream>
#include <complex>
#include <array>
#include <variant>

#include "element_geom.h"
#include "num_param.h"
#include "kernel_param.h"
#include "integral_param.h"
#include "integral_universal.h"
#include "kernel_lib.h"
#include "K0.h"
#include "integrals_analytic.h"



namespace bielcc {
//===========================================================
//-----------Operator K = Rot Rot[surf integral]-------------
//---------in the near zone using smoothing kernel-----------
//===========================================================
/**
    * Fully numerical calculation of the near-field integral
    * operator
    *        K[cell, j_vec](x) = -curv_int(gradx F) + k^2 * surf_int(F),
    * where F = eikr/r.
    * Both integrals are calculated with smoothing.
    * 
    * Calculation at fixed point x. Rectangle formula, 
    * integration over dy.
    * Note: num_param.rs - the integral smoothing radius
    * relative to the second-level cell size (~0.5 - 2)
*/

/**
    * @brief Near-field integral operator K[cell, j_vec] (just smoothing)
    * @param j argument vector [3]
    * @param x calc point
    * @param cell Integration cell [4 / 3][3]
    * @param num_param Numerical parameters for integration
    * @param k Wave number for integrand
    * @param res Result
    * @note Operator K[] calculated by smoothing kernels
*/
template<typename JType, typename KType, size_t CellPoints>
void K_RotRot_Near_Smooth(const JType* j, const double* x, const double (&cell)[CellPoints][3],
                           const NumParam& num_param, const KType k, std::complex<double>* res)
{
    KernelParam<KType> param_F;
    param_F.smoothR = num_param.rs * get_diam(cell) / num_param.n_start;
    param_F.k = k;

    KernelParam<KType> param_GradF;
    param_GradF.smoothR = num_param.rs_seg * get_diam(cell) / num_param.n_start_seg;
    param_GradF.k = k;

    IntegralParam int_parF(1, num_param.n_start, num_param.p_max, num_param.eps);
    IntegralParam int_parGradF(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);

    // k^2
    std::complex<double> cur_res[1];
    IntegralUniversalPnt(x, cell, F_SimplePot_H, param_F, int_parF, cur_res);
    std::complex<double> surf_int = cur_res[0] * k * k / (4. * M_PI);


    // grad div -> -curl integral
    std::complex<double> cur_res3[3];
    K0(j, x, cell, F_GradSimplePot_H, param_GradF, int_parGradF, cur_res3);
    

    // -curl_int + k^2 surf_int j
    for (int i = 0; i < 3; i++) {
        res[i] = surf_int * j[i] + cur_res3[i] / (4. * M_PI);
    }
}    











//===========================================================
//-----------Operator K = Rot Rot[surf integral]-------------
//------------------in the near zone using-------------------
//-----------------highlitghting singularity-----------------
//===========================================================
/**
    * Semi-analytical calculation of the near-field integral
    * operator
    *        K[cell, j_vec](x) = -curv_int(gradx F) + k^2 * surf_int(F),
    * where F = eikr / (4pi r).
    * 
    * curv_int(gradx F) - calculated numerically.
    * surf_int (F) - highlitghting singularity + Semi-analytical
    * 
    * Calculation at fixed point x. Rectangle formula, 
    * integration over dy.
    * Note: num_param.rs - the integral smoothing radius
    * relative to the second-level cell size (~0.5 - 2)
*/

/**
    * @brief Near-field semi-analytical integral operator K[cell, j_vec] (highlight singularity).
    * Smoothing in curv integral
    * @param j argument vector [3]
    * @param x calc point
    * @param cell Integration cell [4 / 3][3]
    * @param num_param Numerical parameters for integration
    * @param k Wave number for integrand
    * @param res Result
    * @note Operator K[] calculated by highlight singularity. Smooth in curv integral.
*/
template<typename JType, typename KType, size_t CellPoints>
void K_RotRot_HS(const JType* j, const double* x, const double (&cell)[CellPoints][3],
                  const NumParam& num_param, const KType k, std::complex<double>* res)
{
    // k^2 with highlight singularity
    IntegralParam int_parF(1, num_param.n_start, num_param.p_max, num_param.eps);
    std::complex<double> cur_res = IntegralSurf_SimplePot_H_HS(cell, x,
                                num_param.analytic_calc, int_parF, k);

    for (int i = 0; i < 3; i++) {
        res[i] = cur_res * k * k * j[i] / (4. * M_PI);
    }


    // -curv
    std::complex<double> res3[3];
    KernelParam<KType> param_seg;
    param_seg.smoothR = num_param.rs_seg * get_diam(cell) / num_param.n_start_seg;
    param_seg.k = k;
    
    IntegralParam  int_parGradF(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);
    K0(j, x, cell, F_GradSimplePot_H, param_seg, int_parGradF, res3);
    for (int i = 0; i < 3; i++) {
        res[i] += res3[i] / (4. * M_PI);
    }
}





//===========================================================
//-----------Operator K = Rot Rot[surf integral]-------------
//------------------in the near zone using-------------------
//-----------------highlitghting singularity-----------------
//===========================================================
/**
    * Semi-analytical calculation of the near-field integral
    * operator
    *        K[cell, j_vec](x) = -curv_int(gradx F) + k^2 * surf_int(F),
    * where F = eikr / (4pi r).
    * 
    * curv_int(gradx F) - calculated numerically.
    * surf_int (F) - highlitghting singularity + Semi-analytical
    * 
    * Calculation at fixed point x. Rectangle formula, 
    * integration over dy.
    * Note: num_param.rs - the integral smoothing radius
    * relative to the second-level cell size (~0.5 - 2)
*/

/**
    * @brief Near-field semi-analytical integral operator K[cell, j_vec] in 
    * collocation points. No smoothing in curv integral.
    * @param j argument vector [3]
    * @param x collocation point
    * @param cell Integration cell [4 / 3][3]
    * @param num_param Numerical parameters for integration
    * @param k Wave number for integrand
    * @param res Result
    * @note Operator K[] calculated by highlight singularity. No smooth
*/
template<typename JType, typename KType, size_t CellPoints>
void K_RotRot_HS_Colloc(const JType* j, const double* x, const double (&cell)[CellPoints][3],
                  const NumParam& num_param, const KType k, std::complex<double>* res)
{
    // k^2 with highlight singularity
    IntegralParam int_parF(1, num_param.n_start, num_param.p_max, num_param.eps);
    std::complex<double> cur_res = IntegralSurf_SimplePot_H_HS(cell, x,
                                num_param.analytic_calc, int_parF, k);
    for (int i = 0; i < 3; i++) {
        res[i] = cur_res * k * k * j[i] / (4. * M_PI);
    }


    // -curv
    std::complex<double> res3[3]{};
    KernelParam<KType> param_seg;
    param_seg.smoothR = 10e-16;
    param_seg.k = k;
    
    IntegralParam int_parGradF(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);
    K0(j, x, cell, F_GradSimplePot_H, param_seg, int_parGradF, res3);
    for (int i = 0; i < 3; i++) {
        res[i] += res3[i] / (4. * M_PI);
    }
}




/**
    * @brief Near-field semi-analytical integral operator K[cell, j_vec] in 
    * collocation points. Highlighting singularity in both terms.
    * @param j argument vector [3]
    * @param x collocation point
    * @param cell Integration cell [4 / 3][3]
    * @param num_param Numerical parameters for integration
    * @param k Wave number for integrand
    * @param res Result
    * @note Operator K[] calculated by highlight singularity. No smooth
*/
template<typename JType, typename KType, size_t CellPoints>
void K_RotRot_HS2_Colloc(const JType* j, const double* x, const double (&cell)[CellPoints][3],
                  const NumParam& num_param, const KType k, std::complex<double>* res)
{
    // k^2 with highlight singularity
    IntegralParam int_parF(1, num_param.n_start, num_param.p_max, num_param.eps);
    std::complex<double> cur_res = IntegralSurf_SimplePot_H_HS(cell, x,
                                num_param.analytic_calc, int_parF, k);
    for (int i = 0; i < 3; i++) {
        res[i] = cur_res * k * k * j[i] / (4. * M_PI);
    }


    // -curv with highlight singularity
    std::complex<double> res3[3]{};
    
    IntegralParam int_parGradF(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);
    K0_GradSimplePot_H_HS(j, x, cell, int_parGradF, num_param.analytic_calc, k, res3);

    for (int i = 0; i < 3; i++) {
        res[i] += res3[i] / (4. * M_PI);
    }
}









//===========================================================
//-----------Operator K = Rot Rot[surf integral]-------------
//--------------------in the far zone------------------------
//===========================================================
/**
    * Fully numerical calculation of the far-field integral
    * operator
    *           K[cell, j_vec](x) = Rot Rot (surf_int)
    *              K[cell, j_vec](x) = surf_int(Kfar)
    * with a special integrand
    *               Kfar = e^(ikr) * (j_vec * f1 +
    *            (x - y)(x - y, j_vec) * f2 / r^2) / 4pi ,
    * where
    *           f1 = -1 / r^3 + (ik) / r^2 + k^2 / r,
    *           f2 = 3 / r^3 - (3ik) / r^2 - k^2 / r.
    * Calculation at fixed point x. Rectangle formula, 
    * integration over dy.
*/



/**
    * @brief Far-field integral operator K[cell, j_vec]. No smoothing.
    * @param j argument vector [3]
    * @param x calc point
    * @param cell Integration cell [4 / 3][3]
    * @param num_param Numerical parameters for integration
    * @param k Wave number for integrand
    * @param res Result
    * @note Far field - x is not on the cell.
*/
template<typename JType, typename KType, size_t CellPoints>
void K_RotRot_Far(const JType* j, const double* x, const double (&cell)[CellPoints][3],
                   const NumParam& num_param, const KType k, std::complex<double>* res)
{
    KernelParam<KType> param_KFar;
    param_KFar.smoothR = 10e-16;
    param_KFar.k = k;
    param_KFar.current[0] = j[0], param_KFar.current[1] = j[1], param_KFar.current[2] = j[2];
    IntegralParam int_parF(3, num_param.n_start, num_param.p_max, num_param.eps);

    IntegralUniversalPnt(x, cell, F_RotRotSimplePot_H, param_KFar, int_parF, res);
    res[0] /= (4. * M_PI), res[1] /= (4. * M_PI), res[2] /= (4. * M_PI);
}


}       // namespace bielcc 
#endif  // _K_H_