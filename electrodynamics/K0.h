#ifndef _K0_H_
#define _K0_H_

#include <iostream>
#include <cstdlib>

#include "kernel_param.h"
#include "integral_param.h"
#include "element_geom.h"
#include "tolerance_constants.h"
#include "integral_universal_seg.h"

constexpr static const int _DIM_MAX_SEG = 3;

namespace bielcc {
//===========================================================
//------Contour integral K0 over a cell by F(xk, y)----------
//===========================================================
/**
    * @brief Contour integral K0[sigma, j_vec] applied to F(x, y) (sign "-")
    * @details Calculates the contour integral by integral operator
    * K0[sigma, j_vec] over the sides of a cell from some
    * function F(x_pnt, y):
    *              K0 = -int_curv { (j, nu) * F(x - y) dy}.
    * Calculation at fixed point x. Rectangle formula, 
    * integration over dy.
    * @param j argument vector [3] (constant on cell)
    * @param x calc point
    * @param cell Integration cell [4 / 3][3]
    * @param f_0 Integrand (F(x_pnt - y))
    * @param ker_param Integrand parameters
    * @param int_param Integration parameters
    * @param res Result
    * @note Rectangle formula, integration by dy
*/
template <typename P, typename KType, size_t CellPoints>
void K0(const P* j, const double* x, const double (&cell)[CellPoints][3],
        void (*f_0)(const double*, const double*, const KernelParam<KType>&, P*),
           KernelParam<KType>& ker_param, const IntegralParam& int_param, P* res)
{
    int idim = int_param.GetIDim();
    P ks[_DIM_MAX_SEG]{};
    double len_nu;
    double nu[3], diff[3];

    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
        ks[g] = static_cast<P>(0);
    }


    double norm[3];
    norm_func(cell, norm);

    for (int s = 0; s < (int)CellPoints; s++) {
        for (int g = 0; g < idim; g++) {
            ks[g] = static_cast<P>(0);
            diff[g] = 0.;
        }
        int Apnt, Bpnt;
        if (s != (int)CellPoints - 1) {
            Apnt = s, Bpnt = s + 1;
        } else {
            Apnt = (int)CellPoints - 1, Bpnt = 0;
        }
        IntegralUniversalSegPnt(cell[Apnt], cell[Bpnt], x, f_0,
                                        ker_param, int_param, ks);
        for (int g = 0; g < idim; g++) {
            diff[g] = cell[Bpnt][g] - cell[Apnt][g]; // b - a
        }

        vec_prod(diff, norm, nu);
        len_nu = vec_length(nu);

        if (len_nu > Tolerance_Constants::MACHINE_ZERO) {
            for (int g = 0; g < idim; g++) {
                nu[g] /= len_nu;
            }
        }

        for (int g = 0; g < idim; g++) {
            ks[g] *= -scal_prod(j, nu); // ks[3] ab
            res[g] += ks[g];
        }
    }
}







//===========================================================
//---Contour integral K0 over a cell by F(xk, y) * nu--------
//===========================================================
/**
    * @brief Contour integral K0[sigma, j_vec] = curl_int [F(x, y) * nu(y)] dy
    * @details Calculates the contour integral over the sides of a cell from some
    * function F(x_pnt, y):
    *              K0_Common = int_curv { F(x-y) * nu } dy
    * Calculation at fixed point x. Rectangle formula, 
    * integration over dy.
    * @param x calc point
    * @param cell integration cell [4 / 3][3]
    * @param f_0 integrand (F(x_pnt - y))
    * @param ker_param integrand parameters
    * @param int_param Integration parameters
    * @param res result
    * @note Rectangle formula, integration by dy
*/
template <typename P, typename KType, size_t CellPoints>
void K0_Common(const double* x, const double (&cell)[CellPoints][3],
        void (*f_0)(const double*, const double*, const KernelParam<KType>&, P*),
           KernelParam<KType>& ker_param, const IntegralParam& int_param, P* res)
{
    int idim = int_param.GetIDim();
    P ks[_DIM_MAX_SEG]{};
    double len_nu;
    double nu[3], diff[3];

    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
        ks[g] = static_cast<P>(0);
    }


    double norm[3];
    norm_func(cell, norm);

    for (int s = 0; s < (int)CellPoints; s++) {
        for (int g = 0; g < idim; g++) {
            ks[g] = static_cast<P>(0);
            diff[g] = 0.;
        }
        int Apnt, Bpnt;
        if (s != (int)CellPoints - 1) {
            Apnt = s, Bpnt = s + 1;
        } else {
            Apnt = (int)CellPoints - 1, Bpnt = 0;
        }


        for (int g = 0; g < idim; g++) {
            diff[g] = cell[Bpnt][g] - cell[Apnt][g]; // b - a
        }
        vec_prod(diff, norm, nu);
        len_nu = vec_length(nu);
        if (len_nu > Tolerance_Constants::MACHINE_ZERO) {
            for (int g = 0; g < idim; g++) {
                nu[g] /= len_nu;
            }
        }

        ker_param.vec_dbl[0] = nu[0];
        ker_param.vec_dbl[1] = nu[1];
        ker_param.vec_dbl[2] = nu[2];
        IntegralUniversalSegPnt(cell[Apnt], cell[Bpnt], x, f_0,
                                        ker_param, int_param, ks);
        

        for (int g = 0; g < idim; g++) {
            res[g] += ks[g];
        }
    }
}






//===========================================================
//-----------Semi-analytical contour integral K0-------------
//------------over a cell edge of a function-----------------
//----------gradx F(x-y) = gradx (e^ik|x - y| / |x - y|)-----
//----------------using higlighting singularity--------------
//===========================================================
/**
    * @brief Semi-analytical curvilinear integral of a GradSimplePot_H kernel
    * gradx F(x, y) = gradx (e^{ikr} / r) using hilghlighting singularity (sign "-")
    * @details Calculates the contour integral operator K0[sigma, j_vec]
    * over the cell edges from grad_x F(x_pnt, y):
    *              K0 = -int_curv { (j, nu) * grad_x F(x - y) dy,
    * where grad_x F(x-y) = grad_x (e^ik|x - y| / |x - y|)
    * using higlighting singularity.
    * Calculation at fixed point x. Rectangle formula, 
    * integration over dy.
    * @param j argument vector [3] (constant on cell)
    * @param x calc point
    * @param cell Integration cell [4 / 3][3]
    * @param ker_param Integrand parameters
    * @param int_param Integration parameters
    * @param res Result
    * @note Rectangle formula, integration by dy
*/
template <typename P, typename KType, size_t CellPoints>
void K0_GradSimplePot_H_HS(const P* j, const double* x,
                    const double (&cell)[CellPoints][3],
                         const IntegralParam& int_param,
                       double analyticR_r, KType k_wave,
                              std::complex<double>* res)
{
    const int idim = 3;
    std::complex<double> ks[3]{};
    double len_nu, nu[3]{}, diff[3]{};
    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
    }


    double norm[3];
    norm_func(cell, norm);

    for (int s = 0; s < (int)CellPoints; s++) {
        for (int g = 0; g < idim; g++) {
            ks[g] = static_cast<P>(0.);
            diff[g] = 0.;
        }

        int Apnt, Bpnt;
        if (s != (int)CellPoints - 1) {
            Apnt = s, Bpnt = s + 1;
        } else {
            Apnt = (int)CellPoints - 1, Bpnt = 0;
        }

        IntegralSeg_GradSimplePot_H_HS(cell[Apnt], cell[Bpnt],
                        x, analyticR_r, int_param, k_wave, ks);

        for (int g = 0; g < idim; g++) {
            diff[g] = cell[Bpnt][g] - cell[Apnt][g]; // b - a
        }

        vec_prod(diff, norm, nu);
        len_nu = vec_length(nu);

        if (len_nu > Tolerance_Constants::MACHINE_ZERO) {
            for (int g = 0; g < idim; g++) {
                nu[g] /= len_nu;
            }
        }

        for (int g = 0; g < idim; g++) {
            ks[g] *= -scal_prod(j, nu); // ks[3] ab
            res[g] += ks[g];
        }
    }
}

}       // namespace bielcc
#endif  // _K0_H_