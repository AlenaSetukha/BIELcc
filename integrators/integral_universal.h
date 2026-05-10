#ifndef _INTEGRAL_UNIVERSAL_H_
#define _INTEGRAL_UNIVERSAL_H_

#include <iostream>

#include "kernel_param.h"
#include "integral_param.h"
#include "element_geom.h"



constexpr static const int _DIM_MAX_CELL = 3;

namespace bielcc {
//===========================================================
//-------Integral in point over a cell by F(xk, y)-----------
//===========================================================
/**
    * Calculates the surface integral over a quadrangular
    * (triangular) cell of a function  F(xk, y) at fixed
    * point xk. Rectangle formula, integration over dy.
*/


/**
    * @brief Surface integral over a quadrangular cell from the F(xk, y)
    * @param x Collocation point for F(x, y) (fixed)
    * @param cell Integration cell [4][3]
    * @param f_0 Integrand (F(x_pnt - y))
    * @param ker_param Integrand parameters
    * @param int_param Integration parameters
    * @param res Result
    * @note Rectangle formula, integration by dy
*/
template<typename P, typename KType>
void IntegralUniversalPnt(const double* x, const double (&cell)[4][3],
        void (*f_0)(const double*, const double*, const KernelParam<KType>&, P*),
        KernelParam<KType>& ker_param, const IntegralParam& int_param, P* res)
{
    double p, q, p1, q1, s, a[3], b[3], a1[3], a2[3], a3[3], a4[3], m1[3], m2[3], rc[3], rn[3];
    double delta = 0.;
    int n = int_param.GetNStart();
    int idim = int_param.GetIDim();
    int p_n = 0;


    P ff[_DIM_MAX_CELL]{}, res_prev[_DIM_MAX_CELL]{};

    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
    }


    for (p_n = 0; p_n < int_param.GetPMax(); p_n++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                p = (double)i / (double)n;
                q = (double)j / (double)n;
                p1 = ((double)i + 1.) / (double)n;
                q1 = ((double)j + 1.) / (double)n;

                for (int k = 0; k < 3; k++) {
                    a[k] = q * cell[1][k] + (1 - q) * cell[0][k];
                    b[k] = q * cell[2][k] + (1 - q) * cell[3][k];
                    a1[k] = p * b[k] + (1 - p) * a[k];
                    a4[k] = p1 * b[k] + (1 - p1) * a[k];
                }

                for (int k = 0; k < 3; k++) {
                    a[k] = q1 * cell[1][k] + (1 - q1) * cell[0][k];
                    b[k] = q1 * cell[2][k] + (1 - q1) * cell[3][k];
                    a2[k] = p * b[k] + (1 - p) * a[k];
                    a3[k] = p1 * b[k] + (1 - p1) * a[k];
                }

                for (int k = 0; k < 3; k++) {
                    rc[k] = (a1[k] + a2[k] + a3[k] + a4[k]) / 4.0;
                    m1[k] = ((a2[k] + a3[k]) - (a1[k] + a4[k])) / 2.0;
                    m2[k] = ((a3[k] + a4[k]) - (a1[k] + a2[k])) / 2.0;
                }
                vec_prod(m1, m2, rn);
                s = vec_length(rn);

                f_0(x, rc, ker_param, ff);

                for (int g = 0; g < idim; g++) {
                    res[g] += ff[g] * static_cast<P>(s);
                }
            }
        }

        delta = 0.;
        for (int g = 0; g < idim; g++) {
            delta += std::abs(res[g] - res_prev[g]) * std::abs(res[g] - res_prev[g]);
        }

        if (delta <  int_param.GetEpsAccur() * int_param.GetEpsAccur() && p_n != 0) {
            break;
        }

        n = n * 2; // the step is reduced by 2 times -> smoothing radius divided by 2
        for (int g = 0; g < idim; g++) {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
            ker_param.smoothR /= 2.;
        }
    }
    
    if (p_n == int_param.GetPMax()) {
        for (int g = 0; g < idim; g++) {
            res[g] = res_prev[g];
        }
    }
}





/**
    * @brief Surface integral over a triangular cell from the F(xk, y).
    * @param x Collocation point for F(x, y) (fixed)
    * @param cell Integration cell [3][3]
    * @param f_0 Integrand
    * @param ker_param Integrand parameters
    * @param int_param Integration parameters
    * @param res Result
    * @note Rectangle formula, integration by dy.
*/
template<typename P, typename KType>
void IntegralUniversalPnt(const double* x, const double (&cell)[3][3],
        void (*f_0)(const double*, const double*, const KernelParam<KType>&, P*),
        KernelParam<KType>& ker_param, const IntegralParam& int_param, 
        P* res)
{
    double s, rc[3], A_jP[3], B_j[3], C_jP[3], A_jM[3], C_jM[3];
    double p_vec[3], q_vec[3], r_vec[3];
    double delta = 0.;
    int n = int_param.GetNStart();
    int idim = int_param.GetIDim();
    int p_n = 0;



    P ff[_DIM_MAX_CELL]{}, res_prev[_DIM_MAX_CELL]{};
    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
    }

    s = tr_square(cell[0], cell[1], cell[2]);


    for (p_n = 0; p_n < int_param.GetPMax(); p_n++) {
        for (int k = 0; k < 3; k++) {
            p_vec[k] = (cell[1][k] - cell[0][k]) / n;
            q_vec[k] = (cell[2][k] - cell[0][k]) / n;
            r_vec[k] = (cell[2][k] - cell[1][k]) / n;
        }
        for (int i = 0; i < n; i++) { // m
            for (int j = 0; j < n - i; j++) { // k
                for (int k = 0; k < 3; k++) {
                    B_j[k] = cell[0][k] + (i + 1) * p_vec[k] + j * q_vec[k];
                    A_jM[k] = B_j[k] - p_vec[k];
                    C_jM[k] = A_jM[k] + q_vec[k];
                    rc[k] = (B_j[k] + A_jM[k] + C_jM[k]) / 3.0;
                }

                f_0(x, rc, ker_param, ff);
                for (int g = 0; g < idim; g++) {
                    res[g] += ff[g];
                }

                if (j < n - 1 - i) {
                    for (int k = 0; k < 3; k++) {
                        A_jP[k] = B_j[k] + r_vec[k];
                        C_jP[k] = B_j[k] + q_vec[k];
                        rc[k] = (B_j[k] + A_jP[k] + C_jP[k]) / 3.0;
                    }

                    f_0(x, rc, ker_param, ff);
                    for (int g = 0; g < idim; g++) {
                        res[g] += ff[g];
                    }
                }
            }
        }

	    delta = 0.;
        for (int g = 0; g < idim; g++) {
            res[g] = res[g] * s / (double)n / (double)n;
            delta += std::abs(res[g] - res_prev[g]) * std::abs(res[g] - res_prev[g]);
        }

        if (delta <  int_param.GetEpsAccur() * int_param.GetEpsAccur() && p_n != 0) {
            break;
        }

        n = n * 2; // the step is reduced by 2 times -> smoothing radius divided by 2
        for (int g = 0; g < idim; g++) {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
            ker_param.smoothR /= 2.;
        }
    }

    if (p_n == int_param.GetPMax()) {
        for (int g = 0; g < idim; g++) {
            res[g] = res_prev[g];
        }
    }
}












//===========================================================
//-------------Integral over a cell by F(x)------------------
//===========================================================
/**
    * Calculates the surface integral over a quadrangular
    * (triangular) cell of a function  F(x). Rectangle formula,
    * integration over dx.
*/


/**
    * @brief Surface integral over a quadrangular cell from the F(x)
    * @param cell Integration cell [4][3]
    * @param f_0 Integrand F(x)
    * @param ker_param Integrand parameters
    * @param int_param Integration parameters
    * @param res Result
    * @note Rectangle formula, integration by dx
*/
template<typename P, typename KType>
void IntegralUniversal(const double (&cell)[4][3],
            void (*f_0)(const double*, const KernelParam<KType>&, P*),
                                        KernelParam<KType>& ker_param,
                               const IntegralParam& int_param, P* res)
{
    double p, q, p1, q1, s;
    double a[3], b[3], a1[3], a2[3], a3[3], a4[3], m1[3], m2[3], rc[3], rn[3];
    double delta = 0.;
    int n = int_param.GetNStart();
    int idim = int_param.GetIDim();
    int p_n = 0;

    P ff[_DIM_MAX_CELL]{}, res_prev[_DIM_MAX_CELL]{};
    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
    }


    for (p_n = 0; p_n < int_param.GetPMax(); p_n++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                p = (double)i / (double)n;
                q = (double)j / (double)n;
                p1 = ((double)i + 1.) / (double)n;
                q1 = ((double)j + 1.) / (double)n;

                for (int k = 0; k < 3; k++) {
                    a[k] = q * cell[1][k] + (1. - q) * cell[0][k];
                    b[k] = q * cell[2][k] + (1. - q) * cell[3][k];
                    a1[k] = p * b[k] + (1. - p) * a[k];
                    a4[k] = p1 * b[k] + (1. - p1) * a[k];
                }

                for (int k = 0; k < 3; k++) {
                    a[k] = q1 * cell[1][k] + (1. - q1) * cell[0][k];
                    b[k] = q1 * cell[2][k] + (1. - q1) * cell[3][k];
                    a2[k] = p * b[k] + (1. - p) * a[k];
                    a3[k] = p1 * b[k] + (1. - p1) * a[k];
                }

                for (int k = 0; k < 3; k++) {
                    rc[k] = (a1[k] + a2[k] + a3[k] + a4[k]) / 4.0;
                    m1[k] = ((a2[k] + a3[k]) - (a1[k] + a4[k])) / 2.0;
                    m2[k] = ((a3[k] + a4[k]) - (a1[k] + a2[k])) / 2.0;
                }
                vec_prod(m1, m2, rn);
                s = vec_length(rn);

                f_0(rc, ker_param, ff);

                for (int g = 0; g < idim; g++) {
                    res[g] += ff[g] * static_cast<P>(s);
                }
            }
        }

        delta = 0.;
        for (int g = 0; g < idim; g++) {
            delta += std::abs(res[g] - res_prev[g]) * std::abs(res[g] - res_prev[g]);
        }

        if (delta <  int_param.GetEpsAccur() * int_param.GetEpsAccur() && p_n != 0) {
            break;
        }

        n = n * 2; // the step is reduced by 2 times -> smoothing radius divided by 2
        for (int g = 0; g < idim; g++) {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
            ker_param.smoothR /= 2.;
        }
    }

    if (p_n == int_param.GetPMax()) {
        for (int g = 0; g < idim; g++) {
            res[g] = res_prev[g];
        }
    }
}



/**
    * @brief Surface integral over a triangular cell from the F(x)
    * @param cell Integration cell [3][3]
    * @param f_0 Integrand F(x)
    * @param ker_param Integrand parameters
    * @param int_param Integration parameters
    * @param res Result
    * @note Rectangle formula, integration by dx
*/
template<typename P, typename KType>
void IntegralUniversal(const double (&cell)[3][3],
            void (*f_0)(const double*, const KernelParam<KType>&, P*),
                                        KernelParam<KType>& ker_param,
                               const IntegralParam& int_param, P* res)
{
    double s, rc[3], A_jP[3], B_j[3], C_jP[3], A_jM[3], C_jM[3];
    double p_vec[3], q_vec[3], r_vec[3];
    double delta = 0.;
    int n = int_param.GetNStart();
    int idim = int_param.GetIDim();
    int p_n = 0;



    P ff[_DIM_MAX_CELL]{}, res_prev[_DIM_MAX_CELL]{};
    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
    }

    s = tr_square(cell[0], cell[1], cell[2]);


    for (p_n = 0; p_n < int_param.GetPMax(); p_n++) {
        for (int k = 0; k < 3; k++) {
            p_vec[k] = (cell[1][k] - cell[0][k]) / n;
            q_vec[k] = (cell[2][k] - cell[0][k]) / n;
            r_vec[k] = (cell[2][k] - cell[1][k]) / n;
        }
        for (int i = 0; i < n; i++) { // m
            for (int j = 0; j <= n - 1 - i; j++) { // k
                for (int k = 0; k < 3; k++) {
                    B_j[k] = cell[0][k] + (i + 1) * p_vec[k] + j * q_vec[k];
                    A_jM[k] = B_j[k] - p_vec[k];
                    C_jM[k] = A_jM[k] + q_vec[k];
                    rc[k] = (B_j[k] + A_jM[k] + C_jM[k]) / 3.0;
                }

                f_0(rc, ker_param, ff);
                for (int g = 0; g < idim; g++) {
                    res[g] += ff[g];
                }

                if (j < n - 1 - i) {
                    for (int k = 0; k < 3; k++) {
                        A_jP[k] = B_j[k] + r_vec[k];
                        C_jP[k] = B_j[k] + q_vec[k];
                        rc[k] = (B_j[k] + A_jP[k] + C_jP[k]) / 3.0;
                    }

                    f_0(rc, ker_param, ff);
                    for (int g = 0; g < idim; g++) {
                        res[g] += ff[g];
                    }
                }
            }
        }

	    delta = 0.;
        for (int g = 0; g < idim; g++) {
            res[g] = res[g] * s / (double)n / (double)n;
            delta += std::abs(res[g] - res_prev[g]) * std::abs(res[g] - res_prev[g]);
        }

        if (delta <  int_param.GetEpsAccur() * int_param.GetEpsAccur() && p_n != 0) {
            break;
        }

        n = n * 2; // the step is reduced by 2 times -> smoothing radius divided by 2
        for (int g = 0; g < idim; g++) {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
            ker_param.smoothR /= 2.;
        }
    }

    if (p_n == int_param.GetPMax()) {
        for (int g = 0; g < idim; g++) {
            res[g] = res_prev[g];
        }
    }
}




}       // namespace bielcc
#endif  // _INTEGRAL_UNIVERSAL_H_