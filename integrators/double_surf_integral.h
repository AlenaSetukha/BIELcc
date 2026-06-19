#ifndef _DOUBLE_SURF_INTEGRAL_H_
#define _DOUBLE_SURF_INTEGRAL_H_

#include <iostream>
#include <cmath>

#include "kernel_param.h"
#include "integral_param.h"
#include "integral_universal.h"


namespace bielcc {

constexpr static const int _DIM_MAX_CELL = 3;
constexpr static const int _PMAX_SMOOTH_DOUBLE_SURF = 8;


/**
    * @brief Double surface integral over a quadrangular cell from
    * the F(x, y). Numerically (dy) dx.
    * @param cellOut outer integratoion cell
    * @param cellIn inner integration cell
    * @param f_0 integrand (F(x_pnt - y))
    * @param param integrand parameters
    * @param IP_cell surface integration parameters
    * @param res result
    * @note Rectangle formula
*/
template <size_t CellPoints, typename KType, typename P>
void DoubleSurfIntegral(const double (&cellOut)[4][3], const double (&cellIn)[CellPoints][3],   
        void (*f_0)(const double*, const double*, const KernelParam<KType>&, P*),
        const KernelParam<KType>& ker_param_const, const IntegralParam& IP_cell, P* res)
{
    KernelParam<KType> ker_param = ker_param_const;
    double p, q, p1, q1, s, a[3], b[3], a1[3], a2[3], a3[3], a4[3], m1[3], m2[3], rc[3], rn[3];

    int n = IP_cell.GetNStart();
    int idim = IP_cell.GetIDim();
    const double epsOut = IP_cell.GetEpsAccur();
    const int PMax = IP_cell.GetPMax();

    double start_rs = ker_param.smoothR;


    P ff[_DIM_MAX_CELL]{}, res_prev[_DIM_MAX_CELL]{};
    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
    }

    int p_n = 0;
    for (p_n = 0; p_n < PMax; p_n++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                p = (double)i / (double)n;
                q = (double)j / (double)n;
                p1 = ((double)i + 1.) / (double)n;
                q1 = ((double)j + 1.) / (double)n;

                for (int k = 0; k < 3; k++) {
                    a[k] = q * cellOut[1][k] + (1 - q) * cellOut[0][k];
                    b[k] = q * cellOut[2][k] + (1 - q) * cellOut[3][k];
                    a1[k] = p * b[k] + (1 - p) * a[k];
                    a4[k] = p1 * b[k] + (1 - p1) * a[k];
                }

                for (int k = 0; k < 3; k++) {
                    a[k] = q1 * cellOut[1][k] + (1 - q1) * cellOut[0][k];
                    b[k] = q1 * cellOut[2][k] + (1 - q1) * cellOut[3][k];
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

                // Расчет внутреннего интеграла
                IntegralUniversalPnt(rc, cellIn, f_0, ker_param, IP_cell, ff);

                for (int g = 0; g < idim; g++) {
                    res[g] += ff[g] * static_cast<P>(s);
                }
            }
        }

        double delta = 0.;
        for (int g = 0; g < idim; g++) {
            delta += std::abs(res[g] - res_prev[g]) * std::abs(res[g] - res_prev[g]);
        }

        if (delta < epsOut * epsOut && p_n != 0) {
            break;
        }

        n = n * 2; // the step is reduced by 2 times -> smoothing radius divided by 2
        for (int g = 0; g < idim; g++) {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
            if (p_n < _PMAX_SMOOTH_DOUBLE_SURF && ker_param.smoothR > ker_param.smoothRMin) {
                ker_param.smoothR /= 2.;
            }
        }

    }

    if (p_n == PMax) {
        for (int g = 0; g < idim; g++) {
            res[g] = res_prev[g];
        }
    }
}







/**
    * @brief Double surface integral over a triangular cell from
    * the F(x, y). Numerically (dy) dx.
    * @param cellOut outer integratoion cell
    * @param cellIn inner integration cell
    * @param f_0 integrand (F(x_pnt - y))
    * @param param integrand parameters
    * @param IP_cell surface integration parameters
    * @param res result
    * @note Rectangle formula
*/
template <size_t CellPoints, typename KType, typename P>
void DoubleSurfIntegral(const double (&cellOut)[3][3], const double (&cellIn)[CellPoints][3],   
        void (*f_0)(const double*, const double*, const KernelParam<KType>&, P*),
        const KernelParam<KType>& ker_param_const, const IntegralParam& int_param, P* res)
{
    KernelParam<KType> ker_param = ker_param_const;
    double s, rc[3], A_jP[3], B_j[3], C_jP[3], A_jM[3], C_jM[3];
    double p_vec[3], q_vec[3], r_vec[3];
    double delta = 0.;
    int n = int_param.GetNStart();
    int idim = int_param.GetIDim();
    double eps_acc = int_param.GetEpsAccur();
    const int PMax = int_param.GetPMax();

    double start_rs = ker_param.smoothR;

    P ff[_DIM_MAX_CELL]{}, res_prev[_DIM_MAX_CELL]{};
    for (int g = 0; g < idim; g++) {
        res[g] = static_cast<P>(0);
    }

    s = tr_square(cell[0], cell[1], cell[2]);

    int p_n = 0;
    for (p_n = 0; p_n < PMax; p_n++) {
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


                // Расчет внутреннего интеграла
                IntegralUniversalPnt(rc, cellIn, f_0, ker_param, int_param, ff);

                for (int g = 0; g < idim; g++) {
                    res[g] += ff[g];
                }

                if (j < n - 1 - i) {
                    for (int k = 0; k < 3; k++) {
                        A_jP[k] = B_j[k] + r_vec[k];
                        C_jP[k] = B_j[k] + q_vec[k];
                        rc[k] = (B_j[k] + A_jP[k] + C_jP[k]) / 3.0;
                    }

                    // Расчет внутреннего интеграла
                    IntegralUniversalPnt(rc, cellIn, f_0, param, int_param, ff);
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

        if (delta <  eps_acc * eps_acc && p_n != 0) {
            break;
        }

        n = n * 2; // the step is reduced by 2 times -> smoothing radius divided by 2
        for (int g = 0; g < idim; g++) {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
            if (p_n < _PMAX_SMOOTH_DOUBLE_SURF && ker_param.smoothR > ker_param.smoothRMin) {
                ker_param.smoothR /= 2.;
            }
        }
    }

    if (p_n == PMax) {
        for (int g = 0; g < idim; g++) {
            res[g] = res_prev[g];
        }
    }
}
}           // namespace bielcc     
#endif      // _DOUBLE_SURF_INTEGRAL_H_