#ifndef _DOUBLE_SURF_INTEGRAL_SPH_H_
#define _DOUBLE_SURF_INTEGRAL_SPH_H_

#include <iostream>
#include <cmath>
#include <complex>

#include "kernel_param.h"
#include "integral_param.h"
#include "integral_universal.h"
#include "integrals_analytic.h"


namespace bielcc {

constexpr static const int _PMAX_SMOOTH_DOUBLE_SURF = 8;

/**
    * @brief Calculation of the integral over a cell of the semi-analytical integral
    * over a cell from F(x - y) = e^{ikr} / r. Outer integral - numerically,
    * inner - semi-analytically.
    * @param cellOut outer integration cell
    * @param cellIn inner integration cell
    * @param param integrand parameters
    * @param analyticR_r relativeto h analytical radius ( = 3)
    * @param IP_cell surface integration parameters
    * @note Rectangle formula. Wave number "k_wave" = param.k
*/
template<size_t CellPoints, typename KType>
std::complex<double>
DoubleSurfIntegral_SimplePot_H_HS(const double (&cellOut)[4][3],
                          const double (&cellIn)[CellPoints][3],
                  KernelParam<KType>& param, double analyticR_r,
                                   const IntegralParam& IP_cell)
{
    double p, q, p1, q1, s;
    double a[3], b[3], a1[3], a2[3], a3[3], a4[3], m1[3], m2[3], rc[3], rn[3];
    double start_rs = ker_param.smoothR;

    int n = IP_cell.GetNStart();
    const int PMax = IP_cell.GetPMax();
    const double epsOut = IP_cell.GetEpsAccur();
    KType k_wave = param.k;

    std::complex<double> ff = std::complex<double>(0., 0.);
    std::complex<double> res_prev = std::complex<double>(0., 0.);
    std::complex<double> res = std::complex<double>(0., 0.);

    int p_n = 0;
    for (p_n = 0; p_n < PMax; p_n++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                p = (double)i / (double)n;
                q = (double)j / (double)n;
                p1 = ((double)i + 1.) / (double)n;
                q1 = ((double)j + 1.) / (double)n;

                for (int k = 0; k < 3; k++) {
                    a[k] = q * cellOut[1][k] + (1. - q) * cellOut[0][k];
                    b[k] = q * cellOut[2][k] + (1. - q) * cellOut[3][k];
                    a1[k] = p * b[k] + (1 - p) * a[k];
                    a4[k] = p1 * b[k] + (1 - p1) * a[k];
                }

                for (int k = 0; k < 3; k++) {
                    a[k] = q1 * cellOut[1][k] + (1. - q1) * cellOut[0][k];
                    b[k] = q1 * cellOut[2][k] + (1. - q1) * cellOut[3][k];
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
                ff = IntegralSurf_SimplePot_H_HS(cellIn, rc, analyticR_r, IP_cell, k_wave);
                res += ff * s;
            }
        }

        double delta = std::abs(res - res_prev) * std::abs(res - res_prev);

        if (delta < epsOut * epsOut && p_n != 0) {
            break;
        }

        n = n * 2; // the step is reduced by 2 times -> smoothing radius divided by 2
        res_prev = res;
        res = std::complex<double>(0., 0.);
        if (p_n < _PMAX_SMOOTH_DOUBLE_SURF && ker_param.smoothR > ker_param.smoothRMin) {
            ker_param.smoothR /= 2.;
        }
    }

    if (p_n == PMax) {
        res = res_prev;
    }

    ker_param.smoothR = start_rs;
    return res;
}








template<size_t CellPoints, typename KType>
std::complex<double>
DoubleSurfIntegral_SimplePot_H_HS(const double (&cellOut)[3][3],
                          const double (&cellIn)[CellPoints][3],
                  KernelParam<KType>& ker_param, double analyticR_r,
                                   const IntegralParam& IP_cell)
{
    double s, rc[3], A_jP[3], B_j[3], C_jP[3], A_jM[3], C_jM[3];
    double p_vec[3], q_vec[3], r_vec[3];
    double delta = 0.;
    int n = IP_cell.GetNStart();
    int idim = IP_cell.GetIDim();
    double eps_acc = IP_cell.GetEpsAccur();
    const int PMax = IP_cell.GetPMax();
    KType k_wave = ker_param.k;

    double start_rs = ker_param.smoothR;

    std::complex<double> ff = std::complex<double>(0., 0.);
    std::complex<double> res_prev = std::complex<double>(0., 0.);
    std::complex<double> res = std::complex<double>(0., 0.);

    s = tr_square(cellOut[0], cellOut[1], cellOut[2]);

    int p_n = 0;
    for (p_n = 0; p_n < PMax; p_n++) {
        for (int k = 0; k < 3; k++) {
            p_vec[k] = (cellOut[1][k] - cellOut[0][k]) / n;
            q_vec[k] = (cellOut[2][k] - cellOut[0][k]) / n;
            r_vec[k] = (cellOut[2][k] - cellOut[1][k]) / n;
        }
        for (int i = 0; i < n; i++) { // m
            for (int j = 0; j < n - i; j++) { // k
                for (int k = 0; k < 3; k++) {
                    B_j[k] = cellOut[0][k] + (i + 1) * p_vec[k] + j * q_vec[k];
                    A_jM[k] = B_j[k] - p_vec[k];
                    C_jM[k] = A_jM[k] + q_vec[k];
                    rc[k] = (B_j[k] + A_jM[k] + C_jM[k]) / 3.0;
                }


                // Расчет внутреннего интеграла
                ff = IntegralSurf_SimplePot_H_HS(cellIn, rc, analyticR_r, IP_cell, k_wave);

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
                    ff = IntegralSurf_SimplePot_H_HS(cellIn, rc, analyticR_r, IP_cell, k_wave);
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
        res_prev = res;
        res = std::complex<double>(0., 0.);
        if (p_n < _PMAX_SMOOTH_DOUBLE_SURF && ker_param.smoothR > ker_param.smoothRMin) {
            ker_param.smoothR /= 2.;
        }
    }

    if (p_n == PMax) {
        res = res_prev;
    }

    ker_param.smoothR = start_rs;
    return res;
}



}           // namespace bielcc
#endif      // _DOUBLE_SURF_INTEGRAL_SPH_H_
