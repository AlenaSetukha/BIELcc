#ifndef _INCIDENT_FIELD_H_
#define _INCIDENT_FIELD_H_

#include <iostream>
#include <complex>
#include "material.h"
#include "kernel_param.h"

//=======================================================
//-----------Incident field functions--------------------
//=======================================================
namespace bielcc {
    /**
        * @brief Plane wave electric field function
        * @details Einc(x) = e0 * e^i(k, x)
        * @param x calculation point
        * @param param parameters for plane wave (k_vec as vec_dbl, e0_ort)
        * @param e_inc vector [3] res value
    */ 
    template<typename KType>
    void PlaneWaveFieldE(const double* x, const KernelParam<KType>& param,
                                              std::complex<double>* e_inc)
    {
        std::complex<double> deg =
            std::exp(std::complex<double>(0., 1.) * scal_prod(param.vec_dbl, x));
        e_inc[0] = param.e0[0] * deg;
        e_inc[1] = param.e0[1] * deg;
        e_inc[2] = param.e0[2] * deg;
    }


    /**
        * @brief Plane wave magnetic field function
        * @details Hinc(x) = e^i(k, x) * (k x eo) / (omega * mu * mu0)
        * @param x calculation point
        * @param param HInc function parameters (e0_ort, k_vec as vec_dbl, freq)
        * @param material_param volume material parameters (eps, mu)
        * @param h_inc result
    */ 
    template<typename KType>
    void PlaneWaveFieldH(const double* x, const KernelParam<KType>& param,
                                      const MaterialParam& material_param,
                                              std::complex<double>* h_inc)
    {
        double vp[3];
        vec_prod(param.vec_dbl, param.e0, vp);
        std::complex<double> deg = std::exp(std::complex<double>(0., 1.) *
                                                scal_prod(param.vec_dbl, x));
        deg /= (2. * M_PI * param.freq);

        try {
            deg /= std::get<std::complex<double>>(material_param.GetMu());
        } catch (const std::bad_variant_access& e) {
            deg /= std::get<double>(material_param.GetMu());
        }
        h_inc[0] = vp[0] * deg;
        h_inc[1] = vp[1] * deg;
        h_inc[2] = vp[2] * deg;
    }

}       // namespace bielcc
#endif  //_INCIDENT_FIELD_H_