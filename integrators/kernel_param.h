#ifndef _KERNEL_PARAM_H_
#define _KERNEL_PARAM_H_

#include <iostream>
#include <complex>
#include <optional>
#include <variant>
#include <array>

namespace bielcc {
/** 
    * @brief Numerical parameters of the integration kernel function
    * @param smoothR absolute smoothing radius of the integrand
    * @param smoothRMin minimum absolute smoothing radius of the integrand
    * @param k wave number of the source function 
    * 
    * @param vec_cmplx / vec_dbl - additional vectors (optional)
    * @param n normal vector (optional)
    * @param e0 direction vector of the incident wave (optional)
    * @param ort complementary vector (optional)
    * @param current surface current vector (optional) (complex)
    * @param freq frequency
    * @note smoothR changes during adaptive integration
*/
template<typename KType>
struct KernelParam {
    double smoothR;
    double smoothRMin;
    KType k;

    std::complex<double> vec_cmplx[3]{};
    double vec_dbl[3]{};
    double n[3]{};
    double ort[3]{};
    double e0[3]{};
    std::complex<double> current[3];
    double freq;

    KernelParam() = default;
    /**
        * @brief Integrand parameters
        * @param k_ wave number
        * @param smoothR_ smoothing radius (absolute value)
    */
    KernelParam(KType k_, double smoothR_):
                    smoothR(smoothR), k(k_) {smoothRMin = 0.;};
    /**
        * @brief Integrand parameters
        * @param k_ wave number
        * @param smoothR_ smoothing radius (absolute value)
        * @param smoothRMin_ smoothing radius (absolute value)
    */
    KernelParam(KType k_, double smoothR_, double smoothRMin_):
                    smoothR(smoothR), k(k_), smoothRMin(smoothRMin_) {};
    KernelParam(const KernelParam &obj);
    ~KernelParam() = default;
};


template<typename KType>
KernelParam<KType>::KernelParam(const KernelParam& obj)
{
    smoothR = obj.smoothR;
    k = obj.k;
    for (int i = 0; i < 3; i++) {
        n[i] = obj.n[i];
        vec_cmplx[i] = obj.vec_cmplx[i];
        vec_dbl[i] = obj.vec_dbl[i];
        ort[i] = obj.ort[i];
        e0[i] = obj.e0[i];
        current[i] = obj.current[i];
    }
    freq = obj.freq;
}
}       // namespace bielcc
#endif  // _KERNEL_PARAM_H_