#ifndef _KERNEL_LIB_H_
#define _KERNEL_LIB_H_

#include <iostream>
#include <complex>

#include "element_geom.h"
#include "kernel_param.h"
#include "constants.h"

//===========================================================
//--------------Library of basic integrands------------------
//===========================================================
/**
    * F_SimplePot_L - simple layer potential (Laplace equation)
    * F_GradSimplePot_L - gradient of the simple layer potential (Laplace)
    * F_DoublePot_L - double layer potential (Laplace)
    * F_VectorPot_L - vector potential (Laplace)
    *
    * F_SimplePot_H - simple layer potential (Helmholtz)
    * F_SimplePot_Hm1 - simple layer potential minus 1 (Helmholtz)
    * F_GradSimplePot_H - gradient of the simple layer potential (Helmholtz)
    * F_DoublePot_H - double layer potential (Helmholtz)
    * F_RotRotSimplePot_H - rot rot of simple layer potential (Helmholtz)
    * 
    * F_xmyDivr2 - special function F = (x - y) / |x - y|^2
    * F_EincDotOrt - special function F = E_inc(x) * ort
    * 
    * Ker_AverDiv - surface div averaging kernel (far from edge)
    * Ker_AverDiv_NearEdge - surface div averaging kernel (near from edge)
    * Ker_BasicAverFunc - common averaging function (base)
    * 
    * 
    * Smoothing function: 3t^2 - 2t^3, t = r / eps
    * 
    * 
    * ===============================================================
    *                ALL FUNCTIONS WITHOUT 1 / 4pi
    * ===============================================================
 */




namespace bielcc {

static inline double sqr(double x) { return x * x; }
static inline std::complex<double> sqr(std::complex<double> x) { return x * x; }

//===========================================================
//-------------Simple layer potential kernel-----------------
//---------------for the Laplace equation--------------------
//===========================================================
/**
    * @brief F(x - y) = 1 / |x - y|
    * @details Kernel of the simple layer potential (Laplace)
    * @param x First arg
    * @param y Second arg
    * @param ker_param Parameters of the integrand
    * @param res Scalar res value
*/
template<typename KType>
void F_SimplePot_L(const double *x, const double *y,
            const KernelParam<KType> &ker_param, double *res)
{
    double r = dist(x, y);
    double t = r / ker_param.smoothR;
    res[0] = (t < 1) * (t * (3. - 2. * t) / ker_param.smoothR) + (t >= 1) / r;
}




//===========================================================
//------Gradient of the simple layer potential kernel--------
//--------------for the Laplace equation---------------------
//===========================================================
/**
    * @brief F = (y - x) / |x - y|^3
    * @details Kernel of the gradient of the simple layer potential (Laplace)
    * @param x First arg
    * @param y Second arg
    * @param ker_param Parameters of the integrand
    * @param res Vector [3] res value
    * @note grad_x ( LayerPotLaplace )
*/
template<typename KType>
void F_GradSimplePot_L(const double *x, const double *y,
                    const KernelParam<KType> &ker_param, double *res)
{
    double r = dist(x, y);
    double t = r / ker_param.smoothR;
    for (int i = 0; i < 3; i++) {
        res[i] = (y[i] - x[i]);
        res[i] *= (t < 1) * (3. - 2. * t) / (r * sqr(ker_param.smoothR)) +
                                                    (t >= 1) / (sqr(r) * r);
    }
}




//===========================================================
//--------------Double layer potential kernel----------------
//----------------for the Laplace equation-------------------
//===========================================================
/**
    * @brief F = n(y) x (x - y) / |x - y|^3
    * @details Kernel of the double layer potential (Laplace)
    * @param x First arg
    * @param y Second arg
    * @param ker_param Parameters of the integrand
    * @param res Scalar res value
    * @note d / dn (LayerPotLaplace). Vector n should be
    * passed through ker_param.n[3]
*/
template<typename KType>
void F_DoublePot_L(const double *x, const double *y,
                const KernelParam<KType> &ker_param, double *res)
{
    double r = dist(x, y);
    double diff[3];
    diff[0] = x[0] - y[0], diff[1] = x[1] - y[1], diff[2] = x[2] - y[2];
    
    res[0] = scal_prod(ker_param.n, diff);
    double t = r / ker_param.smoothR;
    res[0] *= (t < 1) * (3. - 2. * t) / (r * sqr(ker_param.smoothR)) +
                                    (t >= 1) / (sqr(r) * r);
}







//===========================================================
//-----Vector potential kernel for the Laplace equation------
//===========================================================
/**
    * @brief F = a * (y - x) / |x - y|^3
    * @details Vector potential for the Laplace equation
    * @param x First arg
    * @param y Second arg
    * @param ker_param Parameters of the integrand
    * @param res Vector [3] res value
    * @note Vector a should be passed through ker_param.vec_dbl[3]
*/
template<typename KType>
void F_VectorPot_L(const double *x, const double *y,
                const KernelParam<KType>& ker_param, double *res)
{
    double r = dist(x, y);
    if (r < Tolerance_Constants::MACHINE_ZERO) {
        for (int i = 0; i < 3; i++) {
            res[i] = 0.;
        }
    } else {
        double f = 1. / r / r / r; // K(x, y)

        double diff[3];
        diff[0] = (y[0] - x[0]) * f;
        diff[1] = (y[1] - x[1]) * f;
        diff[2] = (y[2] - x[2]) * f;

        vec_prod(diff, ker_param.vec_dbl, res);

        if (r >= Tolerance_Constants::MACHINE_ZERO &&
                           r < ker_param.smoothR)
        {
            double t = r / ker_param.smoothR;
            for (int i = 0; i < 3; i++) {
                res[i] = res[i] * (3 * t * t - 2 * t * t * t);
            }
        }
    }
}



















//===========================================================
//--------------Simple layer potential kernel----------------
//----------------for the Helmholtz equation-----------------
//===========================================================
/**
    * @brief F = e^ik|x - y| / |x - y|
    * @details Kernel of the simple layer potential (Helmholtz)
    * @param x First arg
    * @param y Second arg
    * @param ker_param Parameters of the integrand
    * @param res Scalar res value
*/
template<typename KType>
void F_SimplePot_H(const double *x, const double *y,
                const KernelParam<KType> &ker_param, std::complex<double> *res)
{
    double r = dist(x, y);
    res[0] = std::exp(std::complex<double>(0., r) * ker_param.k);
    double t = r / ker_param.smoothR;
    res[0] *= (t < 1) * (t * (3. - 2. * t) / ker_param.smoothR) + (t >= 1) / r;
}


//===========================================================
//--------------Simple layer potential kernel----------------
//--------------for the Helmholtz equation - 1---------------
//===========================================================
/**
    * @brief F = (e^ik|x - y| - 1) / |x - y|
    * @details Kernel for the simple layer potential
    * (Helmholtz) minus 1
    * @param x First arg
    * @param y Second arg
    * @param ker_param Parameters of the integrand
    * @param res Vector [3] res value
    * @note Used in analytical extraction of a feature in SLPH
*/
template<typename KType>
void F_SimplePot_Hm1(const double *x, const double *y,
                         const KernelParam<KType> &ker_param,
                            std::complex<double> *res)
{
    double r = dist(x, y);

    if (std::abs(std::complex<double>(0., r) * ker_param.k) < 1e-7) {
        res[0] = (std::complex<double>(0., 1.) *
                            ker_param.k - ker_param.k * ker_param.k * r / 2.);
    } else {
        res[0] = std::exp(std::complex<double>(0., r) * ker_param.k) - 1.;
        double t = r / ker_param.smoothR;
        res[0] *= (t < 1) * (t * (3. - 2. * t) /
                                ker_param.smoothR) + (t >= 1) / r;
    }
}


//===========================================================
//-----------Gradient of the simple layer potential----------
//------------kernel for the Helmholtz equation--------------
//===========================================================
/**
    * @brief F =  (ikr - 1) * e^ikr * (x - y) / r^3
    * @details Kernel of the gradient of the simple layer
    * potential (Helmholtz)
    * @param x First arg
    * @param y Second arg
    * @param ker_param Parameters of the integrand
    * @param res Vector [3] res value
    * @note grad_x( LayerPotHelmholtz )
*/
template<typename KType>
void F_GradSimplePot_H(const double *x, const double *y,
                           const KernelParam<KType> &ker_param,
                              std::complex<double> *res)
{
    double r = dist(x, y);
    double diff[3];
    diff[0] = x[0] - y[0], diff[1] = x[1] - y[1], diff[2] = x[2] - y[2];

    double t = r / ker_param.smoothR;
    std::complex<double> ikr = std::complex<double>(0., r) * ker_param.k;

    for (int i = 0; i < 3; i++) {
        res[i] = diff[i] * std::exp(ikr);
        res[i] *= (ikr - 1.);
        res[i] *= (t < 1) * (3. - 2. * t) / (r * sqr(ker_param.smoothR)) +
                                                    (t >= 1) / (sqr(r) * r);
    }
}




//===========================================================
//----------Gradient of the simple layer potential-----------
//-----------kernel for the Helmholtz equation - 1-----------
//===========================================================
/**
    * @brief F =  [1 + (ikr - 1) * e^ikr] * (x - y) / r^3
    * @details Kernel of the gradient of the simple layer
    * potential (Helmholtz) - 1
    * @param x First arg
    * @param y Second arg
    * @param ker_param Parameters of the integrand
    * @param res Vector [3] res value
    * @note grad_x( LayerPotHelmholtz - 1)
*/
template<typename KType>
void F_GradSimplePot_Hm1(const double *x, const double *y,
                      const KernelParam<KType> &ker_param,
                                std::complex<double> *res)
{
    double r = dist(x, y);
    double ikr_abs = std::abs(std::complex<double>(0., r) * ker_param.k);

    if (ikr_abs < 1e-6) {
        for (int i = 0; i < 3; i++) {
            res[i] = 0.;
        }
    } else if (ikr_abs < 1e-4) {
        // предельный переход
        for (int i = 0; i < 3; i++) {
            res[i] = (y[i] - x[i]) * ker_param.k * ker_param.k / (2. * r);
        }
    } else {
        // честно
        std::complex<double> ikr = std::complex<double>(0., r) * ker_param.k;
        for (int i = 0; i < 3; i++) {
            res[i] = (x[i] - y[i]) * (1. + (ikr - 1.) * std::exp(ikr)) / (r * r * r);
        }
    }
}



//===========================================================
//--------------Double layer potential kernel----------------
//---------------for the Helmholtz equation------------------
//===========================================================
/**
    * @brief F = (n(y) * (x - y) / |x - y|^3) * (e^ik|x - y|
    * - ikr e^ik|x - y|)
    * @details Double layer potential kernel for the
    * Helmholtz equation
    * @param x First arg
    * @param y Second arg
    * @param res Scalar res value
    * @param ker_param Parameters of the integrand
    * @note d / dn ( LayerPotHelmholtz ). Vector n should be
    * passed through ker_param.n[3]
*/
template<typename KType>
void F_DoublePot_H(const double *x, const double *y,
                        const KernelParam<KType> &ker_param,
                          std::complex<double> *res)
{
    double r = dist(x, y);
    double diff[3];
    diff[0] = x[0] - y[0], diff[1] = x[1] - y[1], diff[2] = x[2] - y[2];
    std::complex<double> ikr = std::complex<double>(0., r) * ker_param.k;

    res[0] = scal_prod(ker_param.n, diff) * exp(ikr) * (1. - ikr);

    double t = r / ker_param.smoothR;
    res[0] *= (t < 1) * (3 - 2 * t) / (r * sqr(ker_param.smoothR)) +
                                                    (t >= 1) / (sqr(r) * r);
}













//===========================================================
//---------Rot_x Rot_x of the simple layer potential---------
//------------kernel for the Helmholtz equation--------------
//===========================================================
/**
    * @brief F = e^(ikr) * (j_vec * f1 + (x - y)(x - y, j_vec) * f2 / r^2)
    * 
    * f1 = -1 / r^3 + (ik) / r^2 + k^2 / r
    * 
    * f2 = 3 / r^3 - (3ik) / r^2 - k^2 / r
    * @details Rot_x Rot_x of a simple layer potential (Helmholtz)
    * @param x First arg
    * @param y Second arg
    * @param ker_param Parameters of the integrand (currents vectors)
    * @param res Vector [3] res value
    * @note Used as an integrand for K[] operator in far field zone (x out of the cell).
    * Is an analogue of the Hertzian dipole potential kernel function.
    * Vector j_vec should be passed through ker_param.current[3]
*/            
template<typename KType>                                           
void F_RotRotSimplePot_H(const double *x, const double *y,
                             const KernelParam<KType> &ker_param,
                               std::complex<double> *res)
{
    double r = dist(x, y);
    if (r < Tolerance_Constants::MACHINE_ZERO) {
        for (int i = 0; i < 3; i++) {
            res[i] = std::complex<double>(0., 0);
        }
    } else {
        std::complex<double> ik = std::complex<double>(0., 1.) * ker_param.k;
        std::complex<double> sqr_k = sqr(ker_param.k);

        double diff[3];
        for (int i = 0; i < 3; i++) {
            diff[i] = (x[i] - y[i]) / r;
        }
        std::complex<double> f1 = (-1.0 / sqr(r) / r) +
                                        (ik / sqr(r)) +
                                            (sqr_k / r);

        std::complex<double> f2 = (3.0 / sqr(r) / r) -
                                 (3.0 * ik / sqr(r)) -
                                            (sqr_k / r);

        std::complex<double> eikr = std::exp(ik * r);
        f2 *=  scal_prod(diff, ker_param.current);
        for (int i = 0; i < 3; i++) {
            res[i] = (diff[i] * f2 +  ker_param.current[i] * f1) * eikr;
        }

        if (r < ker_param.smoothR) {
            double t = r / ker_param.smoothR;
            for (int i = 0; i < 3; i++) {
                res[i] *= (t * t * (3 - 2 * t));
            }
        }
    }
}




//===========================================================
//---------------F(x, y) =  (x - y) / |x - y|^2--------------
//===========================================================
/**
    * @brief F = (x - y) / |x - y|^2
    * @param x First arg
    * @param y Second arg
    * @param ker_param Parameters of the integrand
    * @param res Vector [3] res value
*/ 
template<typename KType>
void F_xmyDivr2(const double *x, const double *y,
            const KernelParam<KType> &ker_param, double *res)
{
    double diff[3];
    for (int i = 0; i < 3; i++) {
        diff[i] = x[i] - y[i];
    }
    if (vec_length(diff) < Tolerance_Constants::MACHINE_ZERO) {
        for (int i = 0; i < 3; i++) {
            res[i] = 0.;
        }
    } else {
        for (int i = 0; i < 3; i++) {
            res[i] = diff[i] / vec_length(diff) / vec_length(diff);
        }
    }
}




//===========================================================
//---------------F(x, y) = Einc(x) * ort (cell edges)--------
//===========================================================
/**
    * @brief F = E0 * e^{ik_vec x} * ort = E_inc * ort
    * @param x Cal point
    * @param ker_param Parameters of the integrand
    * @param res Scalar res value
    * @note k_vec should be passed through ker_param.vec_cmplx[3], 
    * ort - ker_param.ort[3]
*/ 
template<typename KType>
void F_EincDotOrt(const double *x, const KernelParam<KType> &ker_param,
                                      std::complex<double> *res)
{
    std::complex<double> deg = std::exp(std::complex<double>(0., 1.) *
                        scal_prod(ker_param.vec_cmplx, x));
    res[0] = scal_prod(ker_param.e0, ker_param.ort) * deg;
}












//===========================================================
//--------------Kernel averaging for surface-----------------
//-------------divergence (far from the edge)----------------
//===========================================================
/**
    * @brief Kernel averaging for surface divergence (far from edge)
    * @param x First arg
    * @param y Second arg
    * @param res Vector [3] res value
    * @note eps_edge - criterion of proximity to the edge. "Far"
    * means dist(x, edge) > 3 * eps_edge
*/
template<typename KType>
void Ker_AverDiv(const double *x, const double *y,
              const KernelParam<KType> &ker_param, double *res)
{
    double eps_edge = ker_param.eps_edge;
    double reps = dist(x, y) / eps_edge;
    double deg = (-6.0 + 2. * reps * reps) * exp(-reps * reps)
            / M_PI / sqr(eps_edge) / sqr(eps_edge); 
    for (int i = 0; i < 3; i++) {
        res[i] = (x[i] - y[i]) * deg;
    }
}



//===========================================================
//--------------Kernel averaging for surface-----------------
//--------------divergence (near to the edge)----------------
//===========================================================
/**
    * @brief Kernel averaging for shallow divergence (near the edge)
    * @param x First arg
    * @param y Second arg
    * @param res Vector [3] res value
    * @note eps_edge - criterion of proximity to the edge. "Near"
        * means dist(x, edge) < 3 * eps_edge
*/                        
template<typename KType>                                   
void Ker_AverDiv_NearEdge(const double *x, const double *y,
                       const KernelParam<KType> &ker_param, double *res)
{
    double eps_edge = ker_param.eps_edge;
    double reps = dist(x, y) / eps_edge; // eps from formula, not smppthing
    double deg = (-20.0 + 8. * reps * reps) * exp(-reps * reps) / M_PI /
                 eps_edge / eps_edge / eps_edge / eps_edge;
    for (int i = 0; i < 3; i++) {
        res[i] = (x[i] - y[i]) * deg;
    }
}




//===========================================================
//--------------Basic averaging function---------------------
//===========================================================
/**
    * @brief Basic averaging function
    * @details psi(r = |x - y|) = e^(r^2 / eps^2) / (pi * eps^2)
    * @param x First var
    * @param y Second var
    * @param res Vector [3] res value
*/                                                             
inline double Ker_BasicAverFunc(const double *x, const double *y,
                                         const double eps)
{
    double res, deg;

    if (dist(x, y) > eps) {
        return 0.;
    }
    deg = dist(x, y) * dist(x, y) / eps / eps;
    return exp(deg) / M_PI / eps / eps;
}

}       // namespace bielcc
#endif  // _KERNEL_LIB_H_
