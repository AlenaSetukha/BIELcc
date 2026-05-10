#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>


#include "integral_universal.h"
#include "integral_universal_seg.h"
#include "kernel_param.h"
#include "integral_param.h"
#include "element_geom.h"
#include "kernel_lib.h"
#include "integrals_analytic.h"
#include "constants.h"

/**
    * Test to check the calculation of the main integrators for
    * quadrilateral and triangular cells.
    * 
    * Test0: integral of a function f1 = 1 over the surface
    * of a cell. Result - cell area. Function f1 without a singularity
    * - no smoothing needed (smooth_r = 10^{-16}).
    * Result: 5 for rut1
    * 
    * 
    * 
    * Test1: surface integral over a cell of a function
    *       F(x - y) = e^{ikr} / (4 pi r), r = |x - y| with k = 1
    * calculated in three ways:
    *       - numerically by a quadrangular cell rut1 (x lies on rut1) with smoothing.
    *       - semi-analytically with the identification of features (no smothing)
    *       - numerically as the sum of two integrals over two triangular cells:
    *         rut2 + rut3, x lies on rut3. Smoothing only on rut3.
    * What parameters must to be set:
    *       - KernelParam<double> kernel_param0; - see the kernel_param.h
    *       - IntegralParam IntParam(); - see the integral_param.h
    * 
    * 
    * Test2: Surface integral of a vector function f0. No add parameters (no smoothing, no k...)
    * 
    * Test3: Segment integral of a vector function f0. No add parameters (no smoothing, no k...)
*/

using namespace bielcc;

void f0(const double *x, const double *y,
                const KernelParam<double> &param, double *ff) 
{
    ff[0] = x[0] * x[1] + 15.9;
    ff[1] = x[2] + x[1];
    ff[2] = x[0];
}

void f1(const double *x, const KernelParam<double> &param,
                                               double *ff) 
{
    ff[0] = 1.;
}




int main(int argc, char **argv)
{
    //==============================================
    //-------Cells for integrators testing----------
    //==============================================
    const double rut1[4][3] = {{0., 0., 0.},
                                {2., 0., 0.},
                                {3., 2., 0.},
                                {0., 2., 0.}};

    const double rut2[3][3] = {{0., 2., 0.},
                               {0., 0., 0.},
                               {2., 0., 0.}};
    
    
    const double rut3[3][3] = {{3., 2., 0.},
                               {2., 0., 0.},
                               {0., 2., 0.}};


    const double x[3] = {2., 1.5, 0.};




    //==============================================
    //-------Global numerical parameters------------
    //==============================================
    const double IntegralAccuracy = Calculation_Constants::INTEGRAL_ACCURACY;
    const double SmoothingR_r = Calculation_Constants::SMOOTHING_DIST_SURF_INT;         // relative to small cell step h2
    const double SmoothingRSeg_r = Calculation_Constants::SMOOTHING_DIST_SEG_INT;       // relative to small cell step h2
    const double AnalyticCalcR = Calculation_Constants::ANALYTIC_CALC_DIST;             // relative to grid step h
    const int NCellStart = Calculation_Constants::START_CELL_SPLIT;                       
    const int NSegStart = Calculation_Constants::START_SEG_SPLIT;
    const int PMaxCell = Calculation_Constants::PMAX_CELL_SPLIT;                        // 2^{PMax} steps in surface integration
    const int PMaxSeg = Calculation_Constants::PMAX_SEG_SPLIT;                          // 2^{PMaxSeg} steps in segment integration
    
    std::cout << "Numerical values ​​of calculation parameters" << std::endl;
    std::cout << "Integrals calculation accuracy: " << IntegralAccuracy << std::endl;
    std::cout << "Relative smoothing radius (surface, rel. to h2): " << SmoothingR_r << std::endl;
    std::cout << "Relative smoothing radius (segment, rel. to h2): " << SmoothingRSeg_r << std::endl;
    std::cout << "Radius of analytical calculation (rel. to h): " << AnalyticCalcR << std::endl;
    std::cout << "Starting cell split: " << NCellStart << std::endl;
    std::cout << "Starting segment split: " << NSegStart << std::endl;
    std::cout << "Limit cell split (2^{P}): " << PMaxCell << std::endl;
    std::cout << "Limit segment split (2^{P}): " << PMaxSeg << std::endl;
    std::cout << std::endl;
 





    //==============================================
    //-------------------Test0----------------------
    //---------------Area of trapezoid--------------
    //-----1D surface integral with integrand f1----
    //==============================================
    KernelParam<double> kernel_param0;                                   // Integrand parameters (not used for f1 = 1)
    IntegralParam IntParam_0(1, NCellStart, PMaxCell, IntegralAccuracy); // Integral calculation parameters

    double s_trap[1]{};
    IntegralUniversal(rut1, f1, kernel_param0, IntParam_0, s_trap);
    std::cout << "Test0. Area of ​​a trapezoid rut1: ";
    std::cout << s_trap[0] << std::endl;
    std::cout << std::endl;







    //==============================================
    //-------------------Test1---------------------+
    //---------1D surface integral in pnt x---------
    //----------with integrand F_SimplePot_H-------- 
    //==============================================
    std::cout << "Test1. Surface integral of a F_SimplePot_H" << std::endl;
    const double k_wave = 1.;
    IntegralParam Intparam_11(1, NCellStart, PMaxCell * 10, IntegralAccuracy); // Numerical integration parameters


    // 1. Numerically over quad cell
    KernelParam<double> kernel_param11;                                        // Integrand parameters
    kernel_param11.smoothR = SmoothingR_r * (get_diam(rut1) / NCellStart);     // Absolute start value of smoothing radius relative to h2
    kernel_param11.k = k_wave;


    std::complex<double> res11[1];
    IntegralUniversalPnt(x, rut1, F_SimplePot_H, kernel_param11, Intparam_11, res11);
    std::cout << "       Fully numerical ([4][3]): " << res11[0] / 4. / M_PI  << std::endl;




    // 2. Semi-analytically over quad cell
    std::complex<double> res12_analytic =
                 IntegralSurf_SimplePot_H_HS(rut1, x, AnalyticCalcR, Intparam_11, k_wave);
    std::cout << "       Analytical calc ([4][3]) = " << res12_analytic / 4. / M_PI  << std::endl;



    
    // 3. Numerically over two triangle cells
    kernel_param11.smoothR = Calculation_Constants::MACHINE_ZERO; // no smoothing on rut2
    kernel_param11.k = k_wave;
    std::complex<double> res13_rut2[1];
    IntegralUniversalPnt(x, rut2, F_SimplePot_H, kernel_param11, Intparam_11, res13_rut2);


    kernel_param11.smoothR = SmoothingR_r * get_diam(rut3) / NCellStart; // absolute value of smoothing radius
    kernel_param11.k = k_wave;
    std::complex<double> res13_rut3[1];
    IntegralUniversalPnt(x, rut3, F_SimplePot_H, kernel_param11, Intparam_11, res13_rut3);
    std::cout << "       Numerically integrals [3][3] + [3][3]: ";
    std::cout << (res13_rut2[0] + res13_rut3[0]) / 4. / M_PI  << std::endl;

    





    //==============================================
    //-------------------Test2---------------------+
    //---------3D surface integral in pnt x---------
    //---------------with integrand f0-------------- 
    //==============================================
    double res2[3];
    KernelParam<double> kernel_param2;                                   // Integrand parameters (not used for f0)
    IntegralParam IntParam_2(3, NCellStart, PMaxCell, IntegralAccuracy); // Integral calculation parameters

    IntegralUniversalPnt(x, rut1, f0, kernel_param2, IntParam_2, res2);
    std::cout << std::endl;
    std::cout << "Test2. Surface integral of a vector function f0" << std::endl;
    std::cout << "       Integral [4][3]: ";
    std::cout << res2[0] << " " << res2[1] << " " << res2[2] << std::endl;
    



    double res3d_2[3], res3d_3[3];
    IntegralUniversalPnt(x, rut2, f0, kernel_param2, IntParam_2, res3d_2);
    IntegralUniversalPnt(x, rut3, f0, kernel_param2, IntParam_2, res3d_3);
    std::cout << "       Integral [3][3] + [3][3]: " << res3d_2[0] + res3d_3[0] << " " << 
                    res3d_2[1] + res3d_3[1] << " " << res3d_2[2] + res3d_3[2] << std::endl;






    //==============================================
    //-------------------Test3---------------------+
    //---------3D segment integral in pnt x---------
    //---------------with integrand f0--------------
    //==============================================
    KernelParam<double> kernel_param3;                                   // Integrand parameters (not used for f0)
    IntegralParam IntParam_3(3, NSegStart, PMaxSeg, IntegralAccuracy);   // Integral calculation parameters


    double res_seg_3[3];
    IntegralUniversalSegPnt(rut1[0], rut1[1], x, f0, kernel_param3, IntParam_3, res_seg_3);
    std::cout << std::endl;
    std::cout << "Test3. Segment integral of a vector function f0" << std::endl;
    std::cout << "       Segment integral AB: ";
    std::cout << res_seg_3[0] << " " << res_seg_3[1] << " " << res_seg_3[2] << std::endl;
    IntegralUniversalSegPnt(rut1[1], rut1[0], x, f0, kernel_param3, IntParam_3, res_seg_3);
    std::cout << "       Segment integral BA: ";
    std::cout << res_seg_3[0] << " " << res_seg_3[1] << " " << res_seg_3[2] << std::endl;
    
    return 0;
}