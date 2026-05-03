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
    * Тест для проверки вычисления основных интеграторов по 
    * четырехугольным и треугольным ячейкам.
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
    //----------Integrand parameters----------------
    //==============================================
    KernelParam<double> param1, param2, param3;
    param1.smoothR = Calculation_Constants::SMOOTHING_DIST *
                                            get_diam(rut1); // absolute value of smoothing radius
    param2.smoothR = Calculation_Constants::SMOOTHING_DIST *
                                            get_diam(rut2); // absolute value of smoothing radius
    param3.smoothR = Calculation_Constants::SMOOTHING_DIST *
                                            get_diam(rut3); // absolute value of smoothing radius
    param1.k = static_cast<double>(1.0);
    param2.k = static_cast<double>(1.0);
    param3.k = static_cast<double>(1.0);



    //==============================================
    //--------Integral calculation parameters-------
    //==============================================
    IntegralParam int_param_dim1(1, 8, 10, 0.0000001);
    IntegralParam int_param_dim3(3, 8, 10, 0.0000001);
    





    //==============================================
    //-------------------Test0---------------------+
    //-----1D surface integral with integrand f1----
    //==============================================
    double s_trap[1];
    IntegralUniversal(rut1, f1, param1, int_param_dim1, s_trap);
    std::cout << "Test0. Area of ​​a trapezoid 1 (Integral [4][3]): ";
    std::cout << s_trap[0] << std::endl;



    //==============================================
    //-------------------Test1---------------------+
    //---------1D surface integral in pnt x---------
    //----------with integrand F_SimplePot_H-------- 
    //==============================================
    std::complex<double> res1d[1];
    IntegralUniversalPnt(x, rut1, F_SimplePot_H, param1, int_param_dim1, res1d);
    std::cout << "Test1. Surface integral of a F_SimplePot_H" << std::endl;
    std::cout << "       Integral [4][3] = " << res1d[0] / 4. / M_PI  << std::endl;
    


    std::complex<double> res1d_2[1], res1d_3[1];
    IntegralUniversalPnt(x, rut2, F_SimplePot_H, param2, int_param_dim1, res1d_2);
    IntegralUniversalPnt(x, rut3, F_SimplePot_H, param3, int_param_dim1, res1d_3);
    std::cout << "       Integral [3][3] + [3][3]: ";
    std::cout << (res1d_2[0] + res1d_3[0]) / 4. / M_PI  << std::endl;

    





    //==============================================
    //-------------------Test2---------------------+
    //---------3D surface integral in pnt x---------
    //---------------with integrand f0-------------- 
    //==============================================
    double res3d[3];
    IntegralUniversalPnt(x, rut1, f0, param1, int_param_dim3, res3d);
    std::cout << "Test2. Surface integral of a vector function f0" << std::endl;
    std::cout << "       Integral [4][3]: ";
    std::cout << res3d[0] << " " << res3d[1] << " " << res3d[2] << std::endl;
    


    double res3d_2[3], res3d_3[3];
    IntegralUniversalPnt(x, rut2, f0, param2, int_param_dim3, res3d_2);
    IntegralUniversalPnt(x, rut3, f0, param3, int_param_dim3, res3d_3);
    std::cout << "       Integral [3][3] + [3][3]: " << res3d_2[0] + res3d_3[0] << " " << 
                    res3d_2[1] + res3d_3[1] << " " << res3d_2[2] + res3d_3[2] << std::endl;






    //==============================================
    //-------------------Test3---------------------+
    //---------3D segment integral in pnt x---------
    //---------------with integrand f0--------------
    //==============================================
    double res_seg[3];
    int_param_dim3.SetNStart(4);
    int_param_dim3.SetPMax(2);
    int_param_dim3.SetEpsAccur(0.0001);
    param1.smoothR = Calculation_Constants::MACHINE_ZERO; // x does not lies on the cell edge
    IntegralUniversalSegPnt(rut1[0], rut1[1], x, f0, param1, int_param_dim3, res_seg);
    std::cout << "Test3. Segment integral of a vector function f0" << std::endl;
    std::cout << "       Segment integral AB: ";
    std::cout << res_seg[0] << " " << res_seg[1] << " " << res_seg[2] << std::endl;
    IntegralUniversalSegPnt(rut1[1], rut1[0], x, f0, param1, int_param_dim3, res_seg);
    std::cout << "       Segment integral BA: ";
    std::cout << res_seg[0] << " " << res_seg[1] << " " << res_seg[2] << std::endl;





    //==============================================
    //-------------------Test4---------------------+
    //------Semi-Analytical 1D surface integral-----
    //-----------in pnt x with integrand------------
    //------------F(x, y) = e^{ikr} / r-------------
    //==============================================
    double k_wave = 1.;
    double analyticR_r = 3.0;
    std::complex<double> sa = IntegralSurf_SimplePot_H_HS(rut1, x, analyticR_r,
                                                   int_param_dim1, k_wave);
    std::cout << "Test4. Highlighting singularity in e^{ikr} / r" << std::endl;
    std::cout << "Semi analytical: " << sa / 4. / M_PI << std::endl;



    std::complex<double> calc_res[0];
    IntegralUniversalPnt(x, rut1, F_SimplePot_H, param1, int_param_dim1, calc_res);
    std::cout << "Fully numerical: " << calc_res[0] / 4. / M_PI  << std::endl;

    return 0;
}