#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>


#include "kernel_param.h"
#include "R.h"
#include "num_param.h"
#include "element_geom.h"
#include "kernel_lib.h"


/**
    * Тест для проверки вычисления основных интегральных операторов задач
    * рассеяния:
    *               R = rot [surf_int(j * F)]
    *                   F = e^{ikr} / (4 pi r)
    * 
    * Вариации:
    *       - R_Rot_Colloc - оператор R[] без сглаживания в точках коллокации.
    * Оператор равен нулю при совпадении точки коллокации (точки расчета)
    * и центра ячейки, иначе - вычисление без сглаживания. 
    * Использовать: при наборе матрицы в точках коллокации.
    * 
    * Обязательная инициализация следующих параметров: 
    *       const double IntegralAccuracy = 0.00001;
    *       const int NCellStart = 8;
    *       const int NSegStart = 10;
    *       const int PMax = 10;                                // 2^{PMax} steps in surface integration
    *       const int PMaxSeg = 10;                             // 2^{PMaxSeg} steps in segment integration
    * 
    * 
    * 
    * 
    *       - R_Rot_Smooth - оператор R[] со сглаживанием в произвольной точке.
    * Использовать: при расчете поля.
    * Обязательная инициализация следующих параметров:
    *       const double IntegralAccuracy = 0.00001;
    *       const double SmoothingR_r = 0.0005;                 // relative to cell diam h
    *       const int NCellStart = 8;
    *       const int NSegStart = 10;
    *       const int PMax = 10;                                // 2^{PMax} steps in surface integration
    *       const int PMaxSeg = 10;                             // 2^{PMaxSeg} steps in segment integration
*/

using namespace bielcc;

int main(int argc, char **argv)
{
    //==============================================
    //--------Cells for operators testing-----------
    //==============================================
    const double rut1[4][3] = {{0., 0., 0.},
                               {2., 0., 0.},
                               {3., 2., 0.},
                               {0., 2., 0.}};

    const double rut2[3][3] = {{0., 2., 0.},
                               {0., 0., 0.},
                               {2., 0., 0.}};
    
    
    const double rut3[3][3] = {{0., 2., 0.},
                               {2., 0., 0.},
                               {3., 2., 0.}};


    //const double x[3] = {10., 10., 10.}; //+++++
    const double x[3] = {2., 1.5, 0.};




    //==============================================
    //--------Global numerical parameters-----------
    //==============================================
    const double IntegralAccuracy = 0.00001;
    const double SmoothingR_r = 0.005;                // relative to cell diam h
    const double SmoothingRSeg_r = 0.005;             // relative to cell diam h
    const double AnalyticCalcR = 3.0;                  // relative to cell diam h
    const int NCellStart = 8;                       
    const int NSegStart = 10;
    const int PMax = 10;                                // 2^{PMax} steps in surface integration
    const int PMaxSeg = 10;                             // 2^{PMaxSeg} steps in segment integration



    NumParam num_param_smooth(IntegralAccuracy, SmoothingR_r, SmoothingRSeg_r,
                AnalyticCalcR, NCellStart, NSegStart, PMax, PMaxSeg); // параметры со сглаживанием

    NumParam num_param(IntegralAccuracy, 10e-16, 10e-16,
                AnalyticCalcR, NCellStart, NSegStart, PMax, PMaxSeg); // параметры без сглаживания



    const std::complex<double> j[3] = {1., 1., 1.,};
    //const std::complex<double> j[3] = {1., 0., 0.};
    const std::complex<double> k = std::complex<double>(1., 0.);


    //==============================================
    //-------------------Test1---------------------+
    //-------------R[] в произвольной точке---------
    //---------------со сглаживанием----------------
    //-----------------(расчет поля)----------------
    //==============================================
    std::complex<double> res[3]{};
    R_Rot_Smooth(j, x, rut1, num_param_smooth, k, res); // Сглаживание нужно, точка на ячейке
    std::cout << "Test1. R в произвольной точке, сглаживание." << std::endl;
    std::cout << "      [4][4]: " << res[0] << " " << res[1] << " " << res[2] << std::endl;


    std::complex<double> res1[3]{}, res2[3]{};
    R_Rot_Smooth(j, x, rut2, num_param, k, res1); // Сглаживание не нужно, точка вне ячейки
    R_Rot_Smooth(j, x, rut3, num_param_smooth, k, res2); // Сглаживание нужно, точка на ячейке
    std::cout << "      [3][3] + [3][3]: " << res1[0] + res2[0] << " " << res1[1] + res2[1] << " " << res1[2] + res2[2] << std::endl;




    //==============================================
    //-------------------Test2---------------------+
    //-----------R[] в точках коллокации------------
    //------------внутри сглаживания нет------------
    //==============================================
    const double x1_colloc[3] = {1.25, 1., 0.};
    R_Rot_Colloc(j, x1_colloc, rut1, num_param, k, res);
    std::cout << "Test2. R в точке коллокации" << std::endl;
    std::cout << "      [4][4]: " << res[0] << " " << res[1] << " " << res[2] << std::endl;


    const double x2_colloc[3] = {2. / 3., 2. / 3., 0.};
    R_Rot_Colloc(j, x2_colloc, rut2, num_param, k, res1);
    std::cout << "      [3][3]: " << res1[0] << " " << res1[1] << " " << res1[2] << std::endl;

    R_Rot_Colloc(j, x2_colloc, rut3, num_param, k, res1);
    std::cout << "      [3][3] (not colloc pnt): " << res1[0] << " " << res1[1] << " " << res1[2] << std::endl;


    return 0;
}
