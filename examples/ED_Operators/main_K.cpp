#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>


#include "kernel_param.h"
#include "K.h"
#include "K0.h"
#include "R.h"
#include "num_param.h"
#include "element_geom.h"
#include "kernel_lib.h"


/**
    * Тест для проверки вычисления основных интегральных операторов задач
    * рассеяния:
    *               K = rot rot [surf_int(j * F)]
    *                   F = e^{ikr} / (4 pi r)
    * Ток должен лежить в плоскости ячейки(если в операторе присутствует
    * контурный интеграл).
    * 
    * 
    * Вариации:
    *       1) K_RotRot_Near_Smooth - оператор K[] в ближней зоне со сглаживанием
    * в обоих внутренних интегралах. Точка может лежать прямо на ячейке.
    * Использовать: при расчете поля в ближней зоне.
    * Обязательная инициализация следующих параметров: 
    *       const double IntegralAccuracy = 0.00001;
    *       const double SmoothingR_r = 0.0005;                // relative to cell diam h
    *       const double SmoothingRSeg_r = 0.0005;             // relative to cell diam h
    *       const int NCellStart = 8; 
    *       const int NSegStart = 10;
    *       const int PMax = 5;                                // 2^{PMax} steps in surface integration
    *       const int PMaxSeg = 5;                             // 2^{PMaxSeg} steps in segment integration
    * 
    * 
    * 
    *       2) K_RotRot_HS - оператор K[] в ближней зоне с выделением особенности.
    * Точка лежит на ячейке или просто на r < 3h от ячейки. Сглаживание в контурном
    * интеграле.
    * Использовать: при расчете поля в ближней зоне.
    * Обязательная инициализация следующих параметров: 
    *       const double IntegralAccuracy = 0.00001;
    *       const double SmoothingRSeg_r = 0.0005;             // relative to cell diam h
    *       const double AnalyticCalcR = 3.0;                  // relative to cell diam h
    *       const int NCellStart = 8; 
    *       const int NSegStart = 10;
    *       const int PMax = 5;                                
    *       const int PMaxSeg = 5;                             
    * 
    * 
    * 
    *       3) K_RotRot_HS_Colloc - оператор K[] в точках коллокации с выделением
    * особенности. Контурный интеграл без сглаживания.
    * Использовать: при наборе матрицы системы в точках коллокации.
    * Обязательная инициализация следующих параметров:
    *       const double IntegralAccuracy = 0.00001;
    *       const double AnalyticCalcR = 3.0;                  // relative to cell diam h
    *       const int NCellStart = 8; 
    *       const int NSegStart = 10;
    *       const int PMax = 5;                                
    *       const int PMaxSeg = 5;   
    * 
    * 
    * 
    * 
    *       3) K_RotRot_Far - оператор K[] в дальней зоне со специальным ядром.
    * Точка лежит далеко от поверхности (например, > 3h). Быстрый вариант для 
    * вычисления поля вдали от поверхности.
    * Обязательная инициализация следующих параметров: 
    *       const double IntegralAccuracy = 0.00001;
    *       const int NCellStart = 8; 
    *       const int NSegStart = 10;
    *       const int PMax = 5;                                
    *       const int PMaxSeg = 5;                             
    * 
    * 
    * 
    * 
    * 
    * 
    * Notes:
    * SmoothingR_r важно задавать ненулевым (0.0005), если точка x попадает на ячейку.
    * SmoothingRSeg_r важно задавать ненулевым (0.0005), если точка может попасть
    * на край ячейки (при наборе матрицы - не может, следовательно SmoothingRSeg_r = 10e-16).
    * 
    * При решении задачи все численные параметры задаются универсальными один раз.
    * Внутри каждого из вариантов оператора K[] происходит обращение к
    * необходимым численным параметрам.
    *        
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
    
    
    const double rut3[3][3] = {{3., 2., 0.},
                               {2., 0., 0.},
                               {0., 2., 0.}};


    const double x[3] = {10., 10., 10.};      // all three options will match
    //const double x[3] = {2., 1.5, 0.};      // K_RotRot_Far won't work right in this case



    //==============================================
    //--------Global numerical parameters-----------
    //==============================================
    const double IntegralAccuracy = 0.00001;
    const double SmoothingR_r = 0.0005;                // relative to cell diam h
    const double SmoothingRSeg_r = 0.0005;             // relative to cell diam h
    const double AnalyticCalcR = 3.0;                  // relative to cell diam h
    const int NCellStart = 8;                       
    const int NSegStart = 10;
    const int PMax = 5;                                // 2^{PMax} steps in surface integration
    const int PMaxSeg = 5;                             // 2^{PMaxSeg} steps in segment integration
    NumParam num_param(IntegralAccuracy, SmoothingR_r, SmoothingRSeg_r,
                AnalyticCalcR, NCellStart, NSegStart, PMax, PMaxSeg);


    const std::complex<double> k = std::complex<double>(1., 0.);
    const std::complex<double> j[3] = {1., 0., 0.,};




    //==============================================
    //-------------------Test1---------------------+
    //-------K[] в ближней зоне со сглаживанием-----
    //---------------(расчет поля)------------------
    //==============================================
    std::complex<double> res[3];
    K_RotRot_Near_Smooth(j, x, rut1, num_param, k, res);
    std::cout << "Test1. K[] со сглаживанием." << std::endl;
    std::cout << "      [4][4]: " << res[0] << " " << res[1] << " " << res[2] << std::endl;


    std::complex<double> res1[3], res2[3];
    K_RotRot_Near_Smooth(j, x, rut2, num_param, k, res1);
    K_RotRot_Near_Smooth(j, x, rut3, num_param, k, res2);
    std::cout << "      [3][3] + [3][3]: " << res1[0] + res2[0] << " " << res1[1] + res2[1] << " " << res1[2] + res2[2] << std::endl;





    //==============================================
    //-------------------Test2---------------------+
    //-------------K[] в ближней зоне---------------
    //------------с выделением особенности----------
    //---------------(вычисление поля)--------------
    //==============================================
    std::complex<double> res3[3];
    K_RotRot_HS(j, x, rut1, num_param, k, res3);
    std::cout << "Test2. K в ближней зоне с выделением особенности" << std::endl;
    std::cout << "      [4][4]: " << res3[0] << " " << res3[1] << " " << res3[2] << std::endl;



    std::complex<double> res13[3], res23[3];
    K_RotRot_HS(j, x, rut2, num_param, k, res13);
    K_RotRot_HS(j, x, rut3, num_param, k, res23);
    std::cout << "      [3][3] + [3][3]: " << res13[0] + res23[0] << " " << res13[1] + res23[1] << " " << res13[2] + res23[2] << std::endl;
    



    //==============================================
    //-------------------Test3---------------------+
    //--------------K[] в ближней зоне--------------
    //-------------в точках коллокаций--------------
    //------------с выделением особенности----------
    //----------------(набор матрицы )--------------
    //==============================================
    K_RotRot_HS_Colloc(j, x, rut1, num_param, k, res3);
    std::cout << "Test3. K при наборе матрицы в точках коллокаций." << std::endl;
    std::cout << "      [4][4]: " << res3[0] << " " << res3[1] << " " << res3[2] << std::endl;



    K_RotRot_HS_Colloc(j, x, rut2, num_param, k, res13);
    K_RotRot_HS_Colloc(j, x, rut3, num_param, k, res23);
    std::cout << "      [3][3] + [3][3]: " << res13[0] + res23[0] << " " << res13[1] + res23[1] << " " << res13[2] + res23[2] << std::endl;
    



    //==============================================
    //-------------------Test4---------------------+
    //-----------K[] в дальней зоне со--------------
    //-------------специальным ядром----------------
    //-------(вычисление поля в дальней зоне)-------
    //==============================================
    K_RotRot_Far(j, x, rut1, num_param, k, res);
    std::cout << "Test4. K в дальней зоне" << std::endl;
    std::cout << "      [4][4]: " << res[0] << " " << res[1] << " " << res[2] << std::endl;

    K_RotRot_Far(j, x, rut2, num_param, k, res1);
    K_RotRot_Far(j, x, rut3, num_param, k, res2);
    std::cout << "      [3][3] + [3][3]: " << res1[0] + res2[0] << " " << res1[1] + res2[1] << " " << res1[2] + res2[2] << std::endl;


    return 0;
}