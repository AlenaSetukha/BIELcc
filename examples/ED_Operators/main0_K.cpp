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
#include "constants.h"


/**
    * Тест для проверки вычисления основных интегральных операторов задач
    * рассеяния:
    *               K = rot rot [surf_int(j * F)]
    *                   F = e^{ikr} / (4 pi r)
    * Ток должен лежать в плоскости ячейки (если в операторе присутствует
    * контурный интеграл).
    * 
    * 
    * Вариации:
    *       1) K_RotRot_Near_Smooth - оператор K[] в ближней зоне со сглаживанием
    * в обоих внутренних интегралах. Точка может лежать прямо на ячейке.
    * Использовать: при расчете поля в ближней зоне.
    * Обязательная инициализация следующих параметров: 
    *       const double IntegralAccuracy = 0.00001;
    *       const double SmoothingR_r = 0.5;                   // relative to small cell diam h2
    *       const double SmoothingRSeg_r = 0.5;                // relative to small cell diam h2
    *       const int NCellStart = 20; 
    *       const int NSegStart = 30;
    *       const int PMax = 2;                                // 2^{PMax} steps in surface integration
    *       const int PMaxSeg = 2;                             // 2^{PMaxSeg} steps in segment integration
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


    // const double x[3] = {10., 10., 10.};      // all three options will match
    const double x[3] = {2., 1.5, 0.};      // K_RotRot_Far will not give the correct result in this case


    
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
    
    std::cout << "Numerical values ​​of calculation parameters BY DEFAULT" << std::endl;
    std::cout << "      Integrals calculation accuracy: " << IntegralAccuracy << std::endl;
    std::cout << "      Relative smoothing radius (surface, rel. to h2): " << SmoothingR_r << std::endl;
    std::cout << "      Relative smoothing radius (segment, rel. to h2): " << SmoothingRSeg_r << std::endl;
    std::cout << "      Radius of analytical calculation (rel. to h): " << AnalyticCalcR << std::endl;
    std::cout << "      Starting cell split: " << NCellStart << std::endl;
    std::cout << "      Starting segment split: " << NSegStart << std::endl;
    std::cout << "      Limit cell split (2^{P}): " << PMaxCell << std::endl;
    std::cout << "      Limit segment split (2^{P}): " << PMaxSeg << std::endl;
    std::cout << std::endl;
 

    
    NumParam num_param(IntegralAccuracy,
                        SmoothingR_r,
                        SmoothingRSeg_r,
                        AnalyticCalcR,
                        NCellStart * 10,
                        NSegStart * 10,
                        PMaxCell * 6,
                        PMaxSeg * 6);

    NumParam num_param_no_smooth(IntegralAccuracy,
                        Calculation_Constants::MACHINE_ZERO,
                        Calculation_Constants::MACHINE_ZERO,
                        AnalyticCalcR,
                        NCellStart * 10,
                        NSegStart * 10,
                        PMaxCell * 6,
                        PMaxSeg * 6);


    const std::complex<double> k = std::complex<double>(1., 0.);
    const std::complex<double> j[3] = {1., 0., 0.,};




    //==============================================
    //-------------------Test1---------------------+
    //-------K[] в ближней зоне со сглаживанием-----
    //---------------(расчет поля)------------------
    //==============================================
    std::complex<double> res[3];
    K_RotRot_Near_Smooth(j, x, rut1, num_param, k, res);             // x on cell -> smooth (both terms)
    std::cout << "Test1. K[] со сглаживанием." << std::endl;
    std::cout << "      [4][4]: " << res[0] << " " << res[1] << " " << res[2] << std::endl;


    std::complex<double> res1[3], res2[3];
    K_RotRot_Near_Smooth(j, x, rut2, num_param_no_smooth, k, res1);  // x not on cell -> no smooth
    K_RotRot_Near_Smooth(j, x, rut3, num_param, k, res2);            // x on cell -> smooth (both terms)
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
    K_RotRot_HS(j, x, rut2, num_param_no_smooth, k, res13);
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
    double rkt1[3]{}, rkt2[3]{}, rkt3[3]{};
    double cell_diam1 = get_diam(rut1);
    get_center_mass(rut1, rkt1);
    get_center_mass(rut2, rkt2);
    get_center_mass(rut3, rkt3);
    if (dist(x, rkt1) > 3. * cell_diam1 &&
        dist(x, rkt2) > 3. * cell_diam1 &&
        dist(x, rkt3) > 3. * cell_diam1)
    {
        K_RotRot_Far(j, x, rut1, num_param, k, res);
        std::cout << "Test4. K в дальней зоне (на ячейке не сработает!)" << std::endl;
        std::cout << "      [4][4]: " << res[0] << " " << res[1] << " " << res[2] << std::endl;

        K_RotRot_Far(j, x, rut2, num_param, k, res1);
        K_RotRot_Far(j, x, rut3, num_param, k, res2);
        std::cout << "      [3][3] + [3][3]: " << res1[0] + res2[0] << " " << res1[1] + res2[1] << " " << res1[2] + res2[2] << std::endl;


    }
    return 0;
}