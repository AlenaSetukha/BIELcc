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
#include "constants.h"


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
    * Использовать: при расчете поля. ВНИМАНИЕ: гиперсингулярный интеграл, считаеть аккуратно.
    * Обязательная инициализация следующих параметров:
    *       const double IntegralAccuracy = 0.00001;
    *       const double SmoothingR_r = 8.;                     // relative to small cell (h2)
    *       const int NCellStart = 40;
    *       const int NSegStart = 10;
    *       const int PMax = 1;                                // 2^{PMax} steps in surface integration
    *       const int PMaxSeg = 1;                             // 2^{PMaxSeg} steps in segment integration
    * 
    * Важно: адаптивное интегрирование лучше не использовать, ячейку сразу бить 
    * подробно, радиус сглаживания - не задавать ммаленьким (~8-10 длин ячеек второго уровня).
    * Не допускать слишком маленького радиуса сглаживания.
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


    // const double x[3] = {10., 10., 10.}; //+++++
    const double x[3] = {2., 1.5, 0.};




    //==============================================
    //--------Global numerical parameters-----------
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
 


    NumParam num_param_smooth(IntegralAccuracy,
                        SmoothingR_r * 1.5,
                        SmoothingRSeg_r * 1.5,
                        AnalyticCalcR,
                        NCellStart * 20,
                        NSegStart * 40,
                        PMaxCell * 4,
                        PMaxSeg * 4);

    NumParam num_param_no_smooth(IntegralAccuracy,
                        Calculation_Constants::MACHINE_ZERO,
                        Calculation_Constants::MACHINE_ZERO,
                        AnalyticCalcR,
                        NCellStart * 10,
                        NSegStart * 10,
                        PMaxCell * 4,
                        PMaxSeg * 4);

    const std::complex<double> j[3] = {1., 1., 1.,};
    // const std::complex<double> j[3] = {1., 2., 0.};
    const std::complex<double> k = std::complex<double>(1., 0.);


    //==============================================
    //-------------------Test1---------------------+
    //-------------R[] в произвольной точке---------
    //---------------со сглаживанием----------------
    //-----------------(расчет поля)----------------
    //==============================================
    std::complex<double> res[3]{};
    R_Rot_Smooth(j, x, rut1, num_param_smooth, k, res);           // Сглаживание нужно, точка на ячейке
    std::cout << "Test1. R в произвольной точке, сглаживание." << std::endl;
    std::cout << "      [4][4]: " << res[0] << " " << res[1] << " " << res[2] << std::endl;

    
    std::complex<double> res1[3]{}, res2[3]{};
    R_Rot_Smooth(j, x, rut2, num_param_no_smooth, k, res1);        // Сглаживание не нужно, точка вне ячейки
    R_Rot_Smooth(j, x, rut3, num_param_smooth, k, res2);           // Сглаживание нужно, точка на ячейке
    std::cout << "      [3][3] + [3][3]: " << res1[0] + res2[0] << " " << res1[1] + res2[1] << " " << res1[2] + res2[2] << std::endl;

    
    R_Rot_Near(j, x, rut1, num_param_smooth, k, res);
    std::cout << "      [4][4] (curl + surf): " << res[0] << " " << res[1] << " " << res[2] << std::endl;




    //==============================================
    //-------------------Test2---------------------+
    //-----------R[] в точках коллокации------------
    //------------(внутри сглаживания нет)------------
    //==============================================
    const double x1_colloc[3] = {1.25, 1., 0.};
    R_Rot_Colloc(j, x1_colloc, rut1, num_param_no_smooth, k, res);
    std::cout << "Test2. R в точке коллокации" << std::endl;
    std::cout << "      [4][4]: " << res[0] << " " << res[1] << " " << res[2] << std::endl;


    const double x2_colloc[3] = {2. / 3., 2. / 3., 0.};
    R_Rot_Colloc(j, x2_colloc, rut2, num_param_no_smooth, k, res1);
    std::cout << "      [3][3]: " << res1[0] << " " << res1[1] << " " << res1[2] << std::endl;




    //==============================================
    //-------------------Test3---------------------+
    //-----------R[] в точках коллокации------------
    //-----------(внутри сглаживания нет)-----------
    //==============================================
    std::cout << std::endl;
    std::cout << "Test3. R_Rot_Colloc не в точке коллокации" << std::endl;
    R_Rot_Colloc(j, x2_colloc, rut3, num_param_no_smooth, k, res1); // not colloc point for rut3 -> calculation with no smoothing
    std::cout << "      [3][3] (not colloc pnt): " << res1[0] << " " << res1[1] << " " << res1[2] << std::endl;
    R_Rot_Smooth(j, x2_colloc, rut3, num_param_no_smooth, k, res1);
    std::cout << "      [3][3] (numerically): " << res1[0] << " " << res1[1] << " " << res1[2] << std::endl;

    return 0;
}
