#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>
#include <omp.h>

#include "radiation_source.h"
#include "surf_mesh.h"
#include "io_aereco.h"
#include "io_vtk.h"
#include "K.h"
#include "R.h"
#include "material.h"
#include "num_param.h"
#include "incident_field.h"


/**
    * Тестовый пример расчета электромагнитного поля идеального проводника
    * в точках некоторой плоскости (без набора матрицы).
    * 
    * 1. Выбирается сетка (ZB: CellQuad) объекта идеального проводника (замкнутая),
    * а также точки плоскости, в которых мы хотим посчитать поле.
    * 2. Электрические и магнитные токи задаются на ячейках разбиения
    * поверхности как:
    *                  j_E(x_i) = -n x H_inc(x_i)
    *                  j_M(x_i) = n x E_inc(x_i)
    * 3. Поле в заданных точках вычисляется двумя способами: 
    *      - напрямую: E(x) = E_inc(x), H(x) = H_inc(x)
    *      - по формулам Стреттона-Чу через построенные токи {j_E, j_M}
    * 
    * 4. Ответ: поле по формулам С-Ч должно получиться E_inc внутри
    * области пересечения с объектом (например, со сферой), вне - ноль.
*/ 


using namespace bielcc;

int main(int argc, char **argv)
{
    const std::string geom_fname = argc > 1 ? argv[1] :
                                    "../examples/data/surface_mesh/sphere_quad.vtk", \
                      result_dir = argc > 2 ? argv[2] :
                                                "../", \
                      pntsForField_fname = argc > 3 ? argv[3] :
                                    "../examples/data/surface_mesh/plate_oxy.vtk";
    auto start = std::chrono::high_resolution_clock::now();



    //===============1. Mesh object creation========================
    bielcc::SurfMesh<CellQuad> surf_mesh(geom_fname);
    double grid_step = surf_mesh.GetGridStep();
    std::cout << "grid step: " << grid_step << std::endl;
    std::cout << "surface square: " << surf_mesh.GetSurfSquare(0) << std::endl;




    //===============2. Generation required cell data===============
    std::vector<std::array<double, 3>> rkt = surf_mesh.FormCollocPoints();
    std::vector<std::array<double, 3>> norm = surf_mesh.FormNormalVectors();
    int NFrm = surf_mesh.GetNFrm();





    //===============3. Collocation points record to file===========
    IO_VTK io_vtk;
    std::string colloc_pnt_fname = result_dir + "body.vtk";
    io_vtk.WriteVecs(colloc_pnt_fname, rkt);
    


    //===============4. Forming EM parameters for one frequency=====
    double freq = 47713465.2095072; // in Hz for k = 1
    double eps_r = 1.;
    double mu_r = 1.;


    MaterialParam vacuum_material = MaterialParam("Vacuum",
                          ChangeMaterialTypeDir::Isotropic,
                     ChangeMaterialTypeFreq::NonDispersive,
                                  eps_r, mu_r, 0., 0., 0.);
    vacuum_material.Calc(freq);
    


    double k_ort[3] = {-1., 0., 0.};
    double e0_ort[3] = {0., 1., 0.}; // VV polar

    RadiationSource src = RadiationSource(RadiationSourceType::PlaneWave,
                                                           k_ort, e0_ort);

    src.CalcK(freq, vacuum_material); // k во внешней среде
    double k = src.GetK();
    double k_vec[3];
    src.GetKVec(k_vec);
    std::cout << "Волновое число в вакууме: " << k << std::endl;
    std::cout << "Волновой вектор в вакууме: " << k_vec[0] << " " <<
                                k_vec[1] << " " << k_vec[2] << std::endl;



    




    //===============6. Points for field reading==================
    std::vector<std::array<double, 3>> points_for_field = 
                            io_vtk.ReadVecs(pntsForField_fname);


    



    //===============7. Field initialization directly===============
    int num_points = points_for_field.size();
    std::cout << "Число точек для расчета поля: " << num_points << std::endl;
    std::vector<std::array<std::complex<double>, 3>> field_E(num_points); 
    std::vector<std::array<std::complex<double>, 3>> field_H(num_points);

    KernelParam<double> EIncParam;
    EIncParam.k = k;
    for (int i = 0; i < 3; i++) {
        EIncParam.vec_dbl[i] = k_vec[i];
        EIncParam.e0[i] = e0_ort[i];
        EIncParam.freq = freq;
    }

    for (int i = 0; i < num_points; i++) {
        PlaneWaveFieldE(points_for_field[i].data(), EIncParam, field_E[i].data());
        PlaneWaveFieldH(points_for_field[i].data(), EIncParam, vacuum_material, field_H[i].data());
    }

    std::string res_fname = result_dir + "field_E.vtk";
    io_vtk.WritePntsWithField_C(res_fname, points_for_field, field_E);
    res_fname = result_dir + "field_H.vtk";
    io_vtk.WritePntsWithField_C(res_fname, points_for_field, field_H);








    //===============8. Initialization of currents on cells ========
    std::vector<std::array<std::complex<double>, 3>> j_E(NFrm);
    std::vector<std::array<std::complex<double>, 3>> j_M(NFrm);
    std::complex<double> H_inc[3]{}, E_inc[3]{}, vp[3]{};
    for (int i = 0; i < NFrm; i++) {
        PlaneWaveFieldH(rkt[i].data(), EIncParam, vacuum_material, H_inc);
        vec_prod(norm[i].data(), H_inc, vp);
        j_E[i][0] = -vp[0], j_E[i][1] = -vp[1], j_E[i][2] = -vp[2]; // n x Hinc


        PlaneWaveFieldE(rkt[i].data(), EIncParam, E_inc);
        vec_prod(norm[i].data(), E_inc, vp);
        j_M[i][0] = vp[0], j_M[i][1] = vp[1], j_M[i][2] = vp[2]; // -n x Einc
    }


    // Writing surface currents on cells in vtk files
    std::string j_fname = result_dir + "j_E.vtk";
    io_vtk.WritePntsWithField_C(j_fname, rkt, j_E);
    j_fname = result_dir + "j_M.vtk";
    io_vtk.WritePntsWithField_C(j_fname, rkt, j_M);








    //===============9. Field by Stratton-Chu formulas==============
    std::complex<double> multE = std::complex<double>(0., 1.) /
                                                (2. * M_PI * freq);
    std::complex<double> multH = std::complex<double>(0., 1.) /
                                                (2. * M_PI * freq); 
    try {
        multE /= std::get<std::complex<double>>(vacuum_material.GetEps());
    } catch (const std::bad_variant_access& e) {
        multE /= std::get<double>(vacuum_material.GetEps());
    }

    try {
        multH /= std::get<std::complex<double>>(vacuum_material.GetMu());
    } catch (const std::bad_variant_access& e) {
        multH /= std::get<double>(vacuum_material.GetMu());
    }


    std::cout << "multE = " << multE << std::endl;
    std::cout << "multH = " << multH << std::endl;

    multE = 1.;
    multH = 1.;

    std::vector<std::array<std::complex<double>, 3>> field_E_SCH(num_points);
    std::vector<std::array<std::complex<double>, 3>> field_H_SCH(num_points);




    //==============================================
    //--------Global numerical parameters-----------
    //==============================================
    const double IntegralAccuracy = 0.000001;
    const double SmoothingR_r = 0.5;                // relative to small cell diam h2
    const double SmoothingRSeg_r = 0.5;             // relative to small cell diam h2
    const double AnalyticCalcR = 3.0;                   // relative to cell diam h
    const int NCellStart = 40;                       
    const int NSegStart = 60;
    const int PMax = 1;                                // 2^{PMax} steps in surface integration
    const int PMaxSeg = 1;                             // 2^{PMaxSeg} steps in segment integration
    
    NumParam num_param(IntegralAccuracy, SmoothingR_r, SmoothingRSeg_r,
                                AnalyticCalcR, NCellStart, NSegStart, PMax, PMaxSeg);

    omp_set_dynamic(0);
    omp_set_num_threads(18);
    #pragma omp parallel for
    for (int i = 0; i < num_points; i++) {
        std::complex<double> sum_KE[3]{}, sum_RE[3]{};
        std::complex<double> sum_KH[3]{}, sum_RH[3]{};
        std::complex<double> cur_res3[3]{};

        for (int j = 0; j < NFrm; j++) {
            double pnt_dist = dist(points_for_field[i].data(), rkt[j].data());
            CellQuad cell_j = surf_mesh.GetCell(j);


            if (pnt_dist > AnalyticCalcR * grid_step) { // ДАЛЬНЯЯ ЗОНА (dist > 3h)
                // Если точка далеко от текущей ячейки:
                //          -> K[] спец. вида
                //          -> R_Rot_Colloc или R_Rot_Smooth (оба подходят). No smooth

                // Field E
                K_RotRot_Far(j_E[j].data(), points_for_field[i].data(),
                            cell_j.GetCellCoords(), num_param, k, cur_res3);
                sum_KE[0] += cur_res3[0];
                sum_KE[1] += cur_res3[1];
                sum_KE[2] += cur_res3[2];


                R_Rot_Colloc(j_M[j].data(), points_for_field[i].data(),
                             cell_j.GetCellCoords(), num_param, k, cur_res3);
                sum_RE[0] += cur_res3[0];
                sum_RE[1] += cur_res3[1];
                sum_RE[2] += cur_res3[2];




                // Field H
                K_RotRot_Far(j_M[j].data(), points_for_field[i].data(), 
                                cell_j.GetCellCoords(), num_param, k, cur_res3);
                sum_KH[0] += cur_res3[0];
                sum_KH[1] += cur_res3[1];
                sum_KH[2] += cur_res3[2];

                R_Rot_Colloc(j_E[j].data(), points_for_field[i].data(),
                               cell_j.GetCellCoords(), num_param, k, cur_res3);
                sum_RH[0] += cur_res3[0];
                sum_RH[1] += cur_res3[1];
                sum_RH[2] += cur_res3[2];

            } else { 
                // БЛИЖНЯЯ ЗОНА (dist < 3h)
                // Если близко к текущей ячейке:
                //          -> K[] с выделением особенности и со сглаживанием
                //          -> R[] со сглаживанием
                // Field E
                K_RotRot_HS(j_E[j].data(), points_for_field[i].data(), 
                            cell_j.GetCellCoords(), num_param, k, cur_res3);
                sum_KE[0] += cur_res3[0];
                sum_KE[1] += cur_res3[1];
                sum_KE[2] += cur_res3[2];



                R_Rot_Smooth(j_M[j].data(), points_for_field[i].data(), 
                        cell_j.GetCellCoords(), num_param, k, cur_res3);
                sum_RE[0] += cur_res3[0];
                sum_RE[1] += cur_res3[1];
                sum_RE[2] += cur_res3[2];



                // Field H
                K_RotRot_HS(j_M[j].data(), points_for_field[i].data(), 
                      cell_j.GetCellCoords(), num_param, k, cur_res3);
                sum_KH[0] += cur_res3[0];
                sum_KH[1] += cur_res3[1];
                sum_KH[2] += cur_res3[2];

                R_Rot_Smooth(j_E[j].data(), points_for_field[i].data(), 
                        cell_j.GetCellCoords(), num_param, k, cur_res3);
                sum_RH[0] += cur_res3[0];
                sum_RH[1] += cur_res3[1];
                sum_RH[2] += cur_res3[2];
            }
        }
        field_E_SCH[i][0] = multE * sum_KE[0] - sum_RE[0];
        field_E_SCH[i][1] = multE * sum_KE[1] - sum_RE[1];
        field_E_SCH[i][2] = multE * sum_KE[2] - sum_RE[2];

        field_H_SCH[i][0] = multH * sum_KH[0] + sum_RH[0];
        field_H_SCH[i][1] = multH * sum_KH[1] + sum_RH[1];
        field_H_SCH[i][2] = multH * sum_KH[2] + sum_RH[2];
    }





    //====================Сохранение в формате .vtk====================

    res_fname = result_dir + "field_E_SCH.vtk";
    io_vtk.WritePntsWithField_C(res_fname, points_for_field, field_E_SCH);
    res_fname = result_dir + "field_H_SCH.vtk";
    io_vtk.WritePntsWithField_C(res_fname, points_for_field, field_H_SCH);



    //=========================Очистка памяти===================
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    return 0;
}