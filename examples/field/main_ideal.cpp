#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>

#include <omp.h>
//#include <cblas.h>
//#include <lapacke.h>
//#include <openblas_config.h>



#include "material.h"
#include "radiation_source.h"
#include "surf_mesh.h"
#include "io_aereco.h"
#include "io_vtk.h"
#include "num_param.h"
#include "incident_field.h"
#include "get_field.h"
#include "kernel_param.h"

/**
    * Пример расчета электрического поля вокруг идеального
    * проводника в точках некоторой плоскости.
    * 
    * Расчет поля с помощью оператора K[].
    * E(x) = [i / (omega * eps_ext)] * sum_j K[sigma_j, j_j](x) + E_inc(x)
    * 
    * 
    * Геометрия: 1 объект в формате .dat / .vtk
    * 
*/ 


using namespace bielcc;

int main(int argc, char **argv)
{
    const std::string geom_fname = argc > 1 ? argv[1] :
                                    "../examples/data/grids/plate/geodat_20_20.dat", \
                      j_fname = argc > 2 ? argv[2] :
                                    "../examples/data/currents/E", \
                      result_dir = argc > 3 ? argv[3] :
                                    "../examples/results_ideal/", \
                      pntsForField_fname = argc > 4 ? argv[4] :
                                    "../examples/data/grids/plate/grid.gr";
    auto start = std::chrono::high_resolution_clock::now();

    omp_set_dynamic(0);
    omp_set_num_threads(18);

    //===============1. Mesh object creation========================
    bielcc::SurfMesh<CellQuad> surf_mesh(geom_fname);
    std::cout << "Grid step: " << surf_mesh.GetGridStep() << std::endl;
    std::cout << "Surface square: " << surf_mesh.GetSurfSquare(0) << std::endl;


    


    //===============2. Forming EM parameters for one frequency=====
    double freq = 477134652.095072; // in Hz for k = 10
    double eps_r = 1.;
    double mu_r = 1.;


    MaterialParam vacuum_material = MaterialParam("Vacuum",
                          ChangeMaterialTypeDir::Isotropic,
                     ChangeMaterialTypeFreq::NonDispersive,
                                  eps_r, mu_r, 0., 0., 0.);
    vacuum_material.Calc(freq);


    double eps = std::get<double>(vacuum_material.GetEps());
    double mu = std::get<double>(vacuum_material.GetMu());
    std::cout << "Eps: " << eps << ", Mu: " << mu << std::endl;


    double k_ort[3] = {-1., 0., 0.};
    double e0_ort[3] = {0., -1., 0.};

    RadiationSource src = RadiationSource(RadiationSourceType::PlaneWave,
                                                           k_ort, e0_ort);

    src.CalcK(freq, vacuum_material); // k во внешней среде
    double k = src.GetK();
    double k_vec[3];
    src.GetKVec(k_vec);
    std::cout << "Волновое число в вакууме: " << k << std::endl;
    std::cout << "Волновой вектор в вакууме: " << k_vec[0] << " " <<
                            k_vec[1] << " " << k_vec[2] << std::endl;
    std::cout << "Omega: " << 2. * M_PI * freq << std::endl;



    //==============================================
    //--------Global numerical parameters-----------
    //==============================================
    const double IntegralAccuracy = 0.000001;
    const double SmoothingR_r = 0.2;                   // relative to small cell step h2
    const double SmoothingRSeg_r = 0.2;                // relative to small cell step h2
    const double AnalyticCalcR = 3.0;                  // relative to grid step h
    const int NCellStart = 20;                       
    const int NSegStart = 40;
    const int PMax = 1;                                // 2^{PMax} steps in surface integration
    const int PMaxSeg = 1;                             // 2^{PMaxSeg} steps in segment integration
    
    NumParam num_param(IntegralAccuracy, SmoothingR_r, SmoothingRSeg_r,
                    AnalyticCalcR, NCellStart, NSegStart, PMax, PMaxSeg);




    // External field parameters
    KernelParam<double> EIncParam;
    EIncParam.k = k;
    for (int i = 0; i < 3; i++) {
        EIncParam.vec_dbl[i] = k_vec[i];
        EIncParam.e0[i] = e0_ort[i];
    }





    //===============3. Reading surface currents====================
    IO_AERECO io_aer;
    std::vector<std::array<std::complex<double>, 3>> j_E =
                    io_aer.ReadVecs_C(j_fname + "_real.gv", j_fname + "_image.gv");





    //===============4. Points for field calculation================
    std::vector<std::array<double, 3>> points_for_field = 
                            io_aer.ReadVecs_D(pntsForField_fname);
    



    //===============5. Field in points=============================
    int num_points = points_for_field.size();
    std::vector<std::array<std::complex<double>, 3>> field_E_SCH(num_points);

    // Multiplier for K[] = i / (omega * eps_ext_full)
    std::complex<double> mult = std::complex<double>(0., 1.) /
                                      (2. * M_PI * freq * eps);
    GetFieldE_Ideal(surf_mesh, j_E, mult, EIncParam, PlaneWaveFieldE<double>,
                                    num_param, points_for_field, field_E_SCH);


    //===============6. Record results to file in .vtk format======
    IO_VTK io_vtk;
    // Writing collocation points in vtk file
    std::vector<std::array<double, 3>> rkt = surf_mesh.FormCollocPoints();
    std::string colloc_pnt_fname = result_dir + "body.vtk";
    io_vtk.WriteVecs(colloc_pnt_fname, rkt);


    // Writing electric field in points in vtk file
    std::string fieldE_fname = result_dir + "field/E.vtk";
    io_vtk.WritePntsWithField_C(fieldE_fname, points_for_field, field_E_SCH);








    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    return 0;
}