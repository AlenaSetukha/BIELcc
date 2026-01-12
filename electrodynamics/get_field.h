#ifndef _GET_FIELD_H_
#define _GET_FIELD_H_

#include <iostream>
#include <complex>
#include <vector>
#include <array>
#include <omp.h>

#include "num_param.h"
#include "kernel_param.h"
#include "K.h"


namespace bielcc {

//=======================================================
//----------Ideal conductor Electric field---------------
//=======================================================
/**
    * @brief Electric field around ideal condictor. Random points.
    * @details Calculation of the total Electric field based
    * only on electric currents:
    *      E(x) = mult * sum_j(K[sigma[j]], j_vec[j](x)) + Einc(x)
    * with K = K_Far. 
    * @param surf_mesh object surface mesh
    * @param j_E surface electric currents
    * @param mult multiplier for operator K[]
    * @param EInc_param Einc function parameters (k, k_vec, e0...)
    * @param EInc_Func Einc function
    * @param num_param numerical parameters for integration
    * @param points points for field calculation
    * @param field_E result
*/ 
template<typename CellType, typename KType, typename P>
void GetFieldE_Ideal(const SurfMesh<CellType>& surf_mesh,
         const std::vector<std::array<std::complex<double>, 3>> &j_E,
                                           std::complex<double> mult,
                                const KernelParam<KType>& EInc_param,
      void (*EInc_Func)(const double*, const KernelParam<KType>&, P*),
                                            const NumParam& num_param,
                     const std::vector<std::array<double, 3>> &points,
            std::vector<std::array<std::complex<double>, 3>> &field_E)
{
    int n_points = points.size();
    KType k = EInc_param.k;
    int num_frm = surf_mesh.GetNFrm();

    #pragma omp parallel for
    for (int i = 0; i < n_points; i++) {
        double y[3]{};
        field_E[i][0] = 0., field_E[i][1] = 0., field_E[i][2] = 0.;

        std::complex<double> tmp[3]{};
        for (int j = 0; j < num_frm; j++) {
            CellType cell_j = surf_mesh.GetCell(j);
            get_center_mass(cell_j.GetCellCoords(), y);
            if (dist(points[i].data(), y) < num_param.analytic_calc * get_diam(cell_j.GetCellCoords())) {
                // Near field zone
                K_RotRot_Near_Smooth(j_E[j].data(), points[i].data(), cell_j.GetCellCoords(), num_param, k, tmp);
            } else {
                // Far field zone
                K_RotRot_Far(j_E[j].data(), points[i].data(), cell_j.GetCellCoords(), num_param, k, tmp);
            }
            for (int g = 0; g < 3; g++) {
                field_E[i][g] += tmp[g];
            }
        }
        field_E[i][0] *= mult, field_E[i][1] *= mult, field_E[i][2] *= mult;

        // Outer/ Inner volume
        std::complex<double> e_inc[3]{};
        EInc_Func(points[i].data(), EInc_param, e_inc);
        field_E[i][0] += e_inc[0];
        field_E[i][1] += e_inc[1];
        field_E[i][2] += e_inc[2];
    }
}

}           // namespace bielcc
#endif      // _GET_FIELD_H_