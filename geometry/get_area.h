#ifndef _GET_AREA_H_
#define _GET_AREA_H_

#include <iostream>
#include <vector>
#include <array>

#include "element_geom.h"
#include "surf_mesh.h"

namespace bielcc {
//===========================================================
//------Area where the point is located (in / out)-----------
//===========================================================
/**
    * @brief Location of a point relative to a given closed body
    * @param x point for determining 
    * @param cells cell list of the body
    * @param norms normal vector list
    * @result 0 - inside the object, -1 - outside.
*/
template<size_t CellPoints>
int get_area(const double* x,
    const std::vector<std::array<std::array<double, 3>, CellPoints>> &cells,
                            const std::vector<std::array<double, 3>> &norms)
{
    double teta, vec[3]{};
    double sum = 0.;

    int num_frm = cells.size();
    for (int i = 0; i < num_frm; i++) {
        vec[0] = cells[i][0][0] - x[0];
        vec[1] = cells[i][0][1] - x[1];
        vec[2] = cells[i][0][2] - x[2];

        //123
        teta = solid_angle(cells[i][0].data(), cells[i][1].data(), cells[i][2].data(), x);

        if (scal_prod(vec, norms[i].data()) > 0) {
            sum -= fabs(teta);
        } else {
            sum += fabs(teta);
        }

        //134
        if (CellPoints == 4) {
            teta = solid_angle(cells[i][0].data(), cells[i][2].data(), cells[i][3].data(), x);

            if (scal_prod(vec, norms[i].data()) > 0) {
                sum -= fabs(teta);
            } else {
                sum += fabs(teta);
            }
        }

    }

    //Result
    if (sum > -0.5) {
        return 0;
    } else return -1;
}







/**
    * @brief Location of a point relative to a given closed body
    * @param x point for determining 
    * @param surf_mesh object quadrangular surface mesh
    * @result 0 - inside the object, -1 - outside.
*/
template <typename CellType>
int get_area(const double* x, const SurfMesh<CellType>& surf_mesh)
{

    double teta, vec[3]{};
    double sum = 0.;

    int num_frm = surf_mesh.GetNFrm();
    for (int i = 0; i < num_frm; i++) {
        CellType cell_i = surf_mesh.GetCell(i);
        const double* v0 = cell_i.GetVertex(0);
        const double* v1 = cell_i.GetVertex(1);
        const double* v2 = cell_i.GetVertex(2);

        vec[0] = v0[0] - x[0];
        vec[1] = v0[1] - x[1];
        vec[2] = v0[2] - x[2];

        //123
        teta = solid_angle(v0, v1, v2, x);

        if (scal_prod(vec, cell_i.GetNorm()) > 0) {
            sum -= fabs(teta);
        } else {
            sum += fabs(teta);
        }

        //134
        if (std::is_same<CellType, CellQuad>::value) {
            const double* v3 = cell_i.GetVertex(3);
            teta = solid_angle(v0, v2, v3, x);

            if (scal_prod(vec, cell_i.GetNorm()) > 0) {
                sum -= fabs(teta);
            } else {
                sum += fabs(teta);
            }
        }

    }

    //Result
    if (sum > -0.5) {
        return 0;
    } else return -1;
}



}       // namespace bielcc
#endif  // _GET_AREA_H_