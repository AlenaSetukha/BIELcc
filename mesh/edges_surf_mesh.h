#ifndef _EDGES_SURF_MESH_H_
#define _EDGES_SURF_MESH_H_


#include <iostream>
#include <string>
#include <vector>
#include <array>

#include "cell.h"
#include "rwg_surf_mesh.h"


namespace bielcc {
//NOTE: по идее сработает и с треугольной сеткой
 
// TODO

//===========================================================
//------------Quad surface mesh related to edges-------------
//===========================================================
/**
 * @brief Quadrangular surface mesh class (dim = 3) related to
 * mesh edges ("blue mesh"). Consist of quadrangular cells only.
 * @details Fields:
 *      _cell_listG - a list of quad cells for each edge (CellQuad [num_edges][edge_num_frm[i]])
 *      _muG - list of total area values for each blue cell
*/
class EdgeSurfMesh: public RwgSurfMesh<CellQuad> {
private:
    std::vector<std::vector<CellQuad>> _cell_listG;     // [num_edges][edge_num_frm[i]] - для каждого ребра cell-соседи [4][3]
    std::vector<std::vector<double>> _diag_G;           // [num_edge][edge_num_frm[i]] - диаметр каждой из составных частей G[i]
    std::vector<std::vector<double>> _muG_i;            // [num_edges][edge_num_frm[i]] - площадь каждой из составных частей G[i]
    std::vector<double> _muG;                           // [num_edges] - total-площади голубых ячеек

public:
    ~EdgeSurfMesh() override = default;
    EdgeSurfMesh(const EdgeSurfMesh&) = default;
    EdgeSurfMesh& operator=(const EdgeSurfMesh&) = default;
    EdgeSurfMesh(EdgeSurfMesh&&) noexcept = default;
    EdgeSurfMesh& operator=(EdgeSurfMesh&&) noexcept = default;



    EdgeSurfMesh(const RwgSurfMesh<CellQuad>& rwg_obj):
                    RwgSurfMesh<CellQuad>(rwg_obj) {this->FillMesh();};

    void FillMesh();

    const CellQuad& GetEdgeCellG(int i, int j) const {return _cell_listG[i][j];};
    double GetMuEdgeG(int i) const {return _muG[i];}
    double GetMuEdgeGi(int i, int j) const {return _muG_i[i][j];}
    double GetDiamGPart(int i, int j) const {return _diag_G[i][j];};
};

}           // namespace bielcc
#endif      // _EDGES_SURF_MESH_H