#ifndef _NODES_SURF_MESH_H_
#define _NODES_SURF_MESH_H_


#include <iostream>
#include <string>
#include <vector>
#include <array>

#include "cell.h"
#include "rwg_surf_mesh.h"


namespace bielcc {
    
//NOTE: по идее сработает и с треугольной сеткой

//===========================================================
//------------Quad surface mesh related to nodes-------------
//===========================================================
/**
 * @brief Quadrangular surface mesh class (dim = 3) related to
 * mesh nodes ("grey mesh"). Consist of quadrangular cells only.
 * @details Fields:
 *      _cell_listD - a list of quad cells for each node (CellQuad [num_nodes][nodes_num_frm[i]])
 *      _muD - list of total area values for each grey cell
*/
class NodeSurfMesh: public RwgSurfMesh<CellQuad> {
private:
    std::vector<std::vector<CellQuad>> _cell_listD;     // [num_nodes][nodes_num_frm[i]] - для каждой точки cell-соседи [4][3]
    std::vector<std::vector<double>> _muD_i;            // [num_nodes][nodes_num_frm[i]] - площадь каждой из составных частей D[i]
    std::vector<double> _muD;                           // [num_nodes] - total-площади серых ячеек

public:
    ~NodeSurfMesh() override = default;
    NodeSurfMesh(const NodeSurfMesh&) = default;
    NodeSurfMesh& operator=(const NodeSurfMesh&) = default;
    NodeSurfMesh(NodeSurfMesh&&) noexcept = default;
    NodeSurfMesh& operator=(NodeSurfMesh&&) noexcept = default;

    NodeSurfMesh(const RwgSurfMesh<CellQuad>& rwg_obj):
                    RwgSurfMesh<CellQuad>(rwg_obj) {this->FillMesh();};
    

    void FillMesh();

    const CellQuad& GetNodeCellD(int i, int j) const {return _cell_listD[i][j];};
    double GetMuNodeD(int i) const {return _muD[i];}
    double GetMuNodeDi(int i, int j) const {return _muD_i[i][j];}
};

}           // namespace bielcc
#endif      // _NODES_SURF_MESH_H