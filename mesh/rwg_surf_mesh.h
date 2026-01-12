#ifndef _RWG_SURF_MESH_H_
#define _RWG_SURF_MESH_H_


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <type_traits>

#include "surf_mesh.h"
#include "point.h"
#include "segment.h"





/// ЧИСТКА ПАМЯТИ (ПРАВИЛЬНАЯ):
// std::for_each( orders.begin(), orders.end(), std::default_delete<Order>() );
// edge-ort безотносительно реальному напрвлению!!!!!
// Псевдочетырехугольные ячейки
// Проверить - что быстрее вложенный  vector или ***


namespace bielcc {

//===========================================================
//---------------RWG-type surface mesh class-----------------
//===========================================================
/**
 * @brief Detailed surface mesh class (dim = 3) for rwg-type.
 * Consist of quadrangular or triangular cells.
 * Fields:
 *      _n_points - number of unique grid points
 *      _points - list of unique points coordinates (double [num_points][3])
 *      _point_num_cells - number of cells in which each point lies (int [num_points])
 *      _point_cells - cell numbers for each point (int [num_points][*])
 *      _cell_points - corner points numbers for each cell (int [num_frm][4 / 3]) 
 * 
 *      _n_edges - number of unique edges
 *      _edges - unique edges (Segment [num_edges])
 *      _edge_num_cells - the number of cells that each edge borders (int [num_edges])
 *      _edge_cells - cell numbers for each edge (int [num_edges][*])
 *      _edge_mid - coordinates of the edge centers (double [num_edges][3])
 *      _edge_ort - coordinates of the edge orts (double [num_edges][3])
 * 
 *      _cell_edge_dirIndex - match the direction of the cell side with the edge
 *                  from the list (1 - same, -1 - opposite, 0 - empty edge) (int [num_frm][4 / 3])
 *      _cell_edges - number from the edge list for each side of the cell (int [num_frm][4 / 3])
 *      _cell_nghb_by_edge - side neighbors for each cell (int [num_frm][NVertex][edge_num_frm[k]])
 */

 /**
  * @brief Detailed surface mesh class (dim = 3) for rwg-type.
  */
template<typename CellType>
class RwgSurfMesh: public SurfMesh<CellType> {
private:
    static constexpr std::size_t _NVertex = (std::is_same<CellType, CellQuad>::value) ? 4 : 3;

    int _n_points;
    std::vector<Point> _points;
    std::vector<int> _point_num_cells;
    std::vector<std::vector<int>> _point_cells;
    std::vector<std::array<int, _NVertex>> _cell_points;

    
    int _n_edges;
    std::vector<Segment> _edges;
    std::vector<int> _edge_num_cells;
    std::vector<std::vector<int>> _edge_cells;
    std::vector<Point> _edge_mid, _edge_ort;
    std::vector<std::array<int, _NVertex>> _cell_edges;
    // Важное: в случае псевдочетырехугольной ячейки(2 вершины совпадают),
    //       edge_dir_indx[j][k] = 0, а frm_edges[j][k] = -5 
    std::vector<std::array<int, _NVertex>> _cell_edge_dirIndex;
    std::vector<std::vector<std::vector<int>>> _cell_nghb_by_edge;

    void RwgSurfMeshFill();
    void FillUniquePoints();
    void FillUniqueEdges();
    void FillNeighbors();

public:
    RwgSurfMesh(const RwgSurfMesh<CellType>& rwg_obj): _n_points(rwg_obj._n_points),
                        _points(rwg_obj._points), _point_num_cells(rwg_obj._point_num_cells),
                        _point_cells(rwg_obj._point_cells), _cell_points(rwg_obj._cell_points),
                        _n_edges(rwg_obj._n_edges), _edges(rwg_obj._edges),
                        _edge_num_cells(rwg_obj._edge_num_cells),
                        _edge_cells(rwg_obj._edge_cells),
                        _edge_mid(rwg_obj._edge_mid), _edge_ort(rwg_obj._edge_ort),
                        _cell_edges(rwg_obj._cell_edges),
                        _cell_edge_dirIndex(rwg_obj._cell_edge_dirIndex), 
                        _cell_nghb_by_edge(rwg_obj._cell_nghb_by_edge) {}

    RwgSurfMesh(const SurfMesh<CellType>& mesh_obj): SurfMesh<CellType>(mesh_obj)
                                                        {this->RwgSurfMeshFill();};
    ~RwgSurfMesh() = default;


    int GetNPoints() const {return _n_points;}
    int GetNEdges() const {return _n_edges;}




    const std::vector<Point>& GetPoints() const {return _points;};
    const std::vector<int>& GetPointNumCells() const {return _point_num_cells;};
    const std::vector<std::vector<int>>& GetPointCell() const {return _point_cells;};
    const std::vector<std::array<int, _NVertex>>& GetCellPoints() const {return _cell_points;};

    const std::vector<Segment>& GetEdges() const {return _edges;};
    const std::vector<int>& GetEdgeNumCells() const {return _edge_num_cells;};
    const std::vector<std::vector<int>>& GetEdgeCells() const {return _edge_cells;};
    const std::vector<Point>& GetEdgeMid() const {return _edge_mid;};
    const std::vector<Point>& GetEdgeOrt() const {return _edge_ort;};
    const std::vector<std::array<int, _NVertex>>& GetCellEdges() const {return _cell_edges;};
};










//===========================================================
//-----------------Fill unique points------------------------
//===========================================================
template<typename CellType>
void RwgSurfMesh<CellType>::FillUniquePoints()
{
    //=================================
    //------Forming unique points======
    //=================================
    for (const auto& cell: this->_cell_list) {
        for (int i = 0; i < _NVertex; ++i) {
            _points.push_back(Point(cell.GetVertex(i)));
        }
    }


    std::sort(_points.begin(), _points.end(), [](const Point &a, const Point &b) {
        return a < b;
    });

    _points.resize(std::unique(_points.begin(), _points.end()) - _points.begin());
    _n_points = _points.size();

    //=================================
    //------Forming cell points--------
    //=================================
    _point_num_cells.resize(_n_points);
    for (const auto& cell: this->_cell_list) {
        std::array<int, _NVertex> vertex;
        for (int i = 0; i < _NVertex; ++i) {
            auto it = std::find(_points.begin(), _points.end(), Point(cell.GetVertex(i)));
            int pnt_N = std::distance(_points.begin(), it);
            vertex[i] = pnt_N;
            _point_num_cells[pnt_N]++;
        }
        _cell_points.push_back(vertex);
    }


    //=================================
    //------Forming point frames-------
    //=================================
    for (int i = 0; i < _n_points; i++) {
        std::vector<int> point_frames;
        for (int j = 0; j < _cell_points.size(); j++) {
            auto it = std::find(_cell_points[j].begin(), _cell_points[j].end(), i);
            if (it != _cell_points[j].end()) {
               point_frames.push_back(j);
            }
        }
        _point_cells.push_back(point_frames);
    }
}









//===========================================================
//-------------------Fill unique edges-----------------------
//===========================================================
template<typename CellType>
void RwgSurfMesh<CellType>::FillUniqueEdges()
{
    //=================================
    //------Forming unique Edges======
    //=================================
    for (int i = 0; i < this->_cell_list.size(); i++) {
        for (int j = 0; j < _NVertex; j++) {
            int NA = _cell_points[i][j];
            int Next = (j == (_NVertex - 1)) ? 0 : j + 1;
            int NB = _cell_points[i][Next];
            if (NA < NB) {
                _edges.push_back(Segment(NA, NB));
            } else if (NA > NB) {
                _edges.push_back(Segment(NB, NA));
            }
        }
    }

    std::sort(_edges.begin(), _edges.end(), [](const Segment &a, const Segment &b) {
        return a < b;
    });

    _edges.resize(std::unique(_edges.begin(), _edges.end()) - _edges.begin());
    _n_edges = _edges.size();


    //=================================
    //------Forming edge frames--------
    //=================================
    _edge_num_cells.resize(_n_edges);
    _edge_cells.resize(_n_edges);
    _cell_edge_dirIndex.resize(this->_num_frm);

    for (int i = 0; i < this->_cell_list.size(); i++) {
        std::array<int, _NVertex> edgesNums;
        for (int j = 0; j < _NVertex; j++) {
            int NA = _cell_points[i][j];
            int Next = (j == (_NVertex - 1)) ? 0 : j + 1;
            int NB = _cell_points[i][Next];
            if (NA == NB) {
                // Pseudo-edge
                edgesNums[j] = -5;
                _cell_edge_dirIndex[i][j] = 0;
            } else {
                int SegN;
                if (NA < NB) {
                    auto it = std::find(_edges.begin(), _edges.end(), Segment(NA, NB));
                    SegN = std::distance(_edges.begin(), it);
                    _cell_edge_dirIndex[i][j] = 1;
                } else {
                    auto it = std::find(_edges.begin(), _edges.end(), Segment(NB, NA));
                    SegN = std::distance(_edges.begin(), it);
                    _cell_edge_dirIndex[i][j] = -1;
                }
                edgesNums[j] = SegN;
                _edge_num_cells[SegN]++;
                _edge_cells[SegN].push_back(i);
            }
        }
        _cell_edges.push_back(edgesNums);
    }




    //=================================
    //---Forming edge coordinates------
    //=================================
    _edge_mid.resize(_n_edges);
    _edge_ort.resize(_n_edges);
    for (int i = 0; i < _n_edges; i++) {
        int nA = _edges[i][0], nB = _edges[i][1];
        double A[3], B[3];
        _points[nA].GetCoord(A);
        _points[nB].GetCoord(B);
        _edges[i].SetCoords(A, B);

        _edge_mid[i] = (Point(A) + Point(B)) * 0.5;
        _edge_ort[i] = Point(B) - Point(A);
    }
}





//===========================================================
//-------------------Fill neighbors-----------------------
//===========================================================
template<typename CellType>
void RwgSurfMesh<CellType>::FillNeighbors()
{
    _cell_nghb_by_edge.resize(this->_num_frm);
    for (int i = 0; i < this->_num_frm; i++) {
        _cell_nghb_by_edge[i].resize(_NVertex);
        for (int j = 0; j < _NVertex; j++) {
            int edgeN = _cell_edges[i][j];
            if (edgeN > 0) {
                std::vector<int> cells_numbers = _edge_cells[edgeN];
                auto it = std::find(cells_numbers.begin(), cells_numbers.end(), i);
                if (it != cells_numbers.end()) {
                    cells_numbers.erase(it);
                } else {
                    throw std::runtime_error("Filling RWG Mesh error: FillNeighbors");
                }
                _cell_nghb_by_edge[i][j] = cells_numbers;
            }
        }
    }
}







//===========================================================
//-----------------Filling RWG-Type mesh---------------------
//===========================================================
template<typename CellType>
void RwgSurfMesh<CellType>::RwgSurfMeshFill()
{
    FillUniquePoints();
    FillUniqueEdges();
    FillNeighbors();
}


}           // namespace bielcc
#endif      // _RWG_SURF_MESH_H