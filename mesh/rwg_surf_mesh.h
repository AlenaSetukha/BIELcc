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
#include <iomanip>

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
 *      _point_edges - edges numbers for each point (where point lies)
 * 
 *      _n_edges - number of unique edges
 *      _edges - unique edges (Segment [num_edges])
 *      _edge_num_cells - the number of cells that each edge borders (int [num_edges])
 *      _edge_cells - cell numbers for each edge (int [num_edges][*])
 *      _edge_mid - coordinates of the edge centers (double [num_edges][3])
 *      _edge_ort - coordinates of the edge orts (double [num_edges][3])
 * 
 *      _cell_edge_dirIndex - direction match between cell side and edge from the list  (int [num_frm][4 / 3])
 *      _cell_edges - number from the edge list for each side of the cell (int [num_frm][4 / 3])
 *      _cell_nghb_by_edge - side neighbors for each cell (int [num_frm][NVertex][edge_num_frm[k]])
 * 
 * Notes:
 *      _cell_edge_dirIndex[i][j] = {1 - same, -1 - opposite, 0 - empty edge}
 *      _cell_nghb_by_edge[i][j] = {-5} if the edge ia EMPTY or lies only in one OWN cell
 *      _edge_ort[j] - LOCAL direction, NO mult by dir_index
 * 
 */

 /**
  * @brief Detailed surface mesh class (dim = 3) for rwg-type.
  */
template<typename CellType>
class RwgSurfMesh: public SurfMesh<CellType> {
protected:
    int _n_points;
    std::vector<Point> _points;
    std::vector<int> _point_num_cells;
    std::vector<std::vector<int>> _point_cells;
    std::vector<std::array<int, SurfMesh<CellType>::_NVertex>> _cell_points;
    std::vector<std::vector<int>> _point_edges;

    
    int _n_edges;
    std::vector<Segment> _edges;
    std::vector<int> _edge_num_cells;
    std::vector<std::vector<int>> _edge_cells;
    std::vector<Point> _edge_mid, _edge_ort;
    std::vector<std::array<int, SurfMesh<CellType>::_NVertex>> _cell_edges;
    // Важное: в случае псевдочетырехугольной ячейки(2 вершины совпадают),
    //       _cell_edge_dir_indx[j][k] = 0, а _cell_edges[j][k] = -5 
    std::vector<std::array<int, SurfMesh<CellType>::_NVertex>> _cell_edge_dirIndex;
    std::vector<std::vector<std::vector<int>>> _cell_nghb_by_edge;

    void RwgSurfMeshFill();
    void FillUniquePoints();
    void FillUniqueEdges();
    void FillNeighbors();

public:
    RwgSurfMesh() = default;
    ~RwgSurfMesh() override = default;
    RwgSurfMesh(const RwgSurfMesh&) = default;
    RwgSurfMesh& operator=(const RwgSurfMesh&) = default;
    RwgSurfMesh(RwgSurfMesh&&) noexcept = default;
    RwgSurfMesh& operator=(RwgSurfMesh&&) noexcept = default;


    RwgSurfMesh(const SurfMesh<CellType>& mesh_obj): SurfMesh<CellType>(mesh_obj)
                                                       {this->RwgSurfMeshFill();};


    int GetNPoints() const {return _n_points;}
    int GetNEdges() const {return _n_edges;}





    const std::vector<Point>& GetPoints() const {return _points;};
    const Point& GetPoint(int i) const {return _points[i];}

    int GetPointNumCells(int i) const {return _point_num_cells[i];};

    int GetPointCells(int i, int j) const {return _point_cells[i][j];};

    const std::vector<std::array<int, SurfMesh<CellType>::_NVertex>>& GetCellPoints() const {return _cell_points;};
    const std::array<int, SurfMesh<CellType>::_NVertex>& GetCellPoints(int i) const {return _cell_points[i];}


    const std::vector<int>& GetPointEdges(int nPnt) const {return _point_edges[nPnt];}


    const std::vector<Segment>& GetEdges() const {return _edges;}
    const Segment& GetEdge(int i) const {return _edges[i];}

    int GetEdgeNumCells(int i) const {return _edge_num_cells[i];}

    const std::vector<std::vector<int>>& GetEdgeCells() const {return _edge_cells;}
    int GetEdgeCell(int i, int j) const {return _edge_cells[i][j];}

    const std::vector<Point>& GetEdgeMid() const {return _edge_mid;}
    const Point& GetEdgeMid(int i) const {return _edge_mid[i];}

    const std::vector<Point>& GetEdgeOrt() const {return _edge_ort;}
    const Point& GetEdgeOrt(int i) const {return _edge_ort[i];} 

    const std::array<int, SurfMesh<CellType>::_NVertex>& GetCellEdges(int i) const {return _cell_edges[i];}
    int GetCellEdge(int i, int j) const {return _cell_edges[i][j];}


    int GetCellEdgeDirIndx(int i, int j) const {return _cell_edge_dirIndex[i][j];}




    double GetPointsDist(int pntA, int pntB) const {
        return dist(_points[pntA].GetCoord(), _points[pntB].GetCoord());
    }

    double GetEdgeLen(int i) const {return GetPointsDist(_edges[i][0], _edges[i][1]);}



    /**
        * @brief Export surf mesh to vtk-file format.
        * @param filename vtk-file
    */
    void ExportRwgfMeshToVTK(const std::string& filename);
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
        for (int i = 0; i < this->_NVertex; ++i) {
            _points.emplace_back(Point(cell.GetVertex(i)));
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
    _cell_points.reserve(this->_num_frm);

    for (const auto& cell: this->_cell_list) {
        std::array<int, this->_NVertex> vertex;
        for (int i = 0; i < this->_NVertex; ++i) {
            auto it = std::find(_points.begin(), _points.end(), Point(cell.GetVertex(i)));
            int pnt_N = std::distance(_points.begin(), it);
            vertex[i] = pnt_N;


            bool firstEntrance = true;
            for (int j = 0; j < i; j++) {
                if (vertex[j] == pnt_N) {
                    firstEntrance = false;
                    break;
                }
            }

            if (firstEntrance) { 
                _point_num_cells[pnt_N]++;
            }
        }
        _cell_points.push_back(vertex);
    }


    //=================================
    //------Forming point frames-------
    //=================================
    for (int i = 0; i < _n_points; i++) {
        std::vector<int> point_frames;
        for (int j = 0; j < (int)_cell_points.size(); j++) {
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
    for (int i = 0; i < (int)this->_cell_list.size(); i++) {
        for (int j = 0; j < this->_NVertex; j++) {
            int NA = _cell_points[i][j];
            int Next = (j == (this->_NVertex - 1)) ? 0 : j + 1;
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

    for (int i = 0; i < (int)this->_cell_list.size(); i++) {
        std::array<int, this->_NVertex> edgesNums;
        for (int j = 0; j < this->_NVertex; j++) {
            int NA = _cell_points[i][j];
            int Next = (j == (this->_NVertex - 1)) ? 0 : j + 1;
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
    //---Forming edge mid / ort--------
    //=================================
    _edge_mid.resize(_n_edges);
    _edge_ort.resize(_n_edges);
    _point_edges.resize(_n_points);
    for (int i = 0; i < _n_edges; i++) {
        int nA = _edges[i][0], nB = _edges[i][1];
        _point_edges[nA].push_back(i);
        _point_edges[nB].push_back(i);
        double A[3], B[3];
        for (int g = 0; g < 3; g++) {
            A[g] = _points[nA][g];
            B[g] = _points[nB][g];
        }

        _edge_mid[i] = (Point(A) + Point(B)) * 0.5;
        double dl = GetPointsDist(nA, nB);
        _edge_ort[i] = (Point(B) - Point(A)) * (1. / dl);
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
        _cell_nghb_by_edge[i].resize(this->_NVertex);
        for (int j = 0; j < this->_NVertex; j++) {
            int edgeN = _cell_edges[i][j];
            if (edgeN >= 0) {
                std::vector<int> cells_numbers = _edge_cells[edgeN];
                auto it = std::find(cells_numbers.begin(), cells_numbers.end(), i);
                if (it != cells_numbers.end()) {
                    cells_numbers.erase(it);
                    if (!cells_numbers.empty()) {
                        _cell_nghb_by_edge[i][j] = cells_numbers;
                    } else {
                        // Ребро только в одной (собственной) ячейке
                        _cell_nghb_by_edge[i][j] = {-5};
                    }
                } else {
                    throw std::runtime_error("Filling RWG Mesh error: FillNeighbors");
                }
            } else {
                // Пустое ребро
                _cell_nghb_by_edge[i][j] = {-5};
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








//===========================================================
//----------------Export to vtk-file-------------------------
//===========================================================
template<typename CellType>
void RwgSurfMesh<CellType>::ExportRwgfMeshToVTK(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    file << std::scientific << std::setprecision(16);
    
    // Header
    file << "# vtk DataFile Version 3.0\n";
    file << "RwgSurfMesh export\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n\n";
    

    // Points
    file << "POINTS " << _n_points << " double\n";
    for (int i = 0; i < _n_points; ++i) {
        file << _points[i][0] << " " 
             << _points[i][1] << " " 
             << _points[i][2] << "\n";
    }
    file << "\n";
    

    // Calculate total size for CELLS section
    // Треугольник: `3 0 1 2` - это 4 числа (1 для счетчика + 3 индекса)
    // Четырехугольник: `4 0 1 2 3` - это 5 чисел (1 для счетчика + 4 индекса)
    constexpr int nVertex = this->_NVertex;
    int totalSize = this->_num_frm * (nVertex + 1);
    file << "CELLS " << this->_num_frm << " " << totalSize << "\n";
    for (int i = 0; i < this->_num_frm; ++i) {
        file << nVertex;
        for (int j = 0; j < nVertex; ++j) {
            file << " " << _cell_points[i][j];
        }
        file << "\n";
    }
    file << "\n";
    



    // Cell types
    file << "CELL_TYPES " << this->_num_frm << "\n";
    int vtkCellType = (nVertex == 3) ? 5 : 9;  // 5=VTK_TRIANGLE, 9=VTK_QUAD
    for (int i = 0; i < this->_num_frm; ++i) {
        file << vtkCellType << "\n";
    }
    file << "\n";
    




    // Cell data
    file << "CELL_DATA " << this->_num_frm << "\n";
        
    // Cell edges
    file << "SCALARS cell_id int 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < this->_num_frm; ++i) {
        file << i << "\n";
    }
    file << "\n";
        
    // Edge numbers for each cell
    file << "FIELD CellFields 3\n";
        
    // Cell edges
    file << "edge_indices " << nVertex << " " << this->_num_frm << " int\n";
    for (int i = 0; i < this->_num_frm; ++i) {
        for (int j = 0; j < nVertex; ++j) {
            file << _cell_edges[i][j];
            if (j < nVertex - 1) file << " ";
        }
        file << "\n";
    }
    file << "\n";
        
    // Edge direction indices
    file << "edge_dir_index " << nVertex << " " << this->_num_frm << " int\n";
    for (int i = 0; i < this->_num_frm; ++i) {
        for (int j = 0; j < nVertex; ++j) {
            file << _cell_edge_dirIndex[i][j];
            if (j < nVertex - 1) file << " ";
        }
        file << "\n";
    }
    file << "\n";
        
    // Cell normals (if available - need to compute or access from parent)
    file << "cell_normal 3 " << this->_num_frm << " double\n";
    for (int i = 0; i < this->_num_frm; ++i) {
        const double* normal = this->_cell_list[i].GetNorm();
        file << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
    }
    file << "\n";
    



    // Point data
    file << "POINT_DATA " << _n_points << "\n";
        
    // Number of cells per point
    file << "SCALARS point_num_cells int 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < _n_points; ++i) {
        file << _point_num_cells[i] << "\n";
    }
    file << "\n";
        
    // Number of edges per point
    file << "SCALARS point_num_edges int 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < _n_points; ++i) {
        file << _point_edges[i].size() << "\n";
    }
    file << "\n";
    

    file.close();
}



}           // namespace bielcc
#endif      // _RWG_SURF_MESH_H