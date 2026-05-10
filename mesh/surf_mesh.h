#ifndef _SURF_MESH_H_
#define _SURF_MESH_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <utility>
#include <sstream>
#include <map>
#include <algorithm>
#include <type_traits>

#include "cell.h"
#include "element_geom.h"

//===========================================================
//--------------------Surface mesh class---------------------
//===========================================================
/**
 * Surface mesh class (dim = 3).
 * Consist of quadrangular or triangular cells.
 * Fields:
 *      cell_list[num_frm] - cell list of all objects
 *      num_obj / num_mod / num_frm - total number of objects/modules/cells
 *      num_bl / num_blSeg - total number of tear lines / segments
 *      grid_step - step of the mesh
 *      _surf_square - total surface area of ​​each object
 *      _obj_num_frm - number of cells on each object
 *      x_bound[num_blSeg][2][3] - coordinates of tear segments
 *      beg_endObjFrm[obj][2] - first / last cell on object
 *      beg_endMod[mod][2] - first / last cell in module within one object
 *      beg_endObj[obj][2] - first / last module in object
 * 
 * Note: The order of the cells corresponds to the order of the cells in the file
 */


namespace bielcc {
/**
    * @brief Basic surface mesh class (dim = 3) consisting of quadrangular or triangular cells.
    * @note Cells order corresponds file cells order
*/
template<typename CellType>
class SurfMesh {
protected:
    static constexpr int _NVertex = (std::is_same<CellType, CellQuad>::value) ? 4 : 3;
    std::vector<CellType> _cell_list;

    int _num_obj, _num_mod, _num_frm;
    int _num_bl, _num_blSeg;
    double _grid_step;
    std::vector<double> _surf_square;

    std::vector<int> _obj_num_frm;
    std::vector<std::pair<int, int>> _beg_endMod, _beg_endObj, _beg_endObjFrm;
    std::vector<std::array<std::array<double, 3>, 2>> _x_bound;

    void check_step(const CellType& cell, double cur_grid_step);
    void ReadFromVtkFile(const std::string& vtkFilename);  // пока только один общий объект  
    void ReadFromDatFile(const std::string& datFilename);
    void FillMeshFromVTK(const std::vector<std::vector<int>>& polygons,
                      const std::vector<std::array<double, 3>>& points);

public:
    SurfMesh() = default;
    virtual ~SurfMesh() = default;
    SurfMesh(const SurfMesh&) = default;
    SurfMesh& operator=(const SurfMesh&) = default;
    // Move-constructors
    SurfMesh(SurfMesh&&) noexcept = default;
    SurfMesh& operator=(SurfMesh&&) noexcept = default;


    void ReadFromFile(const std::string& filename);

    /**
        * @brief Reading surf mesh from file
        * @param filename Full path to surface mesh file
        * @note Supported: .dat, .vtk
    */ 
    SurfMesh(const std::string& filename) { this->ReadFromFile(filename);};




    /**
        * @brief Number of all objects
    */
    int GetNObj() const {return _num_obj;}
    /**
        * @brief Number of all modules
    */
    int GetNMod() const {return _num_mod;}
    /**
        * @brief Number of all cells
    */
    int GetNFrm() const {return _num_frm;}
    /**
        * @brief Number of all tear lines
    */
    int GetNBl() const {return _num_bl;}
    /**
        * @brief Number of all tear segments
    */
    int GetNBlSeg() const {return _num_blSeg;}
    /**
        * @brief Total grid step (common)
    */
    double GetGridStep() const {return _grid_step;}
    /**
        * @brief Surface area of ​​an object
    */
    double GetSurfSquare(int objN) const {return _surf_square[objN];};






    /**
        * @brief i-th cell from cell list
        * @return Const copy of CellType cell
    */ 
    const CellType& GetCell(int i) const {return _cell_list[i];}

    /**
        * @brief List of all cells forming surf-mesh
        * @return Const copy of CellType cell list
    */ 
    const std::vector<CellType>& GetTotalCellList() const {return _cell_list;}

    
    /**
        * @brief Common list of pairs <first, last> cell numbers in each object
        * @return Const copy of pairs vector
    */ 
    const std::vector<std::pair<int, int>>& GetBegEndObj() const {return _beg_endObj;}

    /**
        * @brief First and last module in object objN
        * @param objN object number
    */
    const std::pair<int, int>& GetFLModInObj(int objN) const {return _beg_endObj[objN];}



    /**
        * @brief Common list of pairs <first, last> cell numbers in each module
        * @return Const copy of pairs vector
    */ 
    const std::vector<std::pair<int, int>>& GetBegEndMod() const {return _beg_endMod;}
    /**
        * @brief First and last cell in module modN
        * @param objN object number
    */
    const std::pair<int, int>& GetFLFrmInMod(int modN) const {return _beg_endMod[modN];}


    /**
        * @brief First and last cell in object objN
        * @param objN object number
    */
    const std::pair<int, int>& GetFLFrmInObj(int objN) const {return _beg_endObjFrm[objN];}



    /**
        * @brief Common list of pairs <first, last> cell numbers in each object
        * @return Const copy of pairs vector
    */ 
    const std::vector<std::array<std::array<double, 3>, 2>>& GetXBound() const {return _x_bound;}


    /**
        * @brief Cell type
        * @return Vertex count in CellType
    */ 
    constexpr int GetNVertex() const {return _NVertex;}





    /**
        * @brief Formation cell list of certain object
        * @details {cells_j}, j = (first_cell_in_obj, last_cell_in_obj)
        * @param objN object number
        * @return vector of object objN cells (one-typed)
    */
    std::vector<CellType> FormObjCellList(int objN) const;




    /**
        * @brief Formation of all basis functions on each cell in an array.
        * @details {tau_i^1, tau_i^2}, i = (1, num_frm)
        * @param surface_basis result set of all basis coordinates
        * @note memory must be allocated correctly: [num_frm][2][3]
    */ 
    void FormSurfaceBasis(double*** surface_basis) const;
    /**
        * @brief Formation of all basis functions on each cell in an array.
        * @return vector of all basis coordinates
        * @details {tau_i^1, tau_i^2}, i = (1, num_frm)
    */ 
    std::vector<std::array<std::array<double, 3>, 2>> FormSurfaceBasis() const;



    /**
        * @brief Formation of cell basis functions on certain object in an array.
        * @details {tau_i^1, tau_i^2}, i = (1_Obj, last_Obj)
        * @param objN obect number
        * @param surface_basis set of all basis coordinates on certain object
        * @note memory must be allocated correctly: [ObjNFrm[objN]][2][3]
    */ 
    void FormObjSurfaceBasis(int objN, double*** surface_basis) const;
    /**
        * @brief Formation of cell basis functions on certain object in an array.
        * @details {tau_i^1, tau_i^2}, i = (1_Obj, last_Obj)
        * @param objN obect number
        * @return vector of all basis coordinates on certain object
    */ 
    std::vector<std::array<std::array<double, 3>, 2>> FormObjSurfaceBasis(int objN) const;









    /**
        * @brief Formation of all collocation points on each cell in an array.
        * @details {rkt_i}, i = 1, num_frm
        * @param colloc_points set of all collocation points coordinates
        * @note memory must be allocated correctly: [num_frm][3]
    */
    void FormCollocPoints(double** colloc_points) const;
    /**
        * @brief Formation of all collocation points on each cell in an array.
        * @details {rkt_i}, i = 1, num_frm
        * @return vector of all collocation points coordinates
    */
    std::vector<std::array<double, 3>> FormCollocPoints() const;




    /**
        * @brief Formation of all normal vectors on each cell in an array.
        * @details {norm_i}, i = 1, num_frm
        * @param colloc_points set of all normal vectors coordinates
        * @note memory must be allocated correctly: [total_num_frm][3]
    */
    void FormNormalVectors(double** n_vec) const;
    /**
        * @brief Formation of all normal vectors on each cell in an array.
        * @details {norm_i}, i = 1, num_frm
        * @return set of all normal vectors coordinates
    */
    std::vector<std::array<double, 3>> FormNormalVectors() const;



    /**
        * @brief Formation of all vertex coordinates on each cell in an array.
        * @details {cell_i[4][3]}, i = 1, num_frm
        * @param cell_list_coords set of all collocation points coordinates
        * @note memory must be allocated correctly: [num_frm][NVertex][3]
    */
    void FormCellListCoords(double*** cell_list_coords) const;



    /**
        * @brief Formation of all vertex coordinates on each cell in an array.
        * @details {cell_i[4][3]}, i = 1, num_frm
        * @return set of all collocation points coordinates
    */
    auto FormCellListCoords() const {
        const size_t _CellPointsN = (std::is_same<CellType, CellQuad>::value) ? 4 : 3;
        std::vector<std::array<std::array<double, 3>, _CellPointsN>> cells;
        cells.reserve(_num_frm);
    
        for (int i = 0; i < _num_frm; i++) {
            std::array<std::array<double, 3>, _CellPointsN> cell;
            for (int j = 0; j < _CellPointsN; j++) {
                const double* vertex = _cell_list[i].GetVertex(j);
                std::copy(vertex, vertex + 3, cell[j].begin());
            }
            cells.push_back(std::move(cell));
        }
        return cells;
    }
};



















template<typename CellType>
void SurfMesh<CellType>::check_step(const CellType& cell, double cur_grid_step) {
    double diam_cur = get_diam(cell.GetCellCoords());
    if (diam_cur > cur_grid_step) {
        this->_grid_step = diam_cur;
    }
}



static bool endsWith(const std::string& str, const std::string& suffix) {
    if (suffix.length() > str.length()) {
        return false;
    }
    return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
}












//===========================================================
//------------------Reading from file------------------------
//===========================================================
template<typename CellType>
void SurfMesh<CellType>::ReadFromFile(const std::string& filename) {
    if (endsWith(filename, ".dat")) {
        // std::cout << "Чтение .dat файла" << std::endl;
        this->ReadFromDatFile(filename);
    } else if (endsWith(filename, ".vtk")) {
        // std::cout << "Чтение .vtk файла" << std::endl;
        this->ReadFromVtkFile(filename);
    } else {
        throw std::runtime_error("File type error: ReadFromFile");
    }
}





















// Пока что геометрия из одного объекта и одного модуля
template<typename CellType>
void SurfMesh<CellType>::ReadFromVtkFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Не удалось открыть файл: " + filename);

    std::string line;
    std::vector<std::array<double, 3>> points;
    std::vector<std::vector<int>> polygons;
    std::vector<int> objIds;
    std::string dataset_type;


    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string kw;
        iss >> kw;

        // Определяем тип датасета
        if (kw == "DATASET") {
            iss >> dataset_type;
            continue;
        }

        // Чтение точек
        if (kw == "POINTS") {
            int npts;
            std::string dtype;
            iss >> npts >> dtype;
            points.resize(npts);
            for (int i = 0; i < npts; ++i) {
                if (!std::getline(file, line)) {
                    throw std::runtime_error("Unexpected end of file while reading points");
                }
                std::istringstream pss(line);
                if (!(pss >> points[i][0] >> points[i][1] >> points[i][2])) {
                    throw std::runtime_error("Failed to parse point at index " + std::to_string(i));
                }
            }
        } else if (kw == "POLYGONS") { // Чтение полигонов для POLYDATA
            int npoly, total_values;
            iss >> npoly >> total_values;
            polygons.resize(npoly);
            for (int i = 0; i < npoly; ++i) {
                if (!std::getline(file, line)) {
                    throw std::runtime_error("Unexpected end of file while reading polygons");
                }
                std::istringstream pss(line);
                int n;
                if (!(pss >> n)) {
                    throw std::runtime_error("Failed to read polygon vertex count at index " + std::to_string(i));
                }
                polygons[i].resize(n);
                for (int j = 0; j < n; ++j) {
                    if (!(pss >> polygons[i][j])) {
                        throw std::runtime_error("Failed to read polygon vertices at index " + std::to_string(i));
                    }
                }
            }
        } else if (kw == "CELLS") { // Чтение ячеек для UNSTRUCTURED_GRID
            int ncells, total_values;
            iss >> ncells >> total_values;
            polygons.resize(ncells);
            for (int i = 0; i < ncells; ++i) {
                if (!std::getline(file, line)) {
                    throw std::runtime_error("Unexpected end of file while reading cells");
                }
                std::istringstream css(line);
                int n;
                if (!(css >> n)) {
                    throw std::runtime_error("Failed to read cell vertex count at index " + std::to_string(i));
                }
                polygons[i].resize(n);
                for (int j = 0; j < n; ++j) {
                    if (!(css >> polygons[i][j])) {
                        throw std::runtime_error("Failed to read cell vertices at index " + std::to_string(i));
                    }
                }
            }
        } else if (kw == "CELL_TYPES") { // Пропускаем CELL_TYPES
            int ncells;
            iss >> ncells;
            for (int i = 0; i < ncells; ++i) {
                std::getline(file, line); // просто пропускаем
            }
        } else if (kw == "CELL_DATA" || kw == "POINT_DATA") { // Чтение POINT_DATA или CELL_DATA с SCALARS
            int ndata;
            iss >> ndata;
            
            // Ищем следующую строку со SCALARS
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue;
                
                std::istringstream lss(line);
                std::string scalar_kw;
                lss >> scalar_kw;
                
                if (scalar_kw == "SCALARS") {
                    std::string name;
                    lss >> name;
                    
                    if (name == "Object") {
                        std::getline(file, line); // LOOKUP_TABLE
                        objIds.resize(polygons.size());
                        for (size_t i = 0; i < polygons.size(); ++i) {
                            if (!std::getline(file, line)) {
                                throw std::runtime_error("Unexpected end of file while reading object IDs");
                            }
                            try {
                                objIds[i] = std::stoi(line);
                            } catch (...) {
                                throw std::runtime_error("Failed to parse object ID at index " + std::to_string(i));
                            }
                        }
                        break;
                    }
                } else if (scalar_kw == "VECTORS" || scalar_kw == "LOOKUP_TABLE") {
                    continue;
                }
            }
            break; // выходим из основного цикла после обработки данных
        }
    }

    file.close();

    if (points.empty() || polygons.empty()) {
        throw std::runtime_error("VTK file contains no points or cells");
    }

    // Заполняем сетку
    this->_cell_list.clear();
    this->FillMeshFromVTK(polygons, points);

    this->_num_frm = this->_cell_list.size();

    
    // Определяем объекты
    if (objIds.empty()) {
        this->_num_obj = 1;
        this->_obj_num_frm.assign(1, this->_num_frm);
        this->_beg_endObjFrm.push_back({0, this->_num_frm - 1});
    } else {
        if ((int)objIds.size() != this->_num_frm) {
            throw std::runtime_error("Number of object IDs doesn't match number of cells");
        }
        
        int maxObj = *std::max_element(objIds.begin(), objIds.end());
        this->_num_obj = maxObj + 1;
        this->_obj_num_frm.assign(this->_num_obj, 0);
        this->_beg_endObjFrm.assign(this->_num_obj, {-1, -1});
        
        for (int i = 0; i < (int)objIds.size(); ++i) {
            int obj = objIds[i];
            this->_obj_num_frm[obj]++;
            if (this->_beg_endObjFrm[obj].first == -1) {
                this->_beg_endObjFrm[obj].first = i;
            }
            this->_beg_endObjFrm[obj].second = i;
        }
    }

    this->_num_mod = 1;
    this->_beg_endMod.assign(1, {0, this->_num_frm - 1});
    this->_beg_endObj.assign(1, {0, this->_num_obj - 1});
}




















//===========================================================
//-------------Formation of auxiliary data-------------------
//===========================================================
template<typename CellType>
std::vector<CellType> SurfMesh<CellType>::FormObjCellList(int objN) const {
    std::vector<CellType> object_cells;
    object_cells.resize(_obj_num_frm[objN]);
    int frm_start = _beg_endObjFrm[objN].first;
    int frm_end = _beg_endObjFrm[objN].second;
    object_cells.insert(object_cells.end(), _cell_list[frm_start], _cell_list[frm_end]);
    return object_cells;
}







template<typename CellType>
void SurfMesh<CellType>::FormSurfaceBasis(double*** surface_basis) const
{
    for (int i = 0; i < _num_frm; i++) {
        const double* tau1 = _cell_list[i].GetTau1();
        const double* tau2 = _cell_list[i].GetTau2();
        for (int j = 0; j < 3; j++) {
            surface_basis[i][0][j] = tau1[j];
            surface_basis[i][1][j] = tau2[j];
        }
    }
}


template<typename CellType>
std::vector<std::array<std::array<double, 3>, 2>>  SurfMesh<CellType>::FormSurfaceBasis() const
{
    std::vector<std::array<std::array<double, 3>, 2>> surface_basis;
    surface_basis.resize(_num_frm);
    for (int i = 0; i < _num_frm; i++) {
        const double* tau1 = _cell_list[i].GetTau1();
        const double* tau2 = _cell_list[i].GetTau2();
        for (int j = 0; j < 3; j++) {
            surface_basis[i][0][j] = tau1[j];
            surface_basis[i][1][j] = tau2[j];
        }
    }
    return surface_basis;
}


template<typename CellType>
void SurfMesh<CellType>::FormObjSurfaceBasis(int objN,
                                    double*** surface_basis) const
{
    int frm_start = _beg_endObjFrm[objN].first;
    int frm_end = _beg_endObjFrm[objN].second;

    for (int i = frm_start; i <= frm_end; i++) {
        double* tau1 = _cell_list[i].GetTau1();
        double* tau2 = _cell_list[i].GetTau2();
        for (int j = 0; j < 3; j++) {
            surface_basis[i - frm_start][0][j] = tau1[j];
            surface_basis[i - frm_start][1][j] = tau2[j];
        }
    }
}

template<typename CellType>
std::vector<std::array<std::array<double, 3>, 2>>
                    SurfMesh<CellType>::FormObjSurfaceBasis(int objN) const
{
    int ObjNFrm = _obj_num_frm[objN];
    std::vector<std::array<std::array<double, 3>, 2>> obj_surface_basis;
    obj_surface_basis.resize(ObjNFrm);

    int start_frm = _beg_endObjFrm[objN].first; 
    int end_frm = _beg_endObjFrm[objN].second;
    for (int i = start_frm; i <= end_frm; i++) {
        const double* tau1 = _cell_list[i].GetTau1();
        const double* tau2 = _cell_list[i].GetTau2();
        for (int j = 0; j < 3; j++) {
            obj_surface_basis[i - start_frm][0][j] = tau1[j];
            obj_surface_basis[i - start_frm][1][j] = tau2[j];
        }
    }
    return obj_surface_basis;
}













template<typename CellType>
void SurfMesh<CellType>::FormCollocPoints(double** colloc_points) const
{
    for (int i = 0; i < _num_frm; i++) {
        const double* rkt =  _cell_list[i].GetRkt();
        colloc_points[i][0] = rkt[0];
        colloc_points[i][1] = rkt[1];
        colloc_points[i][2] = rkt[2];
    }
}


template<typename CellType>
std::vector<std::array<double, 3>>  SurfMesh<CellType>::FormCollocPoints() const
{
    std::vector<std::array<double, 3>> colloc_points;
    colloc_points.resize(_num_frm);
    for (int i = 0; i < _num_frm; i++) {
        const double* rkt =  _cell_list[i].GetRkt();
        colloc_points[i][0] = rkt[0];
        colloc_points[i][1] = rkt[1];
        colloc_points[i][2] = rkt[2];
    }
    return colloc_points;
}










template<typename CellType>
void SurfMesh<CellType>::FormNormalVectors(double** n_vec) const
{
    for (int i = 0; i < _num_frm; i++) {
        const double* norm =  _cell_list[i].GetNorm();
        n_vec[i][0] = norm[0];
        n_vec[i][1] = norm[1];
        n_vec[i][2] = norm[2];
    }
}


template<typename CellType>
std::vector<std::array<double, 3>>  SurfMesh<CellType>::FormNormalVectors() const
{
    std::vector<std::array<double, 3>> n_vec;
    n_vec.resize(_num_frm);
    for (int i = 0; i < _num_frm; i++) {
        const double* norm =  _cell_list[i].GetNorm();
        n_vec[i][0] = norm[0];
        n_vec[i][1] = norm[1];
        n_vec[i][2] = norm[2];
    }
    return n_vec;
}





template<typename CellType>
void SurfMesh<CellType>::FormCellListCoords(double*** cell_list_coords) const
{
    int NVertex = (std::is_same<CellType, CellQuad>::value) ? 4 : 3;
    for (int i = 0; i < _num_frm; i++) {
        for (int j = 0; j < NVertex; j++) {
            const double* vertex = _cell_list[i].GetVertex(j);
            cell_list_coords[i][j][0] = vertex[0];
            cell_list_coords[i][j][1] = vertex[1];
            cell_list_coords[i][j][2] = vertex[2];
        }
    }
}
}       // namepace bielcc
#endif  // _SURF_MESH_H_