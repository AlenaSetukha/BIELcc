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
 * @brief Surface mesh class (dim = 3) consisting of quadrangular or triangular cells.
 * @note Cells order corresponds file cells order
 */
template<typename CellType>
class SurfMesh {
protected:
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
    SurfMesh(const std::string& filename) { this->ReadFromFile(filename);};
    SurfMesh(): _num_obj(0), _num_mod(0), _num_frm(0), _num_bl(0), _num_blSeg(0), _grid_step(0.0) {}
    SurfMesh(const SurfMesh& obj): _cell_list(obj._cell_list), _num_obj(obj._num_obj), 
                        _num_mod(obj._num_mod), _num_frm(obj._num_frm),
                        _num_bl(obj._num_bl), _num_blSeg(obj._num_blSeg), 
                        _grid_step(obj._grid_step), _surf_square(obj._surf_square),
                        _beg_endMod(obj._beg_endMod), 
                        _beg_endObj(obj._beg_endObj), _x_bound(obj._x_bound) {}



    void ReadFromFile(const std::string& filename);

    SurfMesh& operator=(const SurfMesh& obj);

    ~SurfMesh() = default;



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







    CellType GetCell(int i) {return _cell_list[i];}
    const CellType& GetCell(int i) const {return _cell_list[i];}


    std::vector<CellType> GetTotalCellList() {return _cell_list;}
    const std::vector<CellType>& GetTotalCellList() const {return _cell_list;}


    std::vector<std::pair<int, int>> GetBegEndObj() {return _beg_endObj;}
    const std::vector<std::pair<int, int>>& GetBegEndObj() const {return _beg_endObj;}
    /**
        * @brief First and last module in object objN
        * @param objN object number
    */
    std::pair<int, int> GetFLModInObj(int objN) {return _beg_endObj[objN];}
    const std::pair<int, int>& GetFLModInObj(int objN) const {return _beg_endObj[objN];}


    std::vector<std::pair<int, int>> GetBegEndMod() {return _beg_endMod;}
    const std::vector<std::pair<int, int>>& GetBegEndMod() const {return _beg_endMod;}
    /**
        * @brief First and last cell in module modN
        * @param objN object number
    */
    std::pair<int, int> GetFLFrmInMod(int modN) {return _beg_endMod[modN];}
    const std::pair<int, int>& GetFLFrmInMod(int modN) const {return _beg_endMod[modN];}


    /**
        * @brief First and last cell in object objN
        * @param objN object number
    */
    std::pair<int, int> GetFLFrmInObj(int objN) {return _beg_endObjFrm[objN];}
    const std::pair<int, int>& GetFLFrmInObj(int objN) const {return _beg_endObjFrm[objN];}


    std::vector<std::array<std::array<double, 3>, 2>> GetXBound() {return _x_bound;}
    const std::vector<std::array<std::array<double, 3>, 2>>& GetXBound() const {return _x_bound;}







    /**
        * @brief Formation cell list of certain object
        * @details {cells_j}, j = (first_in_obj, last_in_obj)
        * @param objN object number
    */
    std::vector<CellType> FormObjCellList(int objN) const;




    /**
        * @brief Formation of basis functions on each cell in an array.
        * @details {tau_i^1, tau_i^2}, i = (1, num_frm)
        * @param surface_basis set of all basis coordinates
        * @note memory must be allocated correctly: [num_frm][2][3]
    */ 
    void FormTotalSurfaceBasis(double*** surface_basis) const;
    std::vector<std::array<std::array<double, 3>, 2>> FormTotalSurfaceBasis() const;

    /**
        * @brief Formation of cell basis functions on certain object in an array.
        * @details {tau_i^1, tau_i^2}, i = (1_Obj, last_Obj)
        * @param objN obect number
        * @param surface_basis set of all basis coordinates on certain object
        * @note memory must be allocated correctly: [GetObjNFrm[objN]][2][3]
    */ 
    void FormObjSurfaceBasis(int objN, double*** surface_basis) const;
    std::vector<std::array<std::array<double, 3>, 2>> FormObjSurfaceBasis(int objN) const;









    /**
        * @brief Formation of collocation points on each cell in an array.
        * @details {rkt_i}, i = 1, num_frm
        * @param colloc_points set of all collocation points coordinates
        * @note memory must be allocated correctly: [num_frm][3]
    */
    void FormCollocPoints(double** colloc_points) const;
    std::vector<std::array<double, 3>> FormCollocPoints() const;




    /**
        * @brief Formation of normal vectors on each cell in an array.
        * @details {norm_i}, i = 1, num_frm
        * @param colloc_points set of all normal vectors coordinates
        * @note memory must be allocated correctly: [num_frm][3]
    */
    void FormNormalVectors(double** n_vec) const;
    std::vector<std::array<double, 3>> FormNormalVectors() const;

    /**
        * @brief Formation of cvertex coordinates on each cell in an array.
        * @details {cell_i[4][3]}, i = 1, num_frm
        * @param cell_list_coords set of all collocation points coordinates
        * @note memory must be allocated correctly: [num_frm][NVertex][3]
    */
    void FormCellListCoords(double*** cell_list_coords) const;


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
        std::cout << "Чтение .dat файла" << std::endl;
        this->ReadFromDatFile(filename);
    } else if (endsWith(filename, ".vtk")) {
        std::cout << "Чтение .vtk файла" << std::endl;
        this->ReadFromVtkFile(filename);
    } else {
        throw std::runtime_error("File type error: ReadFromFile");
    }
}







template<>
void SurfMesh<CellTriangle>::ReadFromDatFile(const std::string& datFilename) {
    _grid_step = 0., _num_frm = 0, _num_mod = 0, _num_blSeg = 0;
    double Sobj = 0.;

    int kmod, kseg, sign, nfm, not_used;
    int cur_frm = 0, cur_mod = 0;

    int objNFrm = 0., f_frm_obj, l_frm_obj;

    std::ifstream fin(datFilename);
    if (!fin.is_open()) {
        throw std::runtime_error("Read geodat.dat error: Triangle constructor");
    }


    // Filling mesh
    CellTriangle cell_tmp;
    double root[3][3];

    fin >> _num_obj;
    for (int i = 0; i < _num_obj; i++) {
        objNFrm = 0.;
        Sobj = 0.;

        fin >> kmod;
        _num_mod += kmod;
        std::pair<int, int> beo_i = {cur_mod, 0};
        for (int j = 0; j < kmod; j++) {
            fin >> sign >> nfm >> not_used >> not_used;
            _num_frm += nfm;
            objNFrm += nfm;
            std::pair<int, int> be_i = {cur_frm, 0};
            for (int k = 0; k < nfm; k++) {
                for (int g = 0; g < 3; g++) {
                    fin >> root[g][0] >> root[g][1] >> root[g][2];
                }
                cell_tmp.CellFill(root);
                _cell_list.push_back(cell_tmp);
                Sobj += cell_tmp.GetArea();
                check_step(_cell_list[cur_frm], this->_grid_step);


                if (j == 0 && k == 0) {
                    f_frm_obj = cur_frm;
                }
                if (j == kmod - 1 && k == nfm - 1) {
                    l_frm_obj = cur_frm;
                }
                cur_frm++;
            }
            be_i.second = cur_frm - 1;
            _beg_endMod.push_back(be_i);
            cur_mod++;
        }
        beo_i.second = cur_mod - 1;
        _beg_endObj.push_back(beo_i);
        _surf_square.push_back(Sobj);
        _obj_num_frm.push_back(objNFrm);
        _beg_endObjFrm.push_back(std::pair<int, int> (f_frm_obj, l_frm_obj));
    }


    // Заполнений линий отрыва
    fin >> _num_bl;
    for (int i = 0; i < _num_bl; i++) {
        fin >> kseg;
        _num_blSeg += kseg;
        for (int j = 0; j < kseg; j++) {
            std::array<std::array<double, 3>, 2> xb_i;
            for (int k = 0; k < 2; k++) {
                fin >> xb_i[k][0] >> xb_i[k][1] >> xb_i[k][2];
            }
            _x_bound.push_back(xb_i);
        }
    }
    fin.close();
}






template<>
void SurfMesh<CellQuad>::ReadFromDatFile(const std::string& datFilename) {
    _grid_step = 0., _num_frm = 0, _num_mod = 0, _num_blSeg = 0;

    int kmod, kseg, sign, nfm, not_used;
    int cur_frm = 0, cur_mod = 0;
    double Sobj = 0.;

    int objNFrm = 0., f_frm_obj, l_frm_obj;

    std::ifstream fin(datFilename);
    if (!fin.is_open()) {
        throw std::runtime_error("Read geodat.dat error: Quad constructor");
    }

    // Заполнение сетки
    CellQuad cell_tmp;
    double root[4][3];

    fin >> _num_obj;
    for (int i = 0; i < _num_obj; i++) {
        objNFrm = 0.;
        fin >> kmod;
        _num_mod += kmod;
        std::pair<int, int> beo_i = {cur_mod, 0};

        for (int j = 0; j < kmod; j++) {
            fin >> sign >> nfm >> not_used >> not_used;
            _num_frm += nfm;
            objNFrm += nfm;
            std::pair<int, int> be_i = {cur_frm, 0};
            for (int k = 0; k < nfm; k++) {
                for (int g = 0; g < 4; g++) {
                    fin >> root[g][0] >> root[g][1] >> root[g][2];
                }
                cell_tmp.CellFill(root);
                _cell_list.push_back(cell_tmp);
                Sobj += cell_tmp.GetArea();
                check_step(_cell_list[cur_frm], this->_grid_step);
                if (j == 0 && k == 0) {
                    f_frm_obj = cur_frm;
                }
                if (j == kmod - 1 && k == nfm - 1) {
                    l_frm_obj = cur_frm;
                }
                cur_frm++;
            }
            be_i.second = cur_frm - 1;
            _beg_endMod.push_back(be_i);
            cur_mod++;
        }
        beo_i.second = cur_mod - 1;
        _beg_endObj.push_back(beo_i);
        _surf_square.push_back(Sobj);
        _obj_num_frm.push_back(objNFrm);
        _beg_endObjFrm.push_back(std::pair<int, int> (f_frm_obj, l_frm_obj));
    }


    // Заполнений линий отрыва
    fin >> _num_bl;
    for (int i = 0; i < _num_bl; i++) {
        fin >> kseg;
        _num_blSeg += kseg;
        for (int j = 0; j < kseg; j++) {
            std::array<std::array<double, 3>, 2> xb_i;
            for (int k = 0; k < 2; k++) {
                fin >> xb_i[k][0] >> xb_i[k][1] >> xb_i[k][2];
            }
            _x_bound.push_back(xb_i);
        }
    }
    fin.close();
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


template<>
void SurfMesh<CellQuad>::FillMeshFromVTK(const std::vector<std::vector<int>>& polygons,
                                        const std::vector<std::array<double, 3>>& points)
{
    CellQuad cell_tmp;
    int skipped_count = 0;
    _grid_step = 0.;
    double Sobj = 0.;
    
    for (size_t poly_idx = 0; poly_idx < polygons.size(); ++poly_idx) {
        const auto& poly = polygons[poly_idx];
        
        if (poly.size() == 4) { // четырёхугольник
            double cell[4][3];
            try {
                for (size_t i = 0; i < 4; ++i) {
                    int pt_idx = poly[i];
                    if (pt_idx < 0 || pt_idx >= (int)points.size()) {
                        throw std::runtime_error("Point index " + std::to_string(pt_idx) + 
                                               " out of range [0, " + std::to_string(points.size()-1) + "]");
                    }
                    cell[i][0] = points[pt_idx][0];
                    cell[i][1] = points[pt_idx][1];
                    cell[i][2] = points[pt_idx][2];
                }
                cell_tmp.CellFill(cell);
                this->_cell_list.push_back(cell_tmp);
                check_step(cell_tmp, this->_grid_step);
                Sobj += cell_tmp.GetArea();
            } catch (const std::exception& e) {
                throw std::runtime_error("Error processing cell at index " + std::to_string(poly_idx) + 
                                       ": " + std::string(e.what()));
            }
        } else {
            // Пропускаем ячейки с неправильным количеством вершин (треугольники, пятиугольники и т.д.)
            skipped_count++;
        }
    }
    _surf_square.push_back(Sobj);
    
    if (skipped_count > 0) {
        std::cerr << "Warning: Skipped " << skipped_count << " cells with " 
                  << "unsupported vertex count (expected 4)\n";
    }
    
    if (this->_cell_list.empty()) {
        throw std::runtime_error("No valid quadrilateral cells found in VTK file");
    }
}


template<>
void SurfMesh<CellTriangle>::FillMeshFromVTK(const std::vector<std::vector<int>>& polygons,
                                        const std::vector<std::array<double, 3>>& points)
{
    CellTriangle cell_tmp;
    int skipped_count = 0;
    _grid_step = 0.;
    double Sobj = 0.;
    
    for (size_t poly_idx = 0; poly_idx < polygons.size(); ++poly_idx) {
        const auto& poly = polygons[poly_idx];
        
        if (poly.size() == 3) { // треугольник
            double cell[3][3];
            try {
                for (size_t i = 0; i < 3; ++i) {
                    int pt_idx = poly[i];
                    if (pt_idx < 0 || pt_idx >= (int)points.size()) {
                        throw std::runtime_error("Point index " + std::to_string(pt_idx) + 
                                               " out of range [0, " + std::to_string(points.size()-1) + "]");
                    }
                    cell[i][0] = points[pt_idx][0];
                    cell[i][1] = points[pt_idx][1];
                    cell[i][2] = points[pt_idx][2];
                }
                cell_tmp.CellFill(cell);
                this->_cell_list.push_back(cell_tmp);
                check_step(cell_tmp, this->_grid_step);
                Sobj += cell_tmp.GetArea();
            } catch (const std::exception& e) {
                throw std::runtime_error("Error processing cell at index " + std::to_string(poly_idx) + 
                                       ": " + std::string(e.what()));
            }
        } else {
            // Пропускаем ячейки с неправильным количеством вершин (треугольники, пятиугольники и т.д.)
            skipped_count++;
        }
    }

    _surf_square.push_back(Sobj);
    
    if (skipped_count > 0) {
        std::cerr << "Warning: Skipped " << skipped_count << " cells with " 
                  << "unsupported vertex count (expected 4)\n";
    }
    
    if (this->_cell_list.empty()) {
        throw std::runtime_error("No valid quadrilateral cells found in VTK file");
    }
}













































//===========================================================
//----------------------Operators----------------------------
//===========================================================
template<typename CellType>
SurfMesh<CellType>& SurfMesh<CellType>::operator=(const SurfMesh<CellType>& obj) {
    if (this != &obj) {
        _cell_list = obj._cell_list;

        _num_obj = obj._num_obj;
        _num_mod = obj._num_mod;
        _num_frm = obj._num_frm;
        _num_bl = obj._num_bl;
        _num_blSeg = obj._num_blSeg;
        _grid_step = obj._grid_step;
        _surf_square = obj._surf_square;
        _obj_num_frm = obj._obj_num_frm;

        _beg_endMod = obj._beg_endMod;
        _beg_endObj = obj._beg_endObj;
        _beg_endObjFrm = obj._beg_endObjFrm;
        _x_bound = obj._x_bound;
    }
    return *this;
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
void SurfMesh<CellType>::FormTotalSurfaceBasis(double*** surface_basis) const
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
std::vector<std::array<std::array<double, 3>, 2>>  SurfMesh<CellType>::FormTotalSurfaceBasis() const
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

    
    
    int start_mod = _beg_endObj[objN].first; 
    int end_mod = _beg_endObj[objN].second;

    for (int i = start_mod; i <= end_mod; i++) {
        int f_cell = _beg_endMod[i].first;
        int l_cell = _beg_endMod[i].second;
        for (int j = f_cell; j < l_cell; j++) {

        }
        obj_surface_basis.insert(obj_surface_basis.end(), _cell_list);
    }

    for (int i = 0; i < ObjNFrm; i++) {
        const double* tau1 = _cell_list[i].GetTau1();
        const double* tau2 = _cell_list[i].GetTau2();
        for (int j = 0; j < 3; j++) {
            obj_surface_basis[i][0][j] = tau1[j];
            obj_surface_basis[i][1][j] = tau2[j];
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
    n_vec.reserve(_num_frm);
    for (int i = 0; i < _num_frm; i++) {
        const double* norm =  _cell_list[i].GetNorm();
        n_vec.push_back({norm[0], norm[1], norm[2]});
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


/*
std::vector<std::array<std::array<double, 3>, 4>>
                    SurfMesh<CellQuad>::FormCellListCoords() const
{
    std::vector<std::array<std::array<double, 3>, 4>> cells;
    cells.reserve(_num_frm);
    
    for (int i = 0; i < _num_frm; i++) {
        std::array<std::array<double, 3>, 4> cell;
        for (int j = 0; j < 4; j++) {
            const double* vertex = _cell_list[i].GetVertex(j);
            std::copy(vertex, vertex + 3, cell[j].begin());
        }
        cells.push_back(std::move(cell));
    }
    
    return cells;
}

std::vector<std::array<std::array<double, 3>, 3>>
                    SurfMesh<CellQuad>::FormCellListCoords() const
{
    std::vector<std::array<std::array<double, 3>, 3>> cells;
    cells.reserve(_num_frm);
    
    for (int i = 0; i < _num_frm; i++) {
        std::array<std::array<double, 3>, 3> cell;
        for (int j = 0; j < 3; j++) {
            const double* vertex = _cell_list[i].GetVertex(j);
            std::copy(vertex, vertex + 3, cell[j].begin());
        }
        cells.push_back(std::move(cell));
    }
    return cells;
}*/






}       // namepace bielcc
#endif  // _SURF_MESH_H_