#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <utility>
#include <sstream>
#include <algorithm>

#include "cell.h"
#include "element_geom.h"
#include "surf_mesh.h"



namespace bielcc {
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
}       // namespace bielcc
