#include <iostream>
#include <string>
#include <vector>
#include <array>

#include "nodes_surf_mesh.h"
#include "cell.h"
#include "rwg_surf_mesh.h"
#include "point.h"
#include "element_geom.h"


namespace bielcc {
//===========================================================
//------------Quad surface mesh related to nodes-------------
//===========================================================
void NodeSurfMesh::FillMesh()
{
    _cell_listD.resize(this->_n_points);
    _muD_i.resize(this->_n_points);
    
    
    double root[4][3]{};
    // Заполнение muD, rutD
    for (int i = 0; i < this->_n_points; i++) {
        double mu_sum = 0;
        int n = this->_point_num_cells[i];
        _cell_listD[i].reserve(n);
        _muD_i[i].resize(n);
        for (int j = 0; j < n; j++) {
            int cell_n = this->_point_cells[i][j]; // номер ячейки-соседа 'j' для точки 'i'
            double c_m[3]{};
            CellQuad cellq = this->_cell_list[cell_n];
            get_center_mass(cellq.GetCellCoords(), c_m);
            int n1 = -1, n2 = -1;
            for (int k = 0; k < this->_NVertex; k++) {
                int edge_n = this->_cell_edges[cell_n][k]; // номер текущего отрезка из общего списка
                if (edge_n >= 0) {
                    // ребро не пустое(не совпадающие точки)
                    int node1 = this->_edges[edge_n].GetStartIndx();
                    int node2 = this->_edges[edge_n].GetEndIndx();
                    if (node1 == i || node2 == i) {
                        n2 = n1;      // номер второго ребра, на котором лежит точка
                        n1 = edge_n;  // номер первого ребра, на котором лежит точка
                    }
                }
            }
            if ((n2 == -1) || (n1 == -1)) {
                std::cout << "Ошибка в заполнении nodes_surf_mesh.fill" << std::endl;
                std::cout << "Code 1" << std::endl;
                exit(-1);
            }
        

            // Заполнение rutD = cellD
            for (int g = 0; g < 3; g++) {
                root[0][g] = this->_points[i][g];
                root[1][g] = this->_edge_mid[n1][g];
                root[2][g] = c_m[g];
                root[3][g] = this->_edge_mid[n2][g];
            }
            CellQuad cell_tmp;
            cell_tmp.CellFill(root);
            _cell_listD[i].push_back(cell_tmp);


            // Вычисление площади
            mu_sum += cell_tmp.GetArea();
            _muD_i[i][j] = cell_tmp.GetArea();
        }
        _muD.push_back(mu_sum);
    }
}
}          // namespace bielcc