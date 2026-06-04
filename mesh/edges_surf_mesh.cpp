#include <iostream>
#include <string>
#include <vector>
#include <array>

#include "edges_surf_mesh.h"
#include "cell.h"
#include "rwg_surf_mesh.h"
#include "point.h"
#include "element_geom.h"


namespace bielcc {
//===========================================================
//------------Quad surface mesh related to nodes-------------
//===========================================================
void EdgeSurfMesh::FillMesh()
{
    _cell_listG.resize(this->_n_edges);    
    _diag_G.resize(this->_n_edges); 
    _muG_i.resize(this->_n_edges);


    double root[4][3]{};
    //Заполнение muG, rutG
    for (int i = 0; i < this->_n_edges; i++) {
        double mu_sum = 0.;
        int a = this->_edges[i].GetStartIndx();
        int b = this->_edges[i].GetEndIndx();

        int n = this->_edge_num_cells[i];             // в скольких ячейках лежит ребро
        _cell_listG[i].reserve(n);
        _diag_G[i].reserve(n);
        _muG_i[i].reserve(n);
        for (int j = 0; j < n; j++) {
            int cell_n = this->_edge_cells[i][j];     // номер ячейки-соседа 'j' для ребра 'i'
            int edgeaN = -1, edgebN = -1;

            for (int k = 0; k < this->_NVertex; k++) {
                int edge_n = this->_cell_edges[cell_n][k]; // номер текущего отрезка из общего списка
                if (edge_n >= 0) {
                    // гарантия не пустого ребра(что две вершины не совпадают)
                    int startIndx = this->_edges[edge_n].GetStartIndx();
                    int endIndx = this->_edges[edge_n].GetEndIndx();
                    if ((startIndx == a && endIndx != b) ||
                        (endIndx == a && startIndx != b)) {
                        edgeaN = edge_n;                    //номер отрезка, на котором лежит точка a, но не лежит точка b
                    }

                    if ((startIndx == b && endIndx != a) ||
                        (endIndx == b && startIndx != a)) {
                        edgebN = edge_n; // номер отрезка, на котором лежит точка b, но не лежит точка a
                    }
                }
            }

            if (edgeaN == -1 || edgebN == -1) {
                std::cout << "Ошибка в заполнениии edges_surf_mesh. edgeaN или edgebN = -1" << std::endl;
                exit(-1);
            }


            // Заполнение rutG
            for (int g = 0; g < 3; g++) {
                root[0][g] = this->_points[a][g];
                root[1][g] = this->_points[b][g];
                root[2][g] = this->_edge_mid[edgebN][g];
                root[3][g] = this->_edge_mid[edgeaN][g];
            }
            CellQuad cell_tmp;
            cell_tmp.CellFill(root);
            _cell_listG[i].push_back(cell_tmp);

            // Диагональ 
            _diag_G[i].push_back(get_diam(root));
            _muG_i[i].push_back(cell_tmp.GetArea());

            //Вычисление площади
            mu_sum += cell_tmp.GetArea();   
        }
        _muG.push_back(mu_sum);
    }
}

}          // namespace bielcc