#include <iostream>
#include <fstream>
#include <string>

#include "surf_mesh.h"
#include "edges_surf_mesh.h"
#include "nodes_surf_mesh.h"
#include "cell.h"



namespace bielcc {
//===========================================================
//-------Создание geodat.dat из голубой сетки----------------
//===========================================================
void MakeDatFromEdges(const std::string& res_dir, int NObj,
                                int NMod, const EdgeSurfMesh& b)
{
    std::ofstream fout(res_dir + "geodat_blue.dat");
    if (!fout.is_open()) {
        std::cout << "Open " + res_dir + "geodat_blue.dat " +  " error" << std::endl;
        return;
    }

    fout << NObj << std::endl; //всего объектов
    fout << NMod << std::endl; //всего модулей
    fout << 1 << std::endl;

    int num_edge_sum = 0;
    for (int i = 0; i < b.GetNEdges(); i++) {
        num_edge_sum += b.GetEdgeNumCells(i);
    }
    fout << num_edge_sum << std::endl;
    fout << 1 << " " << 1 << std::endl;


    fout << std::fixed;
    fout.precision(7);
    int n;
    for (int i = 0; i < b.GetNEdges(); i++) {
        n = b.GetEdgeNumCells(i);
        for (int j = 0; j < n; j++) {
            double cell_tmp[4][3];
            b.GetEdgeCellG(i, j).GetCellCoords(cell_tmp);
            for (int k = 0; k < 4; k++) {
                fout << cell_tmp[k][0] << " " << cell_tmp[k][1] << " " << cell_tmp[k][2] << " ";
            }
            fout << std::endl;
        }
    }
    fout.close();
}





//===========================================================
//-------Создание geodat.dat из серой сетки------------------
//===========================================================
void MakeDatFromNodes(const std::string& res_dir, int NObj,
                                int NMod, const NodeSurfMesh& d)
{
    std::ofstream fout(res_dir + "geodat_grey.dat");
    if (!fout.is_open()) {
        std::cout << "Open " + res_dir + "geodat_grey.dat " +  " error" << std::endl;
        return;
    }

    fout << NObj << std::endl; //всего объектов
    fout << NMod << std::endl; //всего модулей
    fout << 1 << std::endl;

    int num_points_sum = 0;
    for (int i = 0; i < d.GetNPoints(); i++) {
        num_points_sum += d.GetPointNumCells(i);
    }
    fout << num_points_sum << std::endl;
    fout << 1 << " " << 1 << std::endl;

    fout << std::fixed;
    fout.precision(7);

    int n;

    for (int i = 0; i < d.GetNPoints(); i++) {
        n = d.GetPointNumCells(i);
        for (int j = 0; j < n; j++) {
            double cell_tmp[4][3];
            d.GetNodeCellD(i, j).GetCellCoords(cell_tmp);
            for (int k = 0; k < 4; k++) {
                fout << cell_tmp[k][0] << " " << cell_tmp[k][1] << " " << cell_tmp[k][2] << " ";
            }
            fout << std::endl;
        }
    }

    fout.close();
    return;
}

}       // namespace bielcc