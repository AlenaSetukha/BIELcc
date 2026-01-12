#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>

#include "cell.h"
#include "surf_mesh.h"
#include "point.h"
#include "rwg_surf_mesh.h"
#include "segment.h"

using namespace bielcc;

int main(int argc, char **argv)
{
    const int NVertex = 3;

    const std::string geom_fname = argc > 1 ? argv[1]:
                  "../examples/data/grids/30_50_triangle.dat";
    auto start = std::chrono::high_resolution_clock::now();

    //==============Mesh object creation================
    bielcc::SurfMesh<CellTriangle> simple_mesh(geom_fname);
    std::cout << "Всего ячеек: " << simple_mesh.GetNFrm() << std::endl;
    bielcc::RwgSurfMesh<CellTriangle> rwg_mesh(simple_mesh);
    std::cout << "Всего уникальных точек: " << rwg_mesh.GetNPoints() << std::endl;
    std::cout << "Всего уникальных ребер: " << rwg_mesh.GetNEdges() << std::endl;




    //========Forming cell coordinates in one 3D vector===
    const std::vector<std::array<int, 3>> cell_pnts = rwg_mesh.GetCellPoints();
    for (int i = 0; i < cell_pnts.size(); i++) {
        std::cout << "Cell: " << i;
        std::cout << "      " << cell_pnts[i][0] << " " << cell_pnts[i][1] << " " << cell_pnts[i][2] << std::endl;
    }







    // Засекаем время окончания выполнения программы
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    return 0;
}