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
    const std::string geom_fname = argc > 1 ? argv[1]:
                  "../examples/data/surface_mesh/sphere_quad.dat";
    const std::string geom_fname_tr = argc > 1 ? argv[1]:
                  "../examples/data/surface_mesh/sphere_triangle.dat";
    auto start = std::chrono::high_resolution_clock::now();


    //==============1. Mesh object creation================
    std::cout << "QUADRANGULAR SURFACE MESH" << std::endl;
    bielcc::SurfMesh<CellQuad> simple_mesh(geom_fname);
    bielcc::RwgSurfMesh<CellQuad> rwg_mesh(simple_mesh);
    std::cout << "Cells number: " << rwg_mesh.GetNFrm() << std::endl;
    std::cout << "Unique points number: " << rwg_mesh.GetNPoints() << std::endl;
    std::cout << "Unique edges number: " << rwg_mesh.GetNEdges() << std::endl;



    //========Forming cell coordinates in one 3D vector===
    std::cout << "SOME CELLS FROM SURFACE QUAD MESH" << std::endl;
    const std::vector<std::array<int, 4>> cells_pnts = rwg_mesh.GetCellsPoints();
    for (int i = 0; i < cells_pnts.size(); i++) {
        if (i % 100 == 0) {
            std::cout << "Cell: " << i;
            std::cout << "      " << cells_pnts[i][0] << " " << cells_pnts[i][1];
            std::cout << " " << cells_pnts[i][2] << " " << cells_pnts[i][3] << std::endl;
        }
    }









    //==============2. Mesh object creation================
    std::cout << std::endl;
    std::cout << "TRIANGULAR SURFACE MESH" << std::endl;
    bielcc::SurfMesh<CellTriangle> simple_mesh_tr(geom_fname_tr);
    bielcc::RwgSurfMesh<CellTriangle> rwg_mesh_tr(simple_mesh_tr);
    std::cout << "Cells number: " << rwg_mesh_tr.GetNFrm() << std::endl;
    std::cout << "Unique points number: " << rwg_mesh_tr.GetNPoints() << std::endl;
    std::cout << "Unique edges number: " << rwg_mesh_tr.GetNEdges() << std::endl;


    //========Forming cell coordinates in one 3D vector===
    std::cout << "SOME CELLS FROM SURFACE TRIANGULAR MESH" << std::endl;
    const std::vector<std::array<int, 3>> cells_pnts_tr = rwg_mesh_tr.GetCellsPoints();
    for (int i = 0; i < cells_pnts_tr.size(); i++) {
        if (i % 100 == 0) {
            std::cout << "Cell: " << i;
            std::cout << "      " << cells_pnts_tr[i][0] << " " << cells_pnts_tr[i][1];
            std::cout << " " << cells_pnts_tr[i][2] << " " << cells_pnts_tr[i][3] << std::endl;
        }
    }









    // Засекаем время окончания выполнения программы
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    return 0;
}