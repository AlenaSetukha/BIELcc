#ifndef _MAKE_GEODAT_H_
#define _MAKE_GEODAT_H_

#include "edges_surf_mesh.h"
#include "nodes_surf_mesh.h"
#include "cell.h"

//Функции, создающие файл geodat.dat(Linux)
namespace bielcc {
    //===========================================================
    //-------Создание geodat.dat из голубой сетки----------------
    //===========================================================
    void MakeDatFromEdges(const std::string& res_dir, int NObj,
                                int NMod, const EdgeSurfMesh& b);


    //===========================================================
    //--------Создание geodat.dat из серой сетки-----------------
    //===========================================================
    void MakeDatFromNodes(const std::string& res_dir, int NObj,
                                int NMod, const NodeSurfMesh& b);
}

#endif