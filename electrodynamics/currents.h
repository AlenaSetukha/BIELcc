#ifndef _CURRENTS_H_
#define _CURRENTS_H_

#include <iostream>
#include <string>
#include <complex>
#include <string>
#include <vector>
#include <array>

namespace bielcc {
//===========================================================
//--------Currents on cells through a local basis------------
//===========================================================
/**
    * Formation of currents on cells using local cell bases
    * @note j(sigma_i) = j_i^1 * tau_i^1 + j_i^1 * tau_i^1
    * @param tau set of coordinates of local bases (2 vectors) within each cell
    * @param b system coefficients (coefficients solution answer)
    * @param j_vec currents on cells
*/ 
template <typename CurType>
void SurfaceCurrents_byBasis(const std::vector<std::array<std::array<double, 3>, 2>> &tau,
                             const CurType* b, std::vector<std::array<CurType, 3>> &j_vec)
{
    int i0, num_frm = tau.size();
    for (int i = 0; i < num_frm; i++) {
        i0 = 2 * i;
        j_vec[i][0] = b[i0] * tau[i][0][0] + b[i0 + 1] * tau[i][1][0];
        j_vec[i][1] = b[i0] * tau[i][0][1] + b[i0 + 1] * tau[i][1][1];
        j_vec[i][2] = b[i0] * tau[i][0][2] + b[i0 + 1] * tau[i][1][2];
    }
}

template <typename CurType>
void SurfaceCurrents_byBasis(int num_frm, const double*** tau, const CurType* b, CurType** j_vec)
{
    int i0;
    for (int i = 0; i < num_frm; i++) {
        i0 = 2 * i;
        j_vec[i][0] = b[i0] * tau[i][0][0] + b[i0 + 1] * tau[i][1][0];
        j_vec[i][1] = b[i0] * tau[i][0][1] + b[i0 + 1] * tau[i][1][1];
        j_vec[i][2] = b[i0] * tau[i][0][2] + b[i0 + 1] * tau[i][1][2];
    }
}












//===========================================================
//-----------Reading a current vector from file--------------
//===========================================================
/**
    * @brief Reading real current vector from file
    * @param filename_real file with currents values (3D)
    * @param j_vec currents on cells
*/ 
void SurfaceCurrentsD_fromFile(const std::string &filename_real,
                      std::vector<std::array<double, 3>> &j_vec);
                      
void SurfaceCurrentsD_fromFile(const std::string &filename_real,
                                                 double** j_vec);

/**
    * @brief Reading complex current vector from file
    * @param filename_real  file with real part currents values (3D)
    * @param filename_image file with image part currents values (3D)
    * @param j_vec currents on cells
*/
void SurfaceCurrentsC_fromFiles(const std::string &filename_real,
                               const std::string &filename_image,
        std::vector<std::array<std::complex<double>, 3>> &j_vec);

void SurfaceCurrentsC_fromFiles(const std::string &filename_real,
                               const std::string &filename_image,
                                    std::complex<double>** j_vec);


}       // namespace bielcc
#endif  // _CURRENTS_H_
