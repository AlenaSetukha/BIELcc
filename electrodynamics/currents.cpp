#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <vector>
#include <array>
#include <iomanip>
#include <limits>
#include <numbers>


#include "currents.h"

namespace bielcc {
void SurfaceCurrentsD_fromFile(const std::string &filename_real,
                      std::vector<std::array<double, 3>> &j_vec)
{
    std::ifstream fin_real(filename_real);

    if (!fin_real.is_open()) {
        throw std::runtime_error("Read currents.txt error: SurfaceCurrentsD_fromFile");
    }

    const auto digits = std::numeric_limits<double>::digits10;
    fin_real >> std::fixed >> std::setprecision(digits);

    int n1_real, n2_real, n;
    fin_real >> n1_real >> n2_real;

    n = n1_real * n2_real;
    double x1_real, x2_real, x3_real;
    for (int i = 0; i < n; i++) {
        fin_real >> x1_real >> x2_real >> x3_real;
        j_vec[i][0] = x1_real;
        j_vec[i][1] = x2_real;
        j_vec[i][2] = x3_real;
    }
    fin_real.close();
}
          

void SurfaceCurrentsD_fromFile(const std::string &filename_real,
                                                 double** j_vec)
{
    std::ifstream fin_real(filename_real);

    if (!fin_real.is_open()) {
        throw std::runtime_error("Read currents.txt error: SurfaceCurrentsD_fromFile");
    }

    const auto digits = std::numeric_limits<double>::digits10;
    fin_real >> std::fixed >> std::setprecision(digits);

    int n1_real, n2_real, n;
    fin_real >> n1_real >> n2_real;

    n = n1_real * n2_real;
    double x1_real, x2_real, x3_real;
    for (int i = 0; i < n; i++) {
        fin_real >> x1_real >> x2_real >> x3_real;
        j_vec[i][0] = x1_real;
        j_vec[i][1] = x2_real;
        j_vec[i][2] = x3_real;
    }
    fin_real.close();
}







void SurfaceCurrentsC_fromFiles(const std::string &filename_real,
                               const std::string &filename_image,
        std::vector<std::array<std::complex<double>, 3>> &j_vec)
{
    std::ifstream fin_real(filename_real);
    std::ifstream fin_image(filename_image);

    if (!fin_real.is_open() || !fin_image.is_open()) {
        throw std::runtime_error("Reading currents.txt error: SurfaceCurrentsC_fromFiles");
    }


    const auto digits = std::numeric_limits<double>::digits10;
    fin_real >> std::fixed >> std::setprecision(digits);
    fin_image >> std::fixed >> std::setprecision(digits);

    int n1_real, n2_real, n1_image, n2_image, n;
    fin_real >> n1_real >> n2_real;
    fin_image >> n1_image >> n2_image;

    if (n1_real != n1_image || n2_real != n2_image) {
        throw std::runtime_error("Reading currents.txt error: SurfaceCurrentsC_fromFiles");
    }
    n = n1_real * n2_real;
    double x1_real, x2_real, x3_real, x1_image, x2_image, x3_image;
    for (int i = 0; i < n; i++) {
        fin_real >> x1_real >> x2_real >> x3_real;
        fin_image >> x1_image >> x2_image >> x3_image;
        j_vec[i][0] = std::complex<double>(x1_real, x1_image);
        j_vec[i][1] = std::complex<double>(x2_real, x2_image);
        j_vec[i][2] = std::complex<double>(x3_real, x3_image);
    }

    fin_real.close();
    fin_image.close();
}

void SurfaceCurrentsC_fromFiles(const std::string &filename_real,
                               const std::string &filename_image,
                                    std::complex<double>** j_vec)
{
    std::ifstream fin_real(filename_real);
    std::ifstream fin_image(filename_image);

    if (!fin_real.is_open() || !fin_image.is_open()) {
        throw std::runtime_error("Reading currents.txt error: SurfaceCurrentsC_fromFiles");
    }


    const auto digits = std::numeric_limits<double>::digits10;
    fin_real >> std::fixed >> std::setprecision(digits);
    fin_image >> std::fixed >> std::setprecision(digits);

    int n1_real, n2_real, n1_image, n2_image, n;
    fin_real >> n1_real >> n2_real;
    fin_image >> n1_image >> n2_image;

    if (n1_real != n1_image || n2_real != n2_image) {
        throw std::runtime_error("Reading currents.txt error: SurfaceCurrentsC_fromFiles");
    }
    n = n1_real * n2_real;
    double x1_real, x2_real, x3_real, x1_image, x2_image, x3_image;
    for (int i = 0; i < n; i++) {
        fin_real >> x1_real >> x2_real >> x3_real;
        fin_image >> x1_image >> x2_image >> x3_image;
        j_vec[i][0] = std::complex<double>(x1_real, x1_image);
        j_vec[i][1] = std::complex<double>(x2_real, x2_image);
        j_vec[i][2] = std::complex<double>(x3_real, x3_image);
    }

    fin_real.close();
    fin_image.close();
}

} // namespace bielcc


