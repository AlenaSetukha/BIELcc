#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <vector>
#include <array>
#include <string>
#include <iomanip>

#include "io_aereco.h"


namespace bielcc {
//===========================================================
//----------Reading data from Aereco-type files--------------
//===========================================================
std::vector<std::array<double, 3>>
IO_AERECO::ReadVecs_D(const std::string& filename)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        throw std::runtime_error("IO_AERECO::ReadVecs_D error. Can't open file: " + filename);
    }

    int n1, n2, num_points;
    fin >> n1 >> n2;
    num_points = n1 * n2;

    std::vector<std::array<double, 3>> points_for_field(num_points);
    for (int i = 0; i < num_points; i++) {
        fin >> points_for_field[i][0] >> points_for_field[i][1] >> points_for_field[i][2];
    }
    fin.close();
    return points_for_field;
}



std::vector<std::array<std::complex<double>, 3>>
IO_AERECO::ReadVecs_C(const std::string& filename)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        throw std::runtime_error("IO_AERECO::ReadVecs_C error. Can't open file: " + filename);
    }

    int n1, n2, num_points;
    fin >> n1 >> n2;
    num_points = n1 * n2;

    std::vector<std::array<std::complex<double>, 3>> points_for_field(num_points);
    for (int i = 0; i < num_points; i++) {
        fin >> points_for_field[i][0] >> points_for_field[i][1] >> points_for_field[i][2];
    }
    fin.close();
    return points_for_field;
}



std::vector<std::array<std::complex<double>, 3>> 
IO_AERECO::ReadVecs_C(const std::string& filenameRe, const std::string& filenameIm)
{
    std::ifstream finRe(filenameRe), finIm(filenameIm);
    if (!finRe.is_open() || !finIm.is_open()) {
        throw std::runtime_error("IO_AERECO::ReadVecs_C error. Can't open files: " + filenameRe + ", " + filenameIm);
    }

    int n1Re, n2Re, n1Im, n2Im, num_points;
    finRe >> n1Re >> n2Re;
    finIm >> n1Im >> n2Im;
    if (n1Re * n2Re != n1Im * n2Im) {
        throw std::runtime_error("Points file mismatch");
    } else {
        num_points = n1Re * n2Re;
    }

    std::vector<std::array<std::complex<double>, 3>>  points_for_field(num_points);
    double r1, r2, r3, im1, im2, im3;
    for (int i = 0; i < num_points; i++) {
        finRe >> r1 >> r2 >> r3;
        finIm >> im1 >> im2 >> im3;
        points_for_field[i][0] = std::complex<double>(r1, im1);
        points_for_field[i][1] = std::complex<double>(r2, im2);
        points_for_field[i][2] = std::complex<double>(r3, im3);
    }

    finRe.close();
    finIm.close();
    return points_for_field;
}





//===========================================================
//-----------Writing data to Aereco-type files---------------
//===========================================================
void IO_AERECO::WriteVecs_C(const std::string& filenameReal,
                            const std::string& filenameImage,
            const std::vector<std::array<std::complex<double>, 3>> &vals)
{
    // Creating files for recording currents
    std::ofstream fout_j_real(filenameReal);
    std::ofstream fout_j_image(filenameImage);

    if (!fout_j_real.is_open() || !fout_j_image.is_open()) {
        throw std::runtime_error("IO_AERECO::WriteVecs_C error. Can't open files:"
                                            + filenameReal + ", " + filenameImage);
    }

    int num_frm = vals.size();
    fout_j_real << 2 << " " << num_frm / 2 << std::endl;
    fout_j_image << 2 << " " << num_frm / 2 << std::endl;

    // Calculation and recording of currents
    for (int i = 0; i < num_frm; i++) {
        fout_j_real << vals[i][0].real() << " " << vals[i][1].real() << " "
                    << vals[i][2].real() << std::endl;
        fout_j_image << vals[i][0].imag() << " " << vals[i][1].imag() << " "
                     << vals[i][2].imag() << std::endl;
    }

    // Closing file streams
    fout_j_real.close();
    fout_j_image.close();
}





void IO_AERECO::WriteField_C(const std::string& filenameReal,
        const std::string& filenameImage, const std::string& filenameAbs,
        const std::vector<std::array<std::complex<double>, 3>> &field_val)
{
    std::ofstream fout_u_real(filenameReal);
    std::ofstream fout_u_image(filenameImage);
    std::ofstream fout_u_abs(filenameAbs);

    if (!fout_u_abs.is_open() ||
        !fout_u_image.is_open() ||
        !fout_u_real.is_open())
    {
        throw std::runtime_error("IO_AERECO::WriteField_C error. Can't open files:"
                        + filenameReal + ", " + filenameImage + ", " + filenameAbs);
    }

    int n1 = field_val.size() / 2;

    fout_u_real << n1 << " " << n1 << std::endl;
    fout_u_image << n1 << " " << n1 << std::endl;
    fout_u_abs << n1 << " " << n1 << std::endl;

    int n = n1 * n1;
    for (int i = 0; i < n; i++) {
        fout_u_real << field_val[i][0].real() << " " <<
                       field_val[i][1].real() << " " <<
                       field_val[i][2].real() << std::endl;
        fout_u_image << field_val[i][0].imag() << " " <<
                        field_val[i][1].imag() << " " <<
                        field_val[i][2].imag() << std::endl;
        fout_u_abs << vec_length(field_val[i].data()) << std::endl;
    }

    // Closing file streams
    fout_u_real.close();
    fout_u_image.close();
    fout_u_abs.close();
}


void IO_AERECO::WriteTauSurf2(const std::vector<std::array<std::array<double, 3>, 2>>& tau,
                                                        const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("IO_AERECO::WriteTauSurf2 error. Can't open file:" + filename);
    }
    file << std::scientific << std::setprecision(10);
    
    // Заголовок
    file << "# i  tau[i][0].x  tau[i][0].y  tau[i][0].z  tau[i][1].x  tau[i][1].y  tau[i][1].z\n";
    
    // Данные
    for (size_t i = 0; i < tau.size(); ++i) {
        file << i << "  ";
        file << tau[i][0][0] << "  " << tau[i][0][1] << "  " << tau[i][0][2] << "  ";
        file << tau[i][1][0] << "  " << tau[i][1][1] << "  " << tau[i][1][2] << "\n";
    }
    file.close();
}



}       // namespace bielcc




