#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <vector>
#include <array>
#include <string>

#include "io_vtk.h"


namespace bielcc {
//===========================================================
//-----------Reading data from .vtk-type files---------------
//===========================================================
std::vector<std::array<double, 3>> IO_VTK::ReadVecs(const std::string& filename)
{
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("ReadVecs .vtk. Cannot open file: " + filename);
    }
    
    std::string line;
    std::vector<std::array<double, 3>> points;
    int line_num = 0;
    
    try {
        // Пропустить заголовок (первые 3 строки)
        for (int i = 0; i < 3; ++i) {
            if (!std::getline(file, line)) {
                throw std::runtime_error("File is too short (missing header)");
            }
            line_num++;
        }
        
        // Прочитать DATASET строку (пропускаем пустые строки и комментарии)
        do {
            if (!std::getline(file, line)) {
                throw std::runtime_error("Failed to read DATASET line");
            }
            line_num++;
        } while (line.empty() || line[0] == '#');
        
        if (line.find("DATASET") == std::string::npos) {
            throw std::runtime_error("Expected DATASET keyword at line " + 
                                   std::to_string(line_num));
        }
        
        // Прочитать POINTS строку (пропускаем пустые строки)
        do {
            if (!std::getline(file, line)) {
                throw std::runtime_error("Failed to read POINTS line");
            }
            line_num++;
        } while (line.empty() || line[0] == '#');
        
        std::istringstream iss(line);
        std::string keyword, data_type;
        int num_points;
        
        iss >> keyword >> num_points >> data_type;
        
        if (keyword != "POINTS") {
            throw std::runtime_error("Expected POINTS keyword at line " + 
                                   std::to_string(line_num));
        }
        
        if (num_points <= 0) {
            throw std::runtime_error("Invalid number of points: " + 
                                   std::to_string(num_points));
        }
        
        // Прочитать координаты точек
        for (int i = 0; i < num_points; ++i) {
            // Пропускаем пустые строки и комментарии
            do {
                if (!std::getline(file, line)) {
                    throw std::runtime_error("Unexpected end of file at point " + 
                                           std::to_string(i));
                }
                line_num++;
            } while (line.empty() || line[0] == '#');
            
            std::istringstream coord_stream(line);
            double x, y, z;
            
            if (!(coord_stream >> x >> y >> z)) {
                throw std::runtime_error("Failed to parse coordinates at line " + 
                                       std::to_string(line_num) + 
                                       " (point " + std::to_string(i) + ")");
            }
            
            points.push_back({x, y, z});
        }
        
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Error reading VTK file: ") + e.what());
    }
    
    file.close();
    return points;
}






//===========================================================
//--------------Writing data to vtk-files--------------------
//===========================================================
void IO_VTK::WriteVecs(const std::string& filename,
              const std::vector<std::array<double, 3>>& vals)
{
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("WriteVec .vtk error. Cannot open file: " + filename);
    }
    
    // Заголовок VTK
    file << "# vtk DataFile Version 3.0\n";
    file << "Points on surface\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    
    // Запись точек
    file << "POINTS " << vals.size() << " double\n";
    for (const auto& pt : vals) {
        file << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }
    
    file.close();
}


void IO_VTK::WriteVecs(const std::string& filename, int N,
                                      const double** vals)
{
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("WriteVec .vtk error. Cannot open file: " + filename);
    }
    
    // Заголовок VTK
    file << "# vtk DataFile Version 3.0\n";
    file << "Points on surface\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    
    // Запись точек
    file << "POINTS " << N << " double\n";
        for (int i = 0; i < N; i++) {
        file << vals[i][0] << " " << vals[i][1] << " " << vals[i][2] << "\n";
    }
    
    file.close();
}








void IO_VTK::WritePntsWithField_D(const std::string& filename,
                    const std::vector<std::array<double, 3>>& points,
                    const std::vector<std::array<double, 3>>& field_vectors)
{
    if (points.size() != field_vectors.size()) {
        throw std::runtime_error("Points and field vectors must have same size");
    }
    
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    // Заголовок
    file << "# vtk DataFile Version 3.0\n";
    file << "Points with vector field\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    
    // Запись точек
    file << "POINTS " << points.size() << " double\n";
    for (const auto& pt : points) {
        file << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }
    
    // Запись векторного поля (например, ток)
    file << "POINT_DATA " << field_vectors.size() << "\n";
    file << "VECTORS current double\n";
    for (const auto& vec : field_vectors) {
        file << vec[0] << " " << vec[1] << " " << vec[2] << "\n";
    }
    
    file.close();
}

void IO_VTK::WritePntsWithField_C(const std::string& filename,
            const std::vector<std::array<double, 3>>& points,
            const std::vector<std::array<std::complex<double>, 3>>& field_vectors)
{
    if (points.size() != field_vectors.size()) {
        throw std::runtime_error("Points and field vectors must have same size");
    }
    
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    // Заголовок
    file << "# vtk DataFile Version 3.0\n";
    file << "Points with complex vector field\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    
    // Запись точек
    file << "POINTS " << points.size() << " double\n";
    for (const auto& pt : points) {
        file << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }
    
    // Запись ТОЧЕК ДАННЫХ
    file << "POINT_DATA " << field_vectors.size() << "\n";
    
    // Записываем действительную часть
    file << "VECTORS real_part double\n";
    for (const auto& vec : field_vectors) {
        file << vec[0].real() << " " << vec[1].real() << " " << vec[2].real() << "\n";
    }
    
    // Записываем мнимую часть
    file << "VECTORS imaginary_part double\n";
    for (const auto& vec : field_vectors) {
        file << vec[0].imag() << " " << vec[1].imag() << " " << vec[2].imag() << "\n";
    }
    
    // Записываем модуль (амплитуду) комплексного вектора
    file << "VECTORS magnitude double\n";
    for (const auto& vec : field_vectors) {
        double mag_x = std::abs(vec[0]);
        double mag_y = std::abs(vec[1]);
        double mag_z = std::abs(vec[2]);
        file << mag_x << " " << mag_y << " " << mag_z << "\n";
    }
    
    // Записываем фазу (угол) комплексного вектора как скалярное поле
    file << "SCALARS phase_x double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& vec : field_vectors) {
        file << std::arg(vec[0]) << "\n";
    }
    
    file << "SCALARS phase_y double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& vec : field_vectors) {
        file << std::arg(vec[1]) << "\n";
    }
    
    file << "SCALARS phase_z double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& vec : field_vectors) {
        file << std::arg(vec[2]) << "\n";
    }
    
    file.close();
}


}       // namespace bielcc