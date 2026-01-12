#ifndef _IO_VTK_H_
#define _IO_VTK_H_

#include <iostream>
#include <complex>
#include <vector>
#include <array>
#include <string>


namespace bielcc {
/**
    * @brief Readers and writers for .vtk files
*/
struct IO_VTK {
    IO_VTK() = default;

    //===========================================================
    //------------Reading data from vtk-files--------------------
    //===========================================================
    /**
        * @brief Reading points coordinates from file
        * @param filename name of file
        * @result list of 3D-vectors.
    */
    std::vector<std::array<double, 3>> ReadVecs(const std::string& filename);






    //===========================================================
    //--------------Writing data to vtk-files--------------------
    //===========================================================
    /**
        * @brief Writing list of 3D-vectors to vtk-file.
        * @param filename name of vtk-file 
        * @param vals writing data
    */
    void WriteVecs(const std::string& filename, 
                    const std::vector<std::array<double, 3>>& vals);

    /**
        * @brief Writing list of N 3D-vectors to vtk-file.
        * @param filename name of vtk-file
        * @param vals writing data
    */
    void WriteVecs(const std::string& filename, int N, const double** vals);


    /**
        * @brief Writing points and double field-vals in vtk-file.
        * @param filename name of vtk-file (path)
        * @param points list of points
        * @param field_vectors field-vals (3D-vectors) in points
    */
    void WritePntsWithField_D(const std::string& filename,
                    const std::vector<std::array<double, 3>>& points,
                    const std::vector<std::array<double, 3>>& field_vectors);

    /**
        * @brief Writing points and complex field-vals in vtk-file.
        * @param filename name of vtk-file (path)
        * @param points list of points
        * @param field_vectors field-vals (3D-vectors) in points
    */
    void WritePntsWithField_C(const std::string& filename,
            const std::vector<std::array<double, 3>>& points,
            const std::vector<std::array<std::complex<double>, 3>>& field_vectors);


    ~IO_VTK() = default;
};

}           // namespace bielcc
#endif      // _IO_VTK_H_

