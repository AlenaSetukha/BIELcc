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
        * @brief Reading 3D-vector list from vtk file
        * @param filename Name of file (full path)
        * @result List of 3D-vectors
        * @note In file should be: POINTS N double {...}
    */
    std::vector<std::array<double, 3>> ReadVecs(const std::string& filename);
    /**
        * @brief Reading 3D-vector list from vtk file
        * @param filename Name of file (full path)
        * @param vecs Double-vals vectors (result)
        * @note Memory must be allocated correctly: [num_vecs][3]
        * @note In file should be: POINTS N double {...}
    */
    void ReadVecs(const std::string& filename, double** vecs);


    //===========================================================
    //------------Reading data from vtk-files--------------------
    //===========================================================
    /**
        * @brief Reading complex 3D-vector list from vtk-file
        * @param filename Name of file (full path)
        * @result List of 3D-vectors (complex-vals)
        * @note In file should be: VECTORS real_part double {...} VECTORS image_part double {...}
    */
    std::vector<std::array<std::complex<double>, 3>> ReadVecs_Cmplx(const std::string& filename);


    /**
        * @brief Reading complex 3D-vector list from vtk-file
        * @param filename Name of file (full path)
        * @param vecs Cmplx-vals vectors (result)
        * @note Memory must be allocated correctly: [num_vecs][3]
        * @note In file should be: VECTORS real_part double {...} VECTORS image_part double {...}
    */
    void ReadVecs_Cmplx(const std::string& filename,
                                std::complex<double>** vecs);







    //===========================================================
    //--------------Writing data to vtk-files--------------------
    //===========================================================
    /**
        * @brief Writing list of 3D-vectors to vtk-file.
        * @param filename Name of file (full path)
        * @param vals Writing data
    */
    void WriteVecs(const std::string& filename, 
                    const std::vector<std::array<double, 3>>& vals);

    /**
        * @brief Writing list of N 3D-vectors to vtk-file.
        * @param filename Name of file (full path)
        * @param N Length of vals
        * @param vals Writing data [N][3]
    */
    void WriteVecs(const std::string& filename, int N, const double** vals);


    /**
        * @brief Writing points and double field-vals in vtk-file.
        * @param filename Name of file (full path)
        * @param points List of points
        * @param field_vectors Field-vals (3D-vectors) in points
    */
    void WritePntsWithField_D(const std::string& filename,
                    const std::vector<std::array<double, 3>>& points,
                    const std::vector<std::array<double, 3>>& field_vectors);


    /**
        * @brief Writing points and double field-vals in vtk-file.
        * @param filename Name of file (full path)
        * @param points List of points
        * @param NF Length of field_vectors
        * @param field_vectors Field-vals (3D-vectors) in points
    */
    void WritePntsWithField_D(const std::string& filename,
                            const std::vector<std::array<double, 3>>& points,
                                        int NF, const double** field_vectors);

    /**
        * @brief Writing points and complex field-vals in vtk-file.
        * @param filename Name of file (full path)
        * @param points List of points
        * @param field_vectors Field-vals (3D-vectors) in points
    */
    void WritePntsWithField_C(const std::string& filename,
            const std::vector<std::array<double, 3>>& points,
            const std::vector<std::array<std::complex<double>, 3>>& field_vectors);



    /**
        * @brief Writing points and complex field-vals in vtk-file.
        * @param filename Name of file (full path)
        * @param points List of points
        * @param NF Field_vectors length
        * @param field_vectors Complex field-vals (3D-vectors) in points
    */
    void WritePntsWithField_C(const std::string& filename,
                    const std::vector<std::array<double, 3>>& points,
                  int NF, const std::complex<double>** field_vectors);


    ~IO_VTK() = default;
};

}           // namespace bielcc
#endif      // _IO_VTK_H_
