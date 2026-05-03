#ifndef _IO_AERECO_H_
#define _IO_AERECO_H_

#include <iostream>
#include <complex>
#include <vector>
#include <array>
#include <string>

#include "element_geom.h"


/**
 * Functions:
        * write_j_toFiles - write complex cell currents to the files "_real.gv", "_image.gv"
        * write_j_toFile - write real cell currents to the file ".gv"
        * write_rkt_toFile - write collocation points to the file "body.gr"
        * write_field_toFiles - write complex cell field values ​​to the files "_real.gv", "_image.gv", "_abs.gdr"
        * write_ED_Params_toFile - write problem ED parameters to the file "ed_param.res"
 */




namespace bielcc {
/**
    * @brief Readers and writers for aereco-type files
*/
struct IO_AERECO {
    IO_AERECO() = default;

    //===========================================================
    //----------Reading data from Aereco-type files--------------
    //===========================================================
    /**
        * @brief Reading list of double 3D-vectors from aereco-file.
        * @param filename name of aereco-file 
        * @result double 3d-vectors
    */
    std::vector<std::array<double, 3>> ReadVecs_D(const std::string& filename);
    /**
        * @brief Reading list of complex 3D-vectors from aereco-file.
        * @param filename name of aereco-file 
        * @result complex 3d-vectors
    */
    std::vector<std::array<std::complex<double>, 3>> ReadVecs_C(const std::string& filename);
    /**
        * @brief Reading list of complex 3D-vectors from aereco-files.
        * @param filename name of aereco-file 
        * @result double 3d-vectors
    */
    std::vector<std::array<std::complex<double>, 3>> 
                                    ReadVecs_C(const std::string& filenameRe,
                                               const std::string& filenameIm);

    





    //===========================================================
    //-----------Writing data to Aereco-type files---------------
    //===========================================================
    /**
        * @brief Writing list of 3D-vectors to an aereco-file.
        * @param filename name of vtk-file 
        * @param vals writing data
    */
    template<typename ValType>
    void WriteVecs(const std::string& filename,
                   const std::vector<std::array<ValType, 3>> &vals)
    {
        std::ofstream fout_j(filename);
        if (!fout_j.is_open()) {
            throw std::runtime_error("IO_AERECO::WriteVecs error. Can't open file: " + filename);
        }

        int num_frm = vals.size();
        fout_j << 2 << " " << num_frm / 2 << std::endl;
        for (int i = 0; i < num_frm; i++) {
            fout_j << vals[i][0] << " " << vals[i][1] << " " << vals[i][2] << std::endl;
        }
        fout_j.close();
    }

    /**
        * @brief Writing list of 3D-vectors to an aereco-file.
        * @param filename name of vtk-file 
        * @param vals writing data 
    */
    template<typename ValType>
    void WriteVecs(const std::string& filename, int N,
                                  const ValType** vals)
    {
        std::ofstream fout_j(filename);
        if (!fout_j.is_open()) {
            throw std::runtime_error("IO_AERECO::WriteVecs error. Can't open file: " + filename);
        }

        fout_j << 2 << " " << N / 2 << std::endl;
        for (int i = 0; i < N; i++) {
            fout_j << vals[i][0] << " " << vals[i][1] << " " << vals[i][2] << std::endl;
        }
        fout_j.close();
    }

    /**
        * @brief Writing list of complex 3D-vectors to an aereco-files.
        * @param filenameReal name of real-part file 
        * @param filenameImage name of image-part file 
        * @param vals writing complex data
    */
    void WriteVecs_C(const std::string& filenameReal,
                     const std::string& filenameImage,
            const std::vector<std::array<std::complex<double>, 3>> &vals);





    /**
        * @brief Writing collocation points to an aereco-file.
        * @param filename name of file
        * @param cells cell coordinates
    */
    template<size_t CellPoints>
    void WriteRkt(const std::string& filename,
            const std::vector<double[CellPoints][3]> &cells)
    {
        std::ofstream fout(filename);
        if (!fout.is_open()) {
            throw std::runtime_error("IO_AERECO::WriteRkt error. Can't open file: " + filename);
        }

        int num_frm = cells.size();
        double cm[3];
        fout << 2 << " " << num_frm / 2 << std::endl;

        for (int i = 0; i < num_frm; i++) {
            get_center_mass(cells[i], cm);
            fout << cm[0] << " " << cm[1] << " " << cm[2] << std::endl;
        }
        fout.close();
    }


    /**
        * @brief Writing field values to an aereco-file.
        * @param filename name of file
        * @param field_val values of field vectors
    */
    template<typename FType>
    void WriteField(const std::string& filename,
                    const std::vector<std::array<FType, 3>> &field_val)
    {
        std::ofstream fout(filename);
        if (!fout.is_open()) {
            throw std::runtime_error("IO_AERECO::WriteField error. Can't open file: " +
                                                                              filename);
        }

        int n1 = field_val.size() / 2;
        fout << n1 << " " << n1 << std::endl;

        int n = n1 * n1;
        for (int i = 0; i < n; i++) {
            fout << field_val[i][0] << " " <<
                    field_val[i][1] << " " <<
                    field_val[i][2] << std::endl;
        }

        // Closing file streams
        fout.close();
    }


    /**
        * @brief Writing complex-field values to an aereco-files.
        * @param filenameReal name of file for real part
        * @param filenameImage name of file for image part
        * @param filenameAbs name of file for abs values
        * @param field_val values of field vectors
    */
    void WriteField_C(const std::string& filenameReal,
        const std::string& filenameImage, const std::string& filenameAbs,
        const std::vector<std::array<std::complex<double>, 3>> &field_val);



    /**
        * @brief Writing surface basis [N][2][3] to an aereco-file.
        * @param tau surface basis coordinates
        * @param filename name of file
    */
    void WriteTauSurf2(const std::vector<std::array<std::array<double, 3>, 2>>& tau,
                                                        const std::string& filename);



    ~IO_AERECO() = default;
};

}           // namespace bielcc
#endif      // _IO_AERECO_H_

