#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iostream>
#include <string>
#include <map>
#include <variant>

#include "material.h"
#include "constants.h"

    

namespace bielcc {
//===========================================================
//--------------Material Parameters--------------------------
//===========================================================
MaterialParam::MaterialParam(const MaterialParam& obj)
{
    _material_name = obj._material_name;
    _eps_r = obj._eps_r;
    _mu_r = obj._mu_r;
    _sigma = obj._sigma;
    _change_type_dir = obj._change_type_dir;
    _change_type_freq = obj._change_type_freq;
    _eps = obj._eps;
    _mu = obj._mu;

    if (obj._n_refract) {
        _n_refract = obj._n_refract;
    } else {
        _n_refract = std::nullopt;
    }
    if (obj._k_reflect) {
        _k_reflect = obj._k_reflect;
    } else {
        _k_reflect = std::nullopt;
    }
    if (obj._skin_depth) {
        _skin_depth = obj._skin_depth;
    } else {
        _skin_depth = std::nullopt;
    }

    _tand = obj._tand;
    _tanm = obj._tanm;
    _eps = obj._eps;
    _mu = obj._mu;
}







void MaterialParam::Calc(double freq) {
    if (_change_type_freq == ChangeMaterialTypeFreq::DispersiveByRefractIndex) {
        if (_n_refract != std::nullopt && _k_reflect != std::nullopt) {
            CalcEps(freq, (*_n_refract)[freq], (*_k_reflect)[freq]); 
        } else if (_n_refract != std::nullopt) {
            CalcEps(freq, (*_n_refract)[freq]);
        } else {
            std::cerr << "MaterialParam::Calc: No matching frequency in refractive indexes" << "\n";
            return;
        }
    } else {
        CalcEps(freq);
    }
    CalcMu(freq);
    CalcSkinDepth(freq);
}







void MaterialParam::CalcEps(double freq) {
    double omega0 = 2. * M_PI * freq;

    if (_sigma == 0.) {
        _eps = ED_Constants::eps0 * _eps_r;
        _tand = 0.;
    } else {
        _eps = std::complex<double>(ED_Constants::eps0 * _eps_r, (_sigma / omega0));
        _tand = _sigma / (omega0 * ED_Constants::eps0 * _eps_r);
    }
}


void MaterialParam::CalcEps(double freq, double n) {
    _eps = n * n * ED_Constants::eps0;
    _tand =  0.;
}



void MaterialParam::CalcEps(double freq, double n, double k) {
    _eps = std::complex<double>(n * n - k * k, 2. * n * k) * ED_Constants::eps0;
    _tand =  2. * n * k * ED_Constants::eps0;
}


void MaterialParam::CalcMu(double freq) {
    /*if (_tanm == 0.) {
        _mu = Constants::mu0 * _mu_r;
    } else {
        _mu = std::complex<double>(Constants::mu0 * _mu_r, Constants::mu0 * _mu_r * _tanm);
    }*/
   _mu = ED_Constants::mu0 * _mu_r;
   _tanm = 0.;
}


void MaterialParam::CalcSkinDepth(double freq)
{
    double omega0 = 2. * M_PI * freq;
    if (_sigma != 0.) {
        _skin_depth = std::sqrt(2. / (omega0 * _mu_r * ED_Constants::mu0 * _sigma));
    } else {
        _skin_depth = std::nullopt;
    }
}
















//===========================================================
//--------------------Shell Element--------------------------
//===========================================================
ShellParam::ShellParam(const ShellParam& obj) {
    _thickness = obj._thickness;
    _ext_material = obj._ext_material;
    _shell_material = obj._shell_material;
    _Z_Matrix[0][0] = obj._Z_Matrix[0][0];
    _Z_Matrix[0][1] = obj._Z_Matrix[0][1];
    _Z_Matrix[1][0] = obj._Z_Matrix[1][0];
    _Z_Matrix[1][1] = obj._Z_Matrix[1][1];
    _Y_Matrix[0][0] = obj._Y_Matrix[0][0];
    _Y_Matrix[0][1] = obj._Y_Matrix[0][1];
    _Y_Matrix[1][0] = obj._Y_Matrix[1][0];
    _Y_Matrix[1][1] = obj._Y_Matrix[1][1];
}

void  ShellParam::Calc_ZYMatrix(const double freq)
{
    std::cout << "Calculating ZY Matrices in Shell Element is not yet written" << std::endl;
    _ext_material.Calc(freq);
    _shell_material.Calc(freq);
    double omega0 = 2. * M_PI * freq;
    std::variant<double, std::complex<double>> muExt = _ext_material.GetMu();
    std::variant<double, std::complex<double>> muShell = _shell_material.GetMu();
    std::variant<double, std::complex<double>> epsExt = _ext_material.GetEps();
    std::variant<double, std::complex<double>> epsShell = _shell_material.GetEps();
    std::variant<double, std::complex<double>> ZExt, ZShell;


/*
    std::visit([&](auto&& arg1) {
        std::visit([&](auto&& arg2) {
            // Здесь описываем логику для каждой пары
            using T1 = std::decay_t<decltype(arg1)>;
            using T2 = std::decay_t<decltype(arg2)>;

            if constexpr (std::is_same_v<T1, int> && std::is_same_v<T2, std::string>) {
                std::cout << "int и string: " << std::to_string(arg1) + arg2 << std::endl;
            } else if constexpr (std::is_same_v<T1, int> && std::is_same_v<T2, int>) {
                std::cout << "Оба int: " << arg1 * arg2 << std::endl;
            }
            // ... и так далее для остальных комбинаций
            else {
                std::cout << "Другие типы" << std::endl;
            }
        }, v2);
    }, v1);



*/

    //std::complex<double> z1 = - std::complex<double>(0., 1.) * omega0 * 
}

void  ShellParam::Get_ZYMatrix(std::complex<double> (&matrZ)[2][2], std::complex<double> (&matrY)[2][2])
{
    matrZ[0][0] = _Z_Matrix[0][0], matrZ[0][1] = _Z_Matrix[0][1];
    matrZ[1][0] = _Z_Matrix[1][0], matrZ[1][1] = _Z_Matrix[1][1];
    matrY[0][0] = _Y_Matrix[0][0], matrY[0][1] = _Y_Matrix[0][1];
    matrY[1][0] = _Y_Matrix[1][0], matrY[1][1] = _Y_Matrix[1][1];
}






















//===========================================================
//---------------Material string converters------------------
//===========================================================
const char* toString(ChangeMaterialTypeDir t) {
    switch (t) {
        case ChangeMaterialTypeDir::Isotropic:   return "Isotropic";
        case ChangeMaterialTypeDir::Anisotropic: return "Anisotropic";
    }
    return "Unknown";
}



ChangeMaterialTypeDir ParseChangeTypeDir(const std::string& s) {
    if (s == "Isotropic") {
        return ChangeMaterialTypeDir::Isotropic;
    } else if (s == "Anisotropic") {
        return ChangeMaterialTypeDir::Anisotropic;
    } else {
        throw std::runtime_error("ParseChangeTypeDir error: Unknown Change Direction Material Type");
    }
}


const char* toString(ChangeMaterialTypeFreq t) {
    switch (t) {
        case ChangeMaterialTypeFreq::Dispersive:  return "Dispersive";
        case ChangeMaterialTypeFreq::NonDispersive:  return "NonDispersive";
        case ChangeMaterialTypeFreq::DispersiveByRefractIndex:  return "DispersiveByRefractInex";
    }
    return "Unknown";
}

ChangeMaterialTypeFreq ParseChangeTypeFreq(const std::string& s) {
    if (s == "Dispersive") {
        return ChangeMaterialTypeFreq::Dispersive;
    } else if (s == "NonDispersive") {
        return ChangeMaterialTypeFreq::NonDispersive;
    } else if (s == "DispersiveByRefractIndex") {
        return ChangeMaterialTypeFreq::DispersiveByRefractIndex;
    } else {
        throw std::runtime_error("ParseChangeTypeFreq error: Unknown Cnange Frequency Material Type");
    }
}


} // namespace bielcc
