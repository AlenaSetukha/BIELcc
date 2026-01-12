#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include <iostream>
#include <complex>
#include <string>
#include <variant>
#include <optional>
#include <vector>
#include <array>
#include <map>

#include "constants.h"

//===========================================================
//----------Electrodynamic parameters of materials-----------
//===========================================================
/**
 * MaterialType fields:
 *      _material_name - unique name of certain material
 *      _type - considered type of boundary condition
 *      _eps_r - relative permittivity
 *      _mu_r  - relative permeability
 *      _sigma - electrical conductivity (set zero for PERFECT conductor/dielectric)
 *      _change_type - behavioral type of material
 *      _n_refract - refractive index (optional for Dispersive materials)
 *      _k_reflect - absorption index (optional for Dispersive materials)
 *      
 *      _tand  - electric loss-tangent (dielectric absorption)
 *      _tanm  - magnetic loss-tangent (magnetization reversal absorption)
 *      _eps - absolute permittivity (in general depends on the frequency)
 *      _mu - absolute permeability (in general depends on the frequency)
 *      _skin_depth - skin layer thickness (optional for conductive materials)
 * 
 * 
 * 
 *  Notes:
 *      Assumed: E = e^{-i omega t} e^{ikr}
 * 
 *      In complex case:
 *              eps_full = eps' + i * eps''
 *              eps' = eps_r * eps0,    eps'' = sigma / omega
 *              tand = eps'' / eps' = sigma / (omega * eps_r * eps_0)
 * 
 * 
 *      Isotropic - independence of parameters from direction
 *      NonDispersive - independence of parameters from frequency
 */


namespace bielcc {

enum class ChangeMaterialTypeDir {
    Isotropic,
    Anisotropic
};
const char* toString(ChangeMaterialTypeDir t);
ChangeMaterialTypeDir ParseChangeTypeDir(const std::string& s);

enum class ChangeMaterialTypeFreq {
    NonDispersive,
    Dispersive,
    DispersiveByRefractIndex
};
const char* toString(ChangeMaterialTypeFreq t);
ChangeMaterialTypeFreq ParseChangeTypeFreq(const std::string& s);




/** 
    * @brief Material parameters.
*/
class MaterialParam {
protected:
    std::string _material_name;
    double _eps_r, _mu_r;
    double _sigma;
    double _tand, _tanm; // tanm not used yet
    ChangeMaterialTypeDir _change_type_dir = ChangeMaterialTypeDir::Isotropic;
    ChangeMaterialTypeFreq _change_type_freq = ChangeMaterialTypeFreq::Dispersive;

    std::optional<std::map<double, double>> _n_refract = std::nullopt;
    std::optional<std::map<double, double>> _k_reflect = std::nullopt; // <freq, n>, <freq, k>
    


    std::variant<double, std::complex<double>> _eps;
    std::variant<double, std::complex<double>> _mu;
    std::optional<double> _skin_depth = std::nullopt;


    void CalcEps(double freq);
    void CalcEps(double freq, double n);
    void CalcEps(double freq, double n, double k_reflect);
    void CalcMu(double freq);
    void CalcSkinDepth(double freq);

public:
    MaterialParam(): _material_name(""), _eps_r(0.), _mu_r(0.),
                            _sigma(0.), _tand(0.), _tanm(0.) {}

    MaterialParam(std::string name, ChangeMaterialTypeDir change_type_dir,
                    ChangeMaterialTypeFreq change_type_freq, double eps_r,
                     double mu_r, double sigma, double tand, double tanm):
                   _material_name(name), _change_type_dir(change_type_dir),
                    _change_type_freq(change_type_freq), _eps_r(eps_r),
                    _mu_r(mu_r), _sigma(sigma), _tand(tand), _tanm(tanm) {}



    MaterialParam(std::string name, ChangeMaterialTypeDir change_type_dir,
                                  ChangeMaterialTypeFreq change_type_freq,
                                const std::map<double, double>& n_refract,
                                       double mu_r): _material_name(name),
                                        _change_type_dir(change_type_dir),
                                      _change_type_freq(change_type_freq),
                                     _n_refract(n_refract), _mu_r(mu_r) {}

    MaterialParam(std::string name, ChangeMaterialTypeDir change_type_dir,
                                  ChangeMaterialTypeFreq change_type_freq,
                                const std::map<double, double>& n_refract,
                                const std::map<double, double>& k_reflect,
                                       double mu_r): _material_name(name),
                                        _change_type_dir(change_type_dir),
                                      _change_type_freq(change_type_freq),
                             _n_refract(n_refract), _k_reflect(k_reflect),
                                                            _mu_r(mu_r) {}
    

    MaterialParam(const MaterialParam& obj);
    ~MaterialParam() = default;


    /** 
        * @brief Material parameters calculation depending on frequency.
        * @note Calculation eps, mu, skin_depth.
        * @param freq current frequency
    */
    void Calc(double freq);



    std::variant<double, std::complex<double>> GetEps() const {
        return _eps;
    } 
    std::variant<double, std::complex<double>> GetMu() const {
        return _mu;
    }

    std::optional<double> GetSkinDepth() const {
        if (_skin_depth) {
            return _skin_depth.value();
        } else {
            return std::nullopt;
        }
    }

    std::string GetMaterialName() const {return _material_name;}
};





/** 
    * @brief Shell element parameters with twi-sided impedance condtioin.
*/

class ShellParam {
private:
    double _thickness;
    std::complex<double> _Z_Matrix[2][2], _Y_Matrix[2][2];
    MaterialParam _ext_material, _shell_material;

public:
/** 
    * @brief Shell element electrodynamic parameters.
    * @param param_Ext external material parameters
    * @param param_Shell shell material parameters
    * @param thickness shell layer thickness
*/
    ShellParam(const MaterialParam& param_Ext, const MaterialParam& param_Shell,
                  double thickness): _thickness(thickness), _ext_material(param_Ext),
                  _shell_material(param_Shell) {}
    ShellParam(const ShellParam& obj);
    ~ShellParam() = default;

    void Calc_ZYMatrix(const double freq);
    void Get_ZYMatrix(std::complex<double> (&matrZ)[2][2], std::complex<double> (&matrY)[2][2]);

};
}       // namespace bielcc
#endif  // _MATERIAL_H_
