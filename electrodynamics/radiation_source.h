#ifndef _RADIATION_SOURCE_H_
#define _RADIATION_SOURCE_H_

#include <iostream>
#include <complex>
#include <string>
#include <variant>
#include <limits>
#include <optional>
#include <vector>
#include <array>
#include <map>

#include "constants.h"
#include "material.h"

//===========================================================
//--------------Radiation source parameters------------------
//===========================================================
/**
 * RadiationSource fields:
 *      _phi - polar angle [0, 180] in degrees
 *      _teta - azimuth angle [0, 360] in degrees
 *      _source_type - type of the e/d source
 *      _polar_types - set of required polarizations for calculation
 *      _orts - ort-vectors E according to polar types list
 * 
 *      _k - wave number
 *      _k_vec - wave vector 
 *      _gamma - complex propogation constant
 *      _lambda - wavelength
 *      _omega0 - circular frequency
 * 
 * 
 * 
 *  Notes:
 *      Assumed: E = e^{-i omega t} e^{ikr}
 *      Dont forget to convert into radians:  * M_PI / 180.
 * 
 * 
 *  HH = s (это наша вертикальная)
 *  VV = p (это наша горизонтальная)
 */


namespace bielcc {

enum class RadiationSourceType {
    PlaneWave,
    PointSource,
    SegmentSource,
    LumpedPort
};
const char* toString(RadiationSourceType t);
RadiationSourceType ParseRadiationSourceType(const std::string& s);



enum class PolarizationType {
    HH,
    VV
};
const char* toString(PolarizationType t);
PolarizationType ParsePolarizationType(const std::string& s);




/** 
    * @brief Description for radiation source properties.
*/
class RadiationSource {
private:
    RadiationSourceType _source_type;
    std::vector<std::array<double, 3>> _e0_orts;
    double _k_vec_ort[3];


    // variant for Plane Wave
    double _phi, _theta;
    std::vector<PolarizationType> _polar_types; // {HH}, {VV} or {HH, VV}






    double _k, _k_vec[3];
    std::complex<double> _gamma;
    double _lambda, _omega0;

public:
    /** 
        * @brief Class for describing the radiation source.
        * @param type radiation source type
        * @param k_vec_ort ort of wave-vector k
        * @param e0_orts ort for E0
    */
    RadiationSource(RadiationSourceType type, const double* k_vec_ort,
                           std::vector<std::array<double, 3>> e0_orts);
    RadiationSource(RadiationSourceType type, const double* k_vec_ort,
                                                 const double* e0_ort);
                                             



    /** 
        * @brief Class for describing the radiation source.
        * @param type radiation source type
        * @param phi polar angle [0, 180] in degrees
        * @param theta azimuth angle [0, 360] in degrees
        * @param pl_types polarization types
        * @note The opposite direction must be given. Example:
        * for k = (-1, 0, 0) {phi = 0, theta = 90}. 
    */
    /*
    RadiationSource(RadiationSourceType type, double phi, double theta, 
                                std::vector<PolarizationType> pl_types);*/



    RadiationSource(const RadiationSource& obj);
    ~RadiationSource() = default;


    /** 
        * @brief Wave number and wave vector calculation in a certain medium.
        * @param freq calculation frequency
        * @param material material properties
    */ 
    void CalcK(double freq, MaterialParam& material);


    RadiationSourceType GetSourceType() const {return _source_type;}  

    double GetPhi() const {return _phi;}
    double GetTeta() const {return _theta;}
    
    int GetNPolarTypes() const {return _polar_types.size();}
    PolarizationType GetPolaryzationType(int i) const {return _polar_types[i];}
    void GetOrt(int i, double* ort) const { ort[0] = _e0_orts[i][0];
                                            ort[1] = _e0_orts[i][1];
                                            ort[2] = _e0_orts[i][2];} 



    /** 
        * @brief Wave number in a certain medium.
        * @note Call after Calc()
    */
    double GetK() const {return _k;}; 
    /** 
        * @brief Wave vector in a certain medium.
        * @param k_vec result vector
        * @note Call after Calc()
    */ 
    void GetKVec(double* k_vec) const {k_vec[0] = _k_vec[0];
                                       k_vec[1] = _k_vec[1];
                                       k_vec[2] = _k_vec[2];};
    /** 
        * @brief Complex propogation constant in a certain medium.
        * @note Call after Calc()
    */
    std::complex<double> GetGamma() const {return _gamma;}
};

}       // namespace bielcc
#endif  // _RADIATION_SOURCE_H
