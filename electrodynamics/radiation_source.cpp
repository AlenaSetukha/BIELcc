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

#include "radiation_source.h"
#include "constants.h"

    

namespace bielcc {
//===========================================================
//------------------Radiation Source-------------------------
//===========================================================
const char* toString(RadiationSourceType t) {
    switch (t) {
        case RadiationSourceType::LumpedPort:   return "LumpedPort";
        case RadiationSourceType::PlaneWave:return "PlaneWave";
        case RadiationSourceType::PointSource:  return "PointSource";
        case RadiationSourceType::SegmentSource:  return "SegmentSource";
    }
    return "Unknown";
}

const char* toString(PolarizationType t) {
    switch (t) {
        case PolarizationType::HH: return "HH";
        case PolarizationType::VV: return "VV";
    }
    return "Unknown";
}


RadiationSourceType ParseRadiationSourceType(const std::string& s)
{
    if (s == "Plane Wave") {
        return RadiationSourceType::PlaneWave;
    } else if (s == "Point Source") {
        return RadiationSourceType::PointSource;
    } else if (s == "Segment Source") {
        return RadiationSourceType::SegmentSource;
    } else if (s == "LumpedPort") {
        return RadiationSourceType::LumpedPort;
    } else {
        throw std::runtime_error("ParseRadiationSourceType error: Unknown Radiation Source Type");
    }
}


PolarizationType ParsePolarizationType(const std::string& s)
{
    if (s == "HH") {
        return PolarizationType::HH;
    } else if (s == "VV") {
        return PolarizationType::VV;
    } else {
        throw std::runtime_error("ParsePolarizationType error: Unknown Polarization Type");
    }
}




RadiationSource::RadiationSource(const RadiationSource& obj)
{
    _source_type = obj._source_type; 
    _e0_orts = obj._e0_orts;

    for (int i = 0; i < 3; i++) {
        _k_vec_ort[i] = obj._k_vec_ort[i];
    }


    _phi = obj._phi;
    _theta = obj._theta;
    _polar_types = obj._polar_types;


    _k = obj._k;
    _gamma = obj._gamma;
    _lambda = obj._lambda;
    _omega0 = obj._omega0;
}









RadiationSource::RadiationSource(RadiationSourceType type,
                                  const double* k_vec_ort,
               std::vector<std::array<double, 3>> e0_orts)
{
    for (int i = 0; i < 3; i++) {
        _k_vec_ort[i] = k_vec_ort[i];
    }

    for (int i = 0; i < e0_orts.size(); i++) {
        _e0_orts.push_back(e0_orts[i]);
    }
}




RadiationSource::RadiationSource(RadiationSourceType type,
            const double* k_vec_ort, const double* e0_ort)
{
    for (int i = 0; i < 3; i++) {
        _k_vec_ort[i] = k_vec_ort[i];
    }

    std::array<double, 3> e0 = {e0_ort[0], e0_ort[1], e0_ort[2]};
    _e0_orts.push_back(e0);
}


/*
RadiationSource::RadiationSource(RadiationSourceType type,
                                 double phi, double theta, 
                   std::vector<PolarizationType> pl_types)
{
    _source_type = type; 
    _theta = theta;
    _phi = phi;
    _polar_types = pl_types;

    _k_vec_ort[0] = - cos(_phi) * sin(_theta);
    _k_vec_ort[1] = - sin(_phi) * sin(_theta);
    _k_vec_ort[2] = - cos(_theta);

    

    double po[3], oz[3] = {0., 0., 1.};

    for (int i = 0; i < pl_types.size(); i++) {
        if (pl_types[i] == PolarizationType::HH) {
            
        } else {

        }
    }
}*/






void RadiationSource::CalcK(double freq, MaterialParam& material)
{
    _omega0 = 2. * M_PI * freq;
    material.Calc(freq);
    std::variant<double, std::complex<double>> eps = material.GetEps();
    std::variant<double, std::complex<double>> mu = material.GetMu();

    if (std::holds_alternative<double>(eps) && std::holds_alternative<double>(mu)) {
        _gamma = std::complex<double>(0., 1.) * _omega0 *
                 std::sqrt(std::get<double>(eps) * std::get<double>(mu));
    } else {
        _gamma = std::complex<double>(0., 1.) * _omega0 *
                 std::sqrt(std::get<std::complex<double>>(eps) *
                 std::get<std::complex<double>>(mu));      
    }

    _k = _gamma.imag();
    _k_vec[0] = _k * _k_vec_ort[0];
    _k_vec[1] = _k * _k_vec_ort[1];
    _k_vec[2] = _k * _k_vec_ort[2];
    _lambda = 2. * M_PI / _k;    
}
} // namespace bielcc
