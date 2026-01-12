#ifndef _INTEGRAL_PARAM_H_
#define _INTEGRAL_PARAM_H_

#include <iostream>
//===========================================================
//-------------Intergartion parameters class-----------------
//===========================================================
namespace bielcc {
/** 
     * @brief Class of numerical parameters for integration.
 */
class IntegralParam {
private:
    int _idim;
    int _NStart, _PMax;
    double _epsAccuracy;


public:
    IntegralParam(): _idim(1), _NStart(8), _PMax(1), _epsAccuracy(0.0) {}
    /** 
         * @brief Class of numerical parameters for integration.
         * @param idim Dimension of the integrand
         * @param NStart Initial partition of the area/segment.
         * @param PMax Max partitions number (2^PMax)
         * @param epsAccuracy Accuracy of integral calculation
     */
    IntegralParam(int idim, int NStart, int PMax, double epsAccuracy);
    IntegralParam(const IntegralParam& obj);
    ~IntegralParam() = default;

    int GetIDim() const {return _idim;}
    int GetNStart() const {return _NStart;}
    int GetPMax() const {return _PMax;}
    double GetEpsAccur() const {return _epsAccuracy;}


    void SetNStart(int n_start) {_NStart = n_start;};
    void SetPMax(int p_max) {_PMax = p_max;}
    void SetEpsAccur(double acc) {_epsAccuracy = acc;}
};

}       //namespace bielcc
#endif  // _INTEGRAL_PARAM_H_

