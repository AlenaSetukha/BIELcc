#include <iostream>
#include "integral_param.h"

namespace bielcc {
IntegralParam::IntegralParam(int idim, int NStart, int PMax, double epsAccuracy)
{
    _idim = idim;
    _NStart = NStart;
    _PMax = PMax;
    _epsAccuracy = epsAccuracy;
}

IntegralParam::IntegralParam(const IntegralParam& obj)
{
    _idim = obj._idim;
    _NStart = obj._NStart;
    _PMax = obj._PMax;
    _epsAccuracy = obj._epsAccuracy;
}
} // namespace bielcc