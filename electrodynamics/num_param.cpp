#include <iostream>

#include "num_param.h"

namespace bielcc {

NumParam::NumParam(double _eps, double _rs, double _rs_seg,
                       double _analytic_calc, int _n_start,
               int _n_start_seg, int _p_max, int _p_max_seg)
{
    eps = _eps, rs = _rs, rs_seg = _rs_seg;
    analytic_calc = _analytic_calc;
    n_start = _n_start, n_start_seg = _n_start_seg;
    p_max = _p_max, p_max_seg = _p_max_seg;
    T = 0., dt = 0., kappa = 0., M = 0.;
    k = 0;
}


NumParam::NumParam(const NumParam& obj) {
    eps = obj.eps, rs = obj.rs, rs_seg = obj.rs_seg;
    n_start = obj.n_start, n_start_seg = obj.n_start_seg;
    p_max = obj.p_max, p_max_seg = obj.p_max_seg;
    T = obj.T, dt = obj.dt, kappa = obj.kappa, M = obj.M;
    k = obj.k;
}


}   // namespace bielcc

