#ifndef _NUM_PARAM_H_
#define _NUM_PARAM_H_

#include <iostream>

namespace bielcc {
//===========================================================
//--------Numerical parameters of the problem----------------
//===========================================================
/**
    * eps - integral calculation accuracy
    * rs - smoothing radius for the cell integral kernel function
    * rs_seg - smoothing radius for the segment integral kernel function
    * analytic_calc - dist to a highlighting singularity (relative do cell diam)
    * n_start - starting partition on the cell
    * n_start_seg - starting partition on the segment
    * p_max - limit partition on the cell (2^(p_max) cells)
    * p_max_seg - limit partition on the segment (2^(p_max_seg) segments)
    *
    * 
    * k - number of time steps (optional)
    * T - time interval (optional)
    * dt - time step (optional)
    * kappa - problem parameter (optional)
    * M - problem parameter (optional)
    * 
    * ==========================================================================
    * The smoothing radius can be specified as a fraction of each cell / grid
    * step,  depending on applicability.
 */

 /**
  * @brief Global numerical parameters for certain problem.
  * 
  */
struct NumParam {
    double eps;
    double rs, rs_seg;
    double analytic_calc;
    int n_start, n_start_seg;
    int p_max, p_max_seg;
    
    int k;
    double T, dt;
    double kappa, M;

    NumParam() = default;
    /**
        * @brief Numerical parameters structure
        * @param _eps integrals calculation accuracy
        * @param _rs smoothing radius for the cell integral (relative to cell diam)
        * @param _rs_seg smoothing radius for the segment integral (relative to cell diam)
        * @param _analytic_calc dist to a highlighting singularity (relative do cell diam)
        * @param _n_start starting partition on the cell
        * @param _n_start_seg starting partition on the segment
        * @param _p_max limit partition on the cell (2^(p_max) cells)
        * @param _p_max_seg limit partition on the segment (2^(p_max_seg) segments)
    */
    NumParam(double _eps, double _rs, double _rs_seg, double _analytic_calc,
                 int _n_start, int _n_start_seg, int _p_max, int _p_max_seg);
    NumParam(const NumParam& obj);
    ~NumParam() = default;
};
}       // namespace bielcc
#endif  // _NUM_PARAM_H_