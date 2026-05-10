#ifndef _NUM_PARAM_H_
#define _NUM_PARAM_H_

#include <iostream>

namespace bielcc {
//===========================================================
//--------Numerical parameters of the problem----------------
//===========================================================
/**
    * @brief Global numerical parameters for certain problem.
    * @param eps integral calculation accuracy
    * @param rs smoothing radius for the cell integral kernel function (relative to small cell)
    * @param rs_seg smoothing radius for the segment integral kernel function (relative to small cell)
    * @param analytic_calc dist to a highlighting singularity (relative do cell diam) (relative to mesh step)
    * @param n_start starting partition on the cell
    * @param n_start_seg starting partition on the segment
    * @param p_max limit partition on the cell (2^(p_max) cells)
    * @param p_max_seg limit partition on the segment (2^(p_max_seg) segments)
    * @param k number of time steps (optional)
    * @param T time interval (optional)
    * @param dt time step (optional)
    * @param kappa problem parameter (optional)
    * @param M problem parameter (optional)
    * @note The smoothing radius can be specified as a fraction of
    * each cell / cell 2-nd layer (depending on applicability)
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