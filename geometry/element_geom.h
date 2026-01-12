#ifndef _ELEMENT_GEOM_H_
#define _ELEMENT_GEOM_H_

#include <complex>
//===========================================================
//------------Library of basic geometric operations----------
//===========================================================
/**
 * Description:
 *      scal_prod - scalar product of vectors
 *      vec_prod - vector product of vectors
 *      tr_square - triangle area
 *      quadr_square - area of quadrilateral
 *      solid_angle - solid angle calculation (radians)
 *      norm_func - normal vector to the cell
 *      perp - point projection onto a segment (perpendicular vector)
 *      near_point - base of a perpendicular
 *      get_nu - normal to a segment
 *      get_diam - cell cize (max side / diagonal)
 *      get_center_mass - mass center on the cell
 *      check_points_match - coincidence of two points (with accuracy)
 * 
 *      dist - euclidean distance btw points / vectors
 *      vec_length - vector length
 */


namespace bielcc {
//==Scalar product==
double scal_prod(const double* vec_1, const double* vec_2);
std::complex<double> scal_prod(const std::complex<double>* vec_1, const double* vec_2);
std::complex<double> scal_prod(const double* vec_1, const std::complex<double>* vec_2);
std::complex<double> scal_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2);

//==Vector product==
void vec_prod(const double* vec_1, const double* vec_2, double* res);
void vec_prod(const std::complex<double>* vec_1, const double* vec_2, std::complex<double>* res);
void vec_prod(const double* vec_1, const std::complex<double>* vec_2, std::complex<double>* res);
void vec_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2, std::complex<double>* res);


//==Square of triagnle==
double tr_square(const double* pnt_1, const double* pnt_2, const double* pnt_3);
//==Quadr square==
double quadr_square(const double* a, const double* b, const double* c, const double* d);

//====Cell sqaure=====
double cell_square(const double (&rut0)[4][3]);
double cell_square(const double (&rut0)[3][3]);


//==Solid angle==
double solid_angle(const double* x_a, const double* x_b, const double* x_c, const double* x);

//==Cell normal==
void norm_func(const double (&rut0)[4][3], double* norm_res);
void norm_func(const double (&rut0)[3][3], double* norm_res);

//==Perpendicular to the segment(проекция точки на отрезок)==
void perp(const double* a, const double* b, const double* x, double* n);


//==Nearest point to the edge==
void near_point(const double (&seg)[2][3], const double* x, double* z, double& dist_res);
void near_point(const double* a, const double* b, const double* x, double* z, double& dist_res);


//==Normal to cell edge(нормаль к краю)==
void get_nu(const double (&rut0)[4][3], const double* norm, const int i, const int inext, double* nu);
void get_nu(const double* a, const double* b, const double* norm, double* nu);


//==Cell diameter==
double get_diam(const double (&root_tmp)[4][3]);
double get_diam(const double (&root_tmp)[3][3]);

//==Center of mass==
void get_center_mass(const double* a, const double* b, const double* c, const double* d, double* res);
void get_center_mass(const double* a, const double* b, const double* c, double* res);
void get_center_mass(const double (&root_tmp)[4][3], double* res);
void get_center_mass(const double (&root_tmp)[3][3], double* res);


//==Check two points are different==
int check_points_match(const double* a, const double* b);




//==Distance btw 2 points(vectors)==
template <typename T>
T dist(const T* vec_1, const T* vec_2)
{
    T res = sqrt((vec_1[0] - vec_2[0]) * (vec_1[0] - vec_2[0]) +
        (vec_1[1] - vec_2[1]) * (vec_1[1] - vec_2[1]) +
        (vec_1[2] - vec_2[2]) * (vec_1[2] - vec_2[2]));
    return res;
}


//==Vector length==
template <typename T>
double vec_length(const T* vec_1)
{
    return sqrt(std::abs(vec_1[0]) * std::abs(vec_1[0]) +
                std::abs(vec_1[1]) * std::abs(vec_1[1]) +
                std::abs(vec_1[2]) * std::abs(vec_1[2]));
}

}       // namespace bielcc
#endif  // _ELEMENT_GEOM_H_