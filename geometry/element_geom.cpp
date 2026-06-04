#include <iostream>
#include <cmath>
#include <complex>
#include <limits>

#include "element_geom.h"

namespace bielcc {
//====================================Scalar product=============================================
double scal_prod(const double* vec_1, const double* vec_2) {
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

std::complex<double> scal_prod(const std::complex<double>* vec_1, const double* vec_2) {
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

std::complex<double> scal_prod(const double* vec_1, const std::complex<double>* vec_2) {
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

std::complex<double> scal_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2) {
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}



//====================================Vector product=============================================
void vec_prod(const double* vec_1, const double* vec_2, double* res) {
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
}

void vec_prod(const std::complex<double>* vec_1, const double* vec_2, std::complex<double>* res) {
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
}

void vec_prod(const double* vec_1, const std::complex<double>* vec_2, std::complex<double>* res) {
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
}

void vec_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2, std::complex<double>* res) {
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
}

    
//====================================Square of triagnle=========================================
double tr_square(const double* pnt_1, const double* pnt_2, const double* pnt_3) {
    double ab[3], ac[3];
    for (int i = 0; i < 3; i++) {
        ab[i] = pnt_2[i] - pnt_1[i];
        ac[i] = pnt_3[i] - pnt_1[i];
    }
    double n[3];
    vec_prod(ab, ac, n);
    double res = scal_prod(n, n);
    res = sqrt(res) / 2.;
    return res;
}

//===========================================Quadr square===========================================
double quadr_square(const double* a, const double* b, const double* c, const double* d) {
    double ac[3], bd[3], vec[3];
    for (int i = 0; i < 3; i++) {
        ac[i] = c[i] - a[i];
        bd[i] = d[i] - b[i];
    }
    vec_prod(ac, bd, vec);
    return vec_length(vec) / 2.;
}


double cell_square(const double (&rut0)[4][3]) {
    double ac[3], bd[3], vec[3];
    for (int i = 0; i < 3; i++) {
        ac[i] = rut0[2][i] - rut0[0][i];
        bd[i] = rut0[3][i] - rut0[1][i];
    }
    vec_prod(ac, bd, vec);
    return vec_length(vec) / 2.;
}

double cell_square(const double (&rut0)[3][3]) {
    double ab[3], ac[3];
    for (int i = 0; i < 3; i++) {
        ab[i] = rut0[1][i] - rut0[0][i];
        ac[i] = rut0[2][i] - rut0[0][i];
    }
    double n[3];
    vec_prod(ab, ac, n);
    double res = scal_prod(n, n);
    res = sqrt(res) / 2.;
    return res;
}

//=======================================Solid angle=============================================
double solid_angle(const double* x_a, const double* x_b, const double* x_c, const double* x) {
    double r1[3], r2[3], r3[3], l1, l2, l3, deg, v[3];
    r1[0] = x_a[0] - x[0];
    r1[1] = x_a[1] - x[1];
    r1[2] = x_a[2] - x[2];

    r2[0] = x_b[0] - x[0];
    r2[1] = x_b[1] - x[1];
    r2[2] = x_b[2] - x[2];

    r3[0] = x_c[0] - x[0];
    r3[1] = x_c[1] - x[1];
    r3[2] = x_c[2] - x[2];

    l1 = vec_length(r1);
    l2 = vec_length(r2);
    l3 = vec_length(r3);

    vec_prod(r3, r2, v);

    deg = scal_prod(r1, v) / (l1 * l2 * l3 + scal_prod(r1, r2) * l3 +
                    scal_prod(r2, r3) * l1 + scal_prod(r3, r1) * l2);
    return 2 * atan(deg);
}

//=======================================Cell normal=============================================
void norm_func(const double (&rut0)[4][3], double* norm_res) {
    double ac[3], bd[3];
    for (int i = 0; i < 3; i++) {
        ac[i] = rut0[2][i] - rut0[0][i];
        bd[i] = rut0[3][i] - rut0[1][i]; 
    }
    double v[3], l;
    vec_prod(ac, bd, v);
    l = vec_length(v);
    norm_res[0] = v[0] / l;
    norm_res[1] = v[1] / l;
    norm_res[2] = v[2] / l;
}

void norm_func(const double (&rut0)[3][3], double* norm_res) {
    double ac[3], ab[3];
    for (int i = 0; i < 3; i++) {
        ac[i] = rut0[2][i] - rut0[0][i];
        ab[i] = rut0[1][i] - rut0[0][i]; 
    }
    double v[3], l;
    vec_prod(ab, ac, v);
    l = vec_length(v);
    norm_res[0] = v[0] / l;
    norm_res[1] = v[1] / l;
    norm_res[2] = v[2] / l;
}



//======================Perpendicular to the segment(проекция точки на отрезок)==================
void perp(const double* a, const double* b, const double* x, double* n) {
    double ab[3], ax[3];
    for (int i = 0; i < 3; i++) {
        ab[i] = b[i] - a[i];
        ax[i] = x[i] - a[i];
    }
    double sc_pr = scal_prod(ab, ax) / scal_prod(ab, ab);
    for (int i = 0; i < 3; i++) {
        n[i] = a[i] + (b[i] - a[i]) * sc_pr;
    }  
}



//==================================Nearest point to the segment=================================
void near_point(const double (&seg)[2][3], const double* x, double* z, double& dist_res) {
    double a[3], b[3], n[3];
    double dist_ab, dist_an, dist_bn;
    for (int i = 0; i < 3; i++) {
        a[i] = seg[0][i];
        b[i] = seg[1][i];
    }
    perp(a, b, x, n);
    dist_ab = dist(a,b);
    dist_an = dist(a,n);
    dist_bn = dist(b,n); 
    if (dist_an + dist_bn - dist_ab < std::numeric_limits<float>::epsilon()) {
        z[0] = n[0];
        z[1] = n[1];
        z[2] = n[2];
        dist_res = dist(x, n);
    } else {
        if (dist(x, a) < dist(x, b)) {
            dist_res = dist(x, a);
            z[0] = a[0];
            z[1] = a[1];
            z[2] = a[2];
        } else {
            dist_res = dist(x, b);
            z[0] = b[0];
            z[1] = b[1];
            z[2] = b[2];
        }
    }
}


void near_point(const double* a, const double* b, const double* x, double* z, double& dist_res) {
    double n[3];
    double dist_ab, dist_an, dist_bn;
    perp(a, b, x, n);
    dist_ab = dist(a,b);
    dist_an = dist(a,n);
    dist_bn = dist(b,n); 
    if (dist_an + dist_bn - dist_ab < std::numeric_limits<float>::epsilon()) {
        z[0] = n[0];
        z[1] = n[1];
        z[2] = n[2];
        dist_res = dist(x, n);
    } else {
        if (dist(x, a) < dist(x, b)) {
            dist_res = dist(x, a);
            z[0] = a[0];
            z[1] = a[1];
            z[2] = a[2];
        } else {
            dist_res = dist(x, b);
            z[0] = b[0];
            z[1] = b[1];
            z[2] = b[2];
        }
    }
}



//==============================Normal to cell edge(нормаль к краю)==============================
// rut0 - ячейка с 4 вершинами
// norm - нормаль к ячейке
// i, inext - индексы вершин A,B 
void get_nu(const double (&rut0)[4][3], const double* norm, const int i, const int inext, double* nu) {
    double len_nu, diff[3]; 
    for (int j = 0; j < 3; j++) {
        diff[j] = rut0[inext][j] - rut0[i][j];
    }
    vec_prod(diff, norm, nu);
    len_nu = vec_length(nu);

    if (len_nu > std::numeric_limits<float>::epsilon()) {
        for (int j = 0; j < 3; j++) {
            nu[j] /= len_nu;
        }
    }
}

//a, b - концы отрезка, norm - нормаль к ячейке
void get_nu(const double* a, const double* b, const double* norm, double* nu) {
    double len_nu, diff[3]; 
    for (int j = 0; j < 3; j++) {
        diff[j] = b[j] - a[j];
    }
    vec_prod(diff, norm, nu);
    len_nu = vec_length(nu);

    if (len_nu > std::numeric_limits<float>::epsilon()) {
        for (int j = 0; j < 3; j++) {
            nu[j] /= len_nu;
        }
    }
}



//========================================Cell diameter==========================================
double get_diam(const double (&root_tmp)[4][3]) {
    double diag1[3], diag2[3], res = 0.;
    for (int g = 0; g < 3; g++) {
        diag1[g] = root_tmp[2][g] - root_tmp[0][g];
        diag2[g] = root_tmp[3][g] - root_tmp[1][g];
    }
    if (vec_length(diag1) > vec_length(diag2)) {
        res = vec_length(diag1);
    } else {
        res = vec_length(diag2);
    }

    double ab[3], bc[3], cd[3], da[3];
    for (int g = 0; g < 3; g++) {
        ab[g] = root_tmp[1][g] - root_tmp[0][g];
        bc[g] = root_tmp[2][g] - root_tmp[1][g];
        cd[g] = root_tmp[3][g] - root_tmp[2][g];
        da[g] = root_tmp[0][g] - root_tmp[3][g];
    }

    if (vec_length(ab) > res) {
        res = vec_length(ab);
    }
    if (vec_length(bc) > res) {
        res = vec_length(bc);
    }
    if (vec_length(cd) > res) {
        res = vec_length(cd);
    }
    if (vec_length(da) > res) {
        res = vec_length(da);
    }
    return res;
}


double get_diam(const double (&root_tmp)[3][3]) {
    double len1, len2, len3;
    len1 = dist(root_tmp[0], root_tmp[1]);
    len2 = dist(root_tmp[1], root_tmp[2]);
    len3 = dist(root_tmp[2], root_tmp[0]);
    return std::max(std::max(len1, len2), len3);
}

//========================================Center of mass=========================================
void get_center_mass(const double* a, const double* b, const double* c, const double* d, double* res) {
    for (int i = 0; i < 3; i++) {
        res[i] = (a[i] + b[i] + c[i] + d[i]) * 0.25;
    }
}

void get_center_mass(const double* a, const double* b, const double* c, double* res) {
    for (int i = 0; i < 3; i++) {
        res[i] = (a[i] + b[i] + c[i]) / 3.;
    }
}

void get_center_mass(const double (&root_tmp)[4][3], double* res) {
    for (int i = 0; i < 3; i++) {
        res[i] = (root_tmp[0][i] + root_tmp[1][i] + root_tmp[2][i] + root_tmp[3][i]) * 0.25;
    }
}

void get_center_mass(const double (&root_tmp)[3][3], double* res) {
    for (int i = 0; i < 3; i++) {
        res[i] = (root_tmp[0][i] + root_tmp[1][i] + root_tmp[2][i]) / 3.;
    }
}


//=============================Check two points are different=======================================
int check_points_match(const double* a, const double* b) {
    int res = 0; //если точки различны, вернет 0
    for (int i = 0; i < 3; i++) {
        if (std::abs(a[i] - b[i]) < std::numeric_limits<float>::epsilon()) {
            res = 1;
            break;
        }
    }
    return res;
}



void vecs_orthoganalization(const double* vec1, const double* vec2,
                             double* vec1_orth, double* vec2_orth)
{
    // Normalization + ortogonalization
    for (int i = 0; i < 3; i++) {
        vec1_orth[i] = vec1[i];
        vec2_orth[i] = vec2[i];
    }

    double len1 = vec_length(vec1_orth);
    if (len1 > std::numeric_limits<float>::epsilon()) {
        vec1_orth[0] /= len1;
        vec1_orth[1] /= len1;
        vec1_orth[2] /= len1;
    }


    double sp = scal_prod(vec1_orth, vec2_orth);
    vec2_orth[0] -= sp * vec1_orth[0];
    vec2_orth[1] -= sp * vec1_orth[1];
    vec2_orth[2] -= sp * vec1_orth[2];

    double len2 = vec_length(vec2_orth);
    if (len2 > std::numeric_limits<float>::epsilon()) {
        vec2_orth[0] /= len2;
        vec2_orth[1] /= len2;
        vec2_orth[2] /= len2;
    }
}


} // namespace bielcc