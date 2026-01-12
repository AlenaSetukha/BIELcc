#include <array>
#include <cmath>
#include <complex>
#include <iostream>

#include "integrals_analytic.h"
#include "element_geom.h"

namespace bielcc {
//===========================================================
//-------------------Integral 1/|x-y|------------------------
//------------A triangular cell is assumed below-------------
//===========================================================
// compute rotation
static void lartgp(double x, double y, double &cs, double &sn, double &r)
{
    double tmp = sqrt(x * x + y * y);
    if (tmp == 0) {
        r = 0, cs = 1, sn = 0;
    }
    else {
        cs = x / tmp;
        sn = y / tmp;
        r = tmp;
    }
}

// apply rotation
static void rot(double &x, double &y, double cs, double sn)
{
    double tmp = x * cs + y * sn;
    y = y * cs - x * sn;
    x = tmp;
}

// inner function for integration 1 / |x-y|, x in triangle ABC.
static double
_Integral1Divr_inv(std::array<std::array<double, 3>, 3> const &ABC,
                   std::array<double, 3> const &D)
{
    double p0, l, l_pl, ratio, r, ad, proj;
    double n[3], edges[3][3], Q[3][2], DABC[7][3], c[3], s[3];
    int ed, emn;

    // calculate vectors collinear with to edges of the triangle.

    for (int i = 0; i < 3; ++i) {
        edges[i][0] = ABC[i][1] - ABC[i][0]; // AB
        edges[i][1] = ABC[i][2] - ABC[i][1]; // BC
        edges[i][2] = ABC[i][0] - ABC[i][2]; // CA
    }

    // Compute edge length to select largest one and normalize
    // Also compute direction from D to vertices of triangle
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            DABC[j][i] = ABC[j][i] - D[j]; // DA, DB, DC
        }
        // |DA|, |DB|, |DC|
        DABC[3][i] = sqrt(DABC[0][i] * DABC[0][i] + DABC[1][i] * DABC[1][i] +
                          DABC[2][i] * DABC[2][i]);
        Q[i][0] = edges[i][0];
        Q[i][1] = edges[i][1];
        n[i] = 0;
    }

    DABC[4][0] = DABC[3][1]; // |DB|
    DABC[4][1] = DABC[3][2]; // |DC|
    DABC[4][2] = DABC[3][0]; // |DA|

    // Compute QR with pivoting for triangle edges AB, BC
    lartgp(Q[0][0], Q[1][0], c[0], s[0], Q[0][0]);

    lartgp(Q[0][0], Q[2][0], c[1], s[1], Q[0][0]);
    rot(Q[0][1], Q[1][1], c[0], s[0]);
    rot(Q[0][1], Q[2][1], c[1], s[1]);
    lartgp(Q[1][1], Q[2][1], c[2], s[2], Q[1][1]);

    // S_ABC = Q(1, 1) * Q(2, 2) / 2

    // As edges are linearly dependent Q e_3 is orthogonal to the plane containing
    // triangle This vector is unit and can be used as normal
    n[2] = 1.0;

    // Multiplication of Q by e_3 itself
    rot(n[1], n[2], c[2], -s[2]);
    rot(n[0], n[2], c[1], -s[1]);
    rot(n[0], n[1], c[0], -s[0]);

    // Now we will construct QR decomposition using rotations
    // (in fact, reflections can be used as well, but rotations guarantees correct
    // orientation of columns in matrix Q) for matrices [n edges(:, 1)] = Q_1 R_1,
    // [n edges(:, 2)] = Q_2 R_2 and [n edges(:, 3)] = Q_3 R_3 It will give us
    // orhonormal basises Q_1, Q_2, Q_3 as it needs Q_1^T DABC(:, 1), Q_2^T
    // DABC(:, 2), Q_3^T DABD(:, 3) and [Q_1^T DABC(:, 2)]_2, [Q_2^T DABC(:, 3)]_2
    // and [Q_3^T DABC(:, 1)]_2 It is easy to note, that each decomposition
    // requires 3 rotations and first two of them is the same for all the three
    // matrices

    // Compute rotation to zero 2nd element of normal
    lartgp(n[0], n[1], c[0], s[0], n[0]); // First rotation to save

    for (int i = 0; i < 3; ++i)
    {
        // Apply this rotation to edges
        // call drot(3, edges(1, 1), 3, edges(2, 1), 3, c(1), s(1))
        rot(edges[0][i], edges[1][i], c[0], s[0]); // AB, BC, CA

        // Also we can apply it to target vectors
        // call drot(3, DABC(1, 1), 7, DABC(2, 1), 7, c(1), s(1))
        rot(DABC[0][i], DABC[1][i], c[0], s[0]); // XA, XB, XC
    }
    // The same here to zero 3rd element of normal
    lartgp(n[0], n[2], c[0], s[0], n[0]); // Second rotation to save
    for (int i = 0; i < 3; ++i)
    {
        rot(edges[0][i], edges[2][i], c[0], s[0]);
        rot(DABC[0][i], DABC[2][i], c[0], s[0]);
    }

    // Now compute last rotation for each edges
    lartgp(edges[1][0], edges[2][0], c[0], s[0],
           edges[1][0]); // Third rotation for AB
    lartgp(edges[1][1], edges[2][1], c[1], s[1],
           edges[1][1]); // Third rotation for BC
    lartgp(edges[1][2], edges[2][2], c[2], s[2],
           edges[1][2]); // Third rotation for CA
    // Compute [Q_1^T DABC(:, 2)]_2, [Q_2^T DABC(:, 3)]_2 and [Q_3^T DABC(:, 1)]_2
    DABC[5][0] = c[0] * DABC[1][1] + s[0] * DABC[2][1];  // (XB, AB)
    DABC[5][1] = c[1] * DABC[1][2] + s[1] * DABC[2][2];  // (XC, BC)
    DABC[5][2] = c[2] * DABC[1][0] + s[2] * DABC[2][0];  // (XA, CA)
    DABC[6][0] = -s[0] * DABC[1][1] + c[0] * DABC[2][1]; // (XB, AB_ort)
    DABC[6][1] = -s[1] * DABC[1][2] + c[1] * DABC[2][2]; // (XC, BC_ort)
    DABC[6][2] = -s[2] * DABC[1][0] + c[2] * DABC[2][0]; // (XA, CA_ort)
    // And now apply last rotations
    rot(DABC[1][0], DABC[2][0], c[0], s[0]);
    rot(DABC[1][1], DABC[2][1], c[1], s[1]);
    rot(DABC[1][2], DABC[2][2], c[2], s[2]);

    // Now DABC contains all the values we want
    //
    // loop by edges of triangle ABC:
    //
    double int_1r = 0;
    proj = 0;
    for (ed = 0; ed < 3; ++ed) {
        ad = 0;
        /*
         * calculate the part of integral for function 1/r that connected
         * with edge(ver)
         */

        r = DABC[3][ed];

        if (std::abs(r) > 0 && std::abs(DABC[4][ed]) > 0) {
            // if p0 is very small all parts are small, so putting them zero
            //  Projection to normal to edge in plane.
            //  Minus because of left orientation (originally normal was in the middle
            //  and now it is swapped with the edge)
            p0 = -DABC[2][ed] / r;

            // Projection to edge
            l = DABC[1][ed] / r;

            // ratio of r_pl and r
            ratio = DABC[4][ed] / r;

            // Projection of next vector to the edge
            l_pl = DABC[5][ed] / DABC[4][ed];
            // Note, that r_pl = 0 means that point D is the same as the next vertice,
            // so direction from the current vertice to D is collinear to the edge.
            // Thus othogonal part is 0, so p0 = 0. r = 0 also means p0 = 0. l = -1 or
            // l_pl = -1 means that vector from point D to the current vertice or to
            // the next vertice is is collinear to the edge connecting those two
            // vertices, so orthogonal part of vector from D to current point is zero,
            // so p_0 = 0 So, problems here are possible only if p0 = 0, and we have
            // already checked this

            if (l > -1 && l_pl > -1) {
                ad = p0 * std::log(ratio * (1 + l_pl) / (1 + l));
            }

            // proj = (XA, XB) * |XC| + (XB, XC) * |XA| + (XC, XA) * |XB|
            emn = (ed + 2) % 3; // emn == ed -1 mod 3
            proj = proj + (DABC[0][ed] * DABC[0][ed] + DABC[1][ed] * DABC[5][ed] +
                           DABC[2][ed] * DABC[6][ed]) *
                              DABC[3][emn];

            int_1r = int_1r + r * ad;
        }
    }

    c[0] = std::abs(DABC[0][0]) * Q[0][0] * Q[1][1]; // (XA, XB, XC)
    c[1] =
        DABC[3][0] * DABC[3][1] * DABC[3][2] + proj; // |XA| * |XB| * |XC| + ...
    if (c[0] > 0 || std::abs(c[1]) > 0) {
        int_1r = int_1r - std::abs(DABC[0][0]) * 2 * std::atan2(c[0], c[1]);
    }
    return int_1r;
}




// integrate 1/|x-y|, y in Triangle t.
static double integral_1Divr(const double (*t)[3], const double *x) {
    std::array<std::array<double, 3>, 3> abc;
    std::array<double, 3> d({x[0], x[1], x[2]});
    abc[0][0] = t[0][0], abc[0][1] = t[1][0], abc[0][2] = t[2][0];
    abc[1][0] = t[0][1], abc[1][1] = t[1][1], abc[1][2] = t[2][1];
    abc[2][0] = t[0][2], abc[2][1] = t[1][2], abc[2][2] = t[2][2];
    return _Integral1Divr_inv(abc, d);
}




double integral1Divr(const double (&rut0)[4][3], const double *x) {
    double tr1[3][3], tr2[3][3];
    for (int i = 0; i < 3; i++)
    {
        tr1[0][i] = rut0[0][i];
        tr1[1][i] = rut0[1][i];
        tr1[2][i] = rut0[2][i];

        tr2[0][i] = rut0[2][i];
        tr2[1][i] = rut0[3][i];
        tr2[2][i] = rut0[0][i];
    }
    return (integral_1Divr(tr1, x) + integral_1Divr(tr2, x));
}

double integral1Divr(const double (&rut0)[3][3], const double *x) {
    return integral_1Divr(rut0, x);
}




void integralGrad1DivrSeg(const double* A, const double* B,
                            const double *x, double* res)
{
    double xmA[3] = {x[0] - A[0],
                     x[1] - A[1],
                     x[2] - A[2]};
    double xmB[3] = {x[0] - B[0],
                     x[1] - B[1],
                     x[2] - B[2]};
    double xA = dist(x, A);
    double xB = dist(x, B);
    double AB = dist(A, B);

    for (int i = 0; i < 3; i++) {
        res[i] = - AB * (xmA[i] / xA + xmB[i] / xB) / 
                  (scal_prod(xmA, xmB) + xA * xB);
    }
}











//===========================================================
//---------------Integral d / dn(1 / |x - y|)----------------
//===========================================================
double integral_ddn_1Divr(const double (&rut0)[4][3], const double *x)
{
    double res = solid_angle(rut0[0], rut0[1], rut0[2], x) +
          solid_angle(rut0[2], rut0[3], rut0[0], x);
    return res;
}







//===========================================================
//-----------Integral (x - y) / (|x - y| ^ 2)----------------
//===========================================================
void integral_xmyDivr2(const double (&rut0)[4][3], const double *x,
                                                        double *res)
{
    double XH, PX[3], PQ[3], QX[3], QP[3], norm[3];
    double e1[3], e2[3];
    double alpha, beta, sina, sinb, len, i1, i2;
    int next;
    for (int i = 0; i < 3; i++) {
        res[i] = 0.;
    }

    for (int i = 0; i < 4; i++) {
        next = (i + 1) % 4; // div_mod(i + 1, 4);
        for (int j = 0; j < 3; j++) {
            PX[j] = x[j] - rut0[next][j];       // X - P
            PQ[j] = rut0[i][j] - rut0[next][j]; // Q - P
            QX[j] = x[j] - rut0[i][j];          // X - Q
            QP[j] = rut0[next][j] - rut0[i][j]; // P - Q
        }
        // углы
        sina = -scal_prod(PX, PQ) / vec_length(PX) / vec_length(PQ);
        sinb = scal_prod(QX, QP) / vec_length(QX) / vec_length(QP);
        alpha = asin(sina);
        beta = asin(sinb);
        // нормаль
        vec_prod(PX, PQ, norm);
        len = vec_length(norm);
        XH = len / vec_length(PQ);
        for (int j = 0; j < 3; j++) {
            norm[j] /= len;
        }
        // ОНБ базис
        vec_prod(norm, PQ, e1);
        len = vec_length(e1);
        for (int j = 0; j < 3; j++) {
            e1[j] /= len;
        }

        vec_prod(e1, norm, e2);
        len = vec_length(e2);
        for (int j = 0; j < 3; j++) {
            e2[j] /= len;
        }
        // I1, I2
        i1 = alpha - beta;
        i2 = log(fabs(cos(beta))) - log(fabs(cos(alpha)));
        for (int j = 0; j < 3; j++) {
            res[j] += XH * (e1[j] * i1 + e2[j] * i2);
        }
    }
}








//===========================================================
//----------Integral nu / |x - y| (curvilinear)--------------
//===========================================================
void integralnu1Divr(const double (&rut0)[4][3], const double *x,
                                    const double eps, double *res)
{
    int next;
    double BX[3], AB[3], BA[3], XA[3], XB[3], AX[3], num, denom, nu[3];
    double norm[3];
    norm_func(rut0, norm);

    for (int i = 0; i < 3; i++) {
        res[i] = 0.;
    }

    for (int i = 0; i < 4; i++) {
        next = (i + 1) % 4; // div_mod(i + 1, 4);
        get_nu(rut0, norm, i, next, nu);
        for (int j = 0; j < 3; j++) {
            BX[j] = x[j] - rut0[next][j];       // x - B
            XB[j] = rut0[next][j] - x[j];       // B - x
            AB[j] = rut0[next][j] - rut0[i][j]; // B - A
            BA[j] = rut0[i][j] - rut0[next][j]; // A - B
            XA[j] = rut0[i][j] - x[j];          // A - x
            AX[j] = x[j] - rut0[i][j];          // x - A
        }

        if (scal_prod(BX, AB) >= 0) { // формула 24, стр 998, Сетуха/Семенова// проверить условие
            num = scal_prod(BA, XA) + vec_length(AX) * vec_length(BA);
            denom = scal_prod(BA, XB) + vec_length(BX) * vec_length(BA);
        } else { // формула 23, стр 997
            num = scal_prod(AB, XB) + vec_length(BX) * vec_length(AB);
            denom = scal_prod(AB, XA) + vec_length(AX) * vec_length(AB);
        }

        if (denom < eps) {
            for (int j = 0; j < 3; j++) {
                res[j] += nu[j] * log(fabs(num / eps));
            }
        } else {
            if (fabs(num / denom) < eps) {
                for (int j = 0; j < 3; j++) {
                    res[j] += nu[j] * log(eps);
                }
            } else {
                for (int j = 0; j < 3; j++) {
                    res[j] += nu[j] * log(fabs(num / denom));
                }
            }
        }
    }
}


} // namespace bielcc