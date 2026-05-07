#include <iostream>

#include "cell.h"
#include "element_geom.h"


namespace bielcc {
void CellQuad::CellFill(const double (&cell_in)[4][3]) {
    for (int i = 0; i < 4; i++) {
        _cell[i][0] = cell_in[i][0];
        _cell[i][1] = cell_in[i][1];
        _cell[i][2] = cell_in[i][2];
    }

    for (int i = 0; i < 3; i++) {
        _rkt[i] = (_cell[0][i] + _cell[1][i] + _cell[2][i] + _cell[3][i]) * 0.25;
    }



    double a[3], b[3], c[3], d[3], ac[3], bd[3];
    for (int i = 0; i < 3; i++) {
        a[i] = (_cell[0][i] + _cell[1][i]) / 2.0;
        b[i] = (_cell[1][i] + _cell[2][i]) / 2.0;
        c[i] = (_cell[2][i] + _cell[3][i]) / 2.0;
        d[i] = (_cell[3][i] + _cell[0][i]) / 2.0;
        ac[i] = c[i] - a[i];
        bd[i] = d[i] - b[i];
    }
    vec_prod(ac, bd, _norm);
    double len = vec_length(_norm);
    for (int k = 0; k < 3; k++) {
        _norm[k] /= len;
    }

    
    _s = len;
    

    // tau[0], tau[1]
    len = vec_length(ac);
    for (int k = 0; k < 3; k++) {
        _tau[0][k] = ac[k] / len;
    }

    vec_prod(_norm, _tau[0], _tau[1]);
    len = vec_length(_tau[1]);
    for (int k = 0; k < 3; k++) {
        _tau[1][k] /= len;
    }
}









void CellTriangle::CellFill(const double (&cell_in)[3][3]) {
    for (int i = 0; i < 3; i++) {
        _cell[i][0] = cell_in[i][0];
        _cell[i][1] = cell_in[i][1];
        _cell[i][2] = cell_in[i][2];
    }

    for (int i = 0; i < 3; i++) {
        _rkt[i] = (_cell[0][i] + _cell[1][i] + _cell[2][i]) / 3.;
    }
    
    double ac[3], ab[3], len;

    for (int k = 0; k < 3; k++) {
        ab[k] = _cell[1][k] - _cell[0][k];
        ac[k] = _cell[2][k] - _cell[0][k];
    }

    vec_prod(ab, ac, _norm);
    len = vec_length(_norm);
    
    for (int k = 0; k < 3; k++) {
        _norm[k] /= len;
    }

    
    _s = len / 2;

    // tau[0] и tau[1]
    len = vec_length(ab);
    for (int k = 0; k < 3; k++) {
        _tau[0][k] = ab[k] / len;
    }

    vec_prod(_norm, _tau[0], _tau[1]);
    len = vec_length(_tau[1]);
    for (int k = 0; k < 3; k++) {
        _tau[1][k] /= len;
    }
}


} //namespace bielcc