#ifndef _CELL_H_
#define _CELL_H_

#include <iostream>

namespace bielcc {

//===========================================================
//-----------------Quadrangular cell class-------------------
//===========================================================
/**
 * Quadrangular cell class.
 * Fields:
 *      _cell - coordinates of 4 vertices
 *      _tau - coordinates of two local basis vectors
 *      _rkt - collocation point coordinates 
 *      _norm - normal vector coordinates 
 *      _s - cell area
 * 
 * @note Cell can be pseudo-quadrangular, two vertices can be the same
*/

/**
 * @brief Quadrangular cell class.
 */

class CellQuad {
protected:
    double _cell[4][3];
    double _tau[2][3];
    double _rkt[3];
    double _norm[3];
    double _s;

public:
    CellQuad() = default;
    CellQuad(const CellQuad& obj);
    ~CellQuad() = default;

    CellQuad& operator=(const CellQuad& obj);

    void CellFill(const double (&cell_in)[4][3]);


    void GetCellCoords(double (&cell)[4][3]) const {
        for (int i = 0; i < 4; i++) {
            cell[i][0] = _cell[i][0];
            cell[i][1] = _cell[i][1];
            cell[i][2] = _cell[i][2];
        }
    }
    const double (&GetCellCoords() const)[4][3] {return _cell;}

    const double (&GetTau() const)[2][3] {return _tau;}

    const double* GetTau1() const {return _tau[0];}
    double* GetTau1() {return _tau[0];}

    const double* GetTau2() const {return _tau[1];}
    double* GetTau2() {return _tau[1];}

    const double* GetRkt() const {return _rkt;}
    double* GetRkt() {return _rkt;}

    const double* GetNorm() const {return _norm;}
    double* GetNorm() {return _norm;}

    const double* GetVertex(int i) const {return _cell[i];} 
    double* GetVertex(int i) {return _cell[i];}

    double GetArea() {return _s;};
};








//===========================================================
//------------------Triangular cell class--------------------
//===========================================================
/**
 * Triangular cell class.
 * Fields:
 *      _cell - coordinates of 3 vertices
 *      _tau - coordinates of two local basis vectors
 *      _rkt - collocation point coordinates 
 *      _norm - normal vector coordinates 
 *      _s - cell area
 */

/**
 * @brief Triangular cell class.
 */
class CellTriangle {
private:
    double _cell[3][3];
    double _tau[2][3];
    double _rkt[3];
    double _norm[3];
    double _s;

public:
    CellTriangle() = default;
    CellTriangle(const CellTriangle& obj);
    ~CellTriangle() = default;

    CellTriangle& operator=(const CellTriangle& obj);

    void CellFill(const double (&cell_in)[3][3]);

    void GetCellCoords(double (&cell)[3][3]) const {
        for (int i = 0; i < 3; i++) {
            cell[i][0] = _cell[i][0];
            cell[i][1] = _cell[i][1];
            cell[i][2] = _cell[i][2];
        }
    }

    const double (&GetCellCoords() const)[3][3] {return _cell;}

    const double (&GetTau() const)[2][3] {return _tau;}

    const double* GetTau1() const {return _tau[0];}
    double* GetTau1() {return _tau[0];}

    const double* GetTau2() const {return _tau[1];}
    double* GetTau2() {return _tau[1];}

    const double* GetRkt() const {return _rkt;}
    double* GetRkt() {return _rkt;}

    const double* GetNorm() const {return _norm;}
    double* GetNorm() {return _norm;}

    const double* GetVertex(int i) const {return _cell[i];} // should not be delated by the caller!!!
    double* GetVertex(int i) {return _cell[i];}
    
    double GetArea() {return _s;};
};


} // namespace bielcc

#endif // _CELL_H_