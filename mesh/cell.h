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
    * Cell can be pseudo-quadrangular, two vertices can be the same.
*/

/**
    * @brief Quadrangular cell class.
    * @note Cell can be pseudo-quadrangular, two vertices can be the same.
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
    ~CellQuad() = default;
    CellQuad(const CellQuad&) = default;
    CellQuad& operator=(const CellQuad&) = default;
    // Move-constructors
    CellQuad(CellQuad&&) noexcept = default;
    CellQuad& operator=(CellQuad&&) noexcept = default;


    /**
        * @brief Cell creation by the vertices coordinates
        * @param cell_in Vertices coordinates
    */
    void CellFill(const double (&cell_in)[4][3]);

    /**
        * @brief Coordinates of cell vertices
        * @param cell Copy of vertices coordinate array
    */
    void GetCellCoords(double (&cell)[4][3]) const {
        for (int i = 0; i < 4; i++) {
            cell[i][0] = _cell[i][0];
            cell[i][1] = _cell[i][1];
            cell[i][2] = _cell[i][2];
        }
    }
    /**
        * @brief Coordinates of cell vertices
        * @return Const reference to the vertices coordinate array
    */
    const double (&GetCellCoords() const)[4][3] {return _cell;}

    /**
        * @brief Coordinates of two local basis vectors
        * @return Const reference to the local basis vectors coordinate array
    */
    const double (&GetTau() const)[2][3] {return _tau;}

    /**
        * @brief Coordinates of first local basis vector
        * @return Const reference to the first local basis vector coordinate array
    */
    const double* GetTau1() const {return _tau[0];}

    /**
        * @brief Coordinates of second local basis vector
        * @return Const reference to the second local basis vector coordinate array
    */
    const double* GetTau2() const {return _tau[1];}

    /**
        * @brief Coordinates of collocation point
        * @return Const reference to the collocation point coordinates
    */
    const double* GetRkt() const {return _rkt;}

    /**
        * @brief Coordinates of cell normal vector
        * @return Const reference to the cell normal vector coordinates
    */
    const double* GetNorm() const {return _norm;}

    /**
        * @brief Coordinates of i-th [0-3] vertex coordinates
        * @param i Vertex number
        * @return Const reference to the i-th [0-3] vertex coordinates
    */
    const double* GetVertex(int i) const {return _cell[i];} 

    /**
        * @brief Cell area
        * @return Cell area
    */
    double GetArea() const {return _s;};
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
    * @brief Triangular cell class
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
    ~CellTriangle() = default;
    CellTriangle(const CellTriangle&) = default;
    CellTriangle& operator=(const CellTriangle&) = default;
    CellTriangle(CellTriangle&&) noexcept = default;
    CellTriangle& operator=(CellTriangle&&) noexcept = default;



    /**
        * @brief Cell creation by the vertices coordinates
        * @param cell_in Vertices coordinates
    */
    void CellFill(const double (&cell_in)[3][3]);

    /**
        * @brief Coordinates of cell vertices
        * @param cell Copy of vertices coordinate array
    */
    void GetCellCoords(double (&cell)[3][3]) const {
        for (int i = 0; i < 3; i++) {
            cell[i][0] = _cell[i][0];
            cell[i][1] = _cell[i][1];
            cell[i][2] = _cell[i][2];
        }
    }

    /**
        * @brief Coordinates of cell vertices
        * @return Const reference to the vertices coordinate array
    */
    const double (&GetCellCoords() const)[3][3] {return _cell;}

    /**
        * @brief Coordinates of two local basis vectors
        * @return Const reference to the local basis vectors coordinate array
    */
    const double (&GetTau() const)[2][3] {return _tau;}

    /**
        * @brief Coordinates of first local basis vector
        * @return Const reference to the first local basis vector coordinate array
    */
    const double* GetTau1() const {return _tau[0];}

    /**
        * @brief Coordinates of second local basis vector
        * @return Const reference to the second local basis vector coordinate array
    */
    const double* GetTau2() const {return _tau[1];}

    /**
        * @brief Coordinates of collocation point
        * @return Const reference to the collocation point coordinates
    */
    const double* GetRkt() const {return _rkt;}

    /**
        * @brief Coordinates of cell normal vector
        * @return Const reference to the cell normal vector coordinates
    */
    const double* GetNorm() const {return _norm;}

    /**
        * @brief Coordinates of i-th [0-3] vertex coordinates
        * @param i Vertex number
        * @return Const reference to the i-th [0-3] vertex coordinates
    */
    const double* GetVertex(int i) const {return _cell[i];} // should not be delated by the caller!!!
    
    /**
        * @brief Cell area
        * @return Cell area
    */
    double GetArea() const {return _s;};
};


} // namespace bielcc
#endif // _CELL_H_