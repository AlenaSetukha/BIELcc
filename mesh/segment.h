#ifndef _SEGMENT_H_
#define _SEGMENT_H_

#include <iostream>

namespace bielcc {

//===========================================================
//--------------------Segment structure----------------------
//===========================================================
/**
 * Segment structure consist of two points.
 * Fields:
 *      nA - starting point number from the list
 *      nB - ending point number from the list
 *      ACoords - starting point coordinates
 *      BCoords - ending point coordinates
 * @note Segment is ordered by point numbers (nA <= nB), comparison on the same principle.
 */

/**
 * @brief Segment structure
 */
struct Segment {
    int nA, nB;
    double ACoords[3]{}, BCoords[3]{};


    Segment(): nA(0), nB(0) {}
    Segment(int _na, int _nb) : nA(_na), nB(_nb) {}
    Segment(int _na, int _nb, const double* _ACoords, const double* _BCoords);
    Segment(const Segment& obj);
    ~Segment() = default;


    bool operator==(const Segment& other) const;
    bool operator<(const Segment& other) const;
    double operator[](int i) const;


    void SetCoords(const double (&_ACoords)[3], const double (&_BCoords)[3]);

    void GetVertIndx(int (&vi)[2]) const {vi[0] = nA;
                                          vi[1] = nB;}
    void GetVertCoord(double (&_ACoords)[3], double (&_BCoords)[3]) const;
};

}       // namespace bielcc
#endif  // _SEGMENT_H