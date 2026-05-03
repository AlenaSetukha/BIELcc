#ifndef _SEGMENT_H_
#define _SEGMENT_H_

#include <iostream>
#include <cmath>

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


    Segment(): nA(-1), nB(-1) {}
    ~Segment() = default;
    Segment(const Segment&) = default;
    Segment& operator=(const Segment&) = default;
    // Move-constructors
    Segment(Segment&&) noexcept = default;
    Segment& operator=(Segment&&) noexcept = default;


    Segment(int _na, int _nb): nA(_na), nB(_nb) {}




    bool operator==(const Segment& other) const;
    bool operator<(const Segment& other) const;
    double operator[](int i) const;



    void GetVertIndx(int (&vi)[2]) const {vi[0] = nA;
                                          vi[1] = nB;}
    int GetStartIndx() const {return nA;}
    int GetEndIndx() const {return nB;}
};

}       // namespace bielcc
#endif  // _SEGMENT_H