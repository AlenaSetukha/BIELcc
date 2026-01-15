#include <iostream>
#include <stdexcept>

#include "segment.h"

namespace bielcc {

Segment::Segment(const Segment& obj) {
    nA = obj.nA;
    nB = obj.nB;
    for (int i = 0; i < 3; i++) {
        ACoords[i] = obj.ACoords[i];
        BCoords[i] = obj.BCoords[i];
    }
}

Segment::Segment(int _na, int _nb, const double* _ACoords, const double* _BCoords) {
    nA = _na;
    nB = _nb;
    for (int i = 0; i < 3; i++) {
        ACoords[i] = ACoords[i];
        BCoords[i] = BCoords[i];
    }
}

bool Segment::operator==(const Segment& other) const {
    return ((nA == other.nA && nB == other.nB) || 
            (nA == other.nB && nB == other.nA));
}

bool Segment::operator<(const Segment& other) const {
    if (nA < other.nA) {return true;}
    if (nA > other.nA) {return false;}
    if (nB < other.nB) {return true;}
    if (nB > other.nB) {return false;}
    return true;  
}


double Segment::operator[](int i) const {
    switch (i) {
        case 0:
            return nA;
        case 1:
            return nB;
        default:
            throw std::runtime_error("Memory error: wrong index in Segment");
    }
}


void Segment::SetCoords(const double (&_ACoords)[3], const double (&_BCoords)[3]) {
    for (int i = 0; i < 3; i++) {
        ACoords[i] = _ACoords[i];
        BCoords[i] = _BCoords[i];
    }
}


void Segment::GetVertCoord(double (&_ACoords)[3], double (&_BCoords)[3]) const {
    for (int i = 0; i < 3; i++) {
        _ACoords[i] = ACoords[i];
        _BCoords[i] = BCoords[i];
    }
}

} // namespace bielcc