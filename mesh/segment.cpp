#include <iostream>
#include <stdexcept>

#include "segment.h"

namespace bielcc {
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

} // namespace bielcc