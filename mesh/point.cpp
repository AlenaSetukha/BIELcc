#include <iostream>
#include <stdexcept>

#include "point.h"
#include "tolerance_constants.h"

namespace bielcc {

bool Point::operator==(const Point& other) const {
    double dx = x - other.x;
    double dy = y - other.y;
    double dz = z - other.z;
    return (dx * dx + dy * dy + dz * dz) < Tolerance_Constants::POINT_TOLERANCE_SQ;
}


double Point::operator[](int i) const {
    switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2: 
            return z;
        default:
            throw std::runtime_error("Memory error: wrong index in Point");
    }
}





Point Point::operator+(const Point& other) const {
    return Point(x + other.x, y + other.y, z + other.z);
}

Point Point::operator-(const Point& other) const {
    return Point(x - other.x, y - other.y, z - other.z);
}

Point Point::operator*(double a) const {
    return Point(x * a, y * a, z * a);
}

// Сравнение в лексиграфическом порядке
bool Point::operator<(const Point& other) const
{
    if (x < other.x) {return true;}
    if (x > other.x) {return false;}
    if (y < other.y) {return true;}
    if (y > other.y) {return false;}
    if (z < other.z) {return true;}
    if (z > other.z) {return false;}
    return false;
}
} // namespace bielcc