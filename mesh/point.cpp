#include <iostream>
#include <stdexcept>

#include "point.h"
#include "tolerance_constants.h"

namespace bielcc {

bool Point::operator==(const Point& other) const {
    const double dx = coords[0] - other.coords[0];
    const double dy = coords[1] - other.coords[1];
    const double dz = coords[2] - other.coords[2];
    return (dx * dx + dy * dy + dz * dz) < Tolerance_Constants::POINT_TOLERANCE_SQ;
}


double Point::operator[](int i) const {
    switch (i) {
        case 0:
            return coords[0];
        case 1:
            return coords[1];
        case 2: 
            return coords[2];
        default:
            throw std::runtime_error("Memory error: wrong index in Point");
    }
}





Point Point::operator+(const Point& other) const {
    return Point(coords[0] + other.coords[0],
                 coords[1] + other.coords[1],
                 coords[2] + other.coords[2]);
}

Point Point::operator-(const Point& other) const {
    return Point(coords[0] - other.coords[0],
                 coords[1] - other.coords[1],
                 coords[2] - other.coords[2]);
}

Point Point::operator*(double a) const {
    return Point(coords[0] * a, coords[1] * a, coords[2] * a);
}

// Сравнение в лексиграфическом порядке
bool Point::operator<(const Point& other) const
{
    if (std::abs(coords[0] - other.coords[0]) > Tolerance_Constants::POINT_TOLERANCE) {
        return coords[0] < other.coords[0];
    }
    if (std::abs(coords[1] - other.coords[1]) > Tolerance_Constants::POINT_TOLERANCE) {
        return coords[1] < other.coords[1];
    }
    if (std::abs(coords[2] - other.coords[2]) > Tolerance_Constants::POINT_TOLERANCE) {
        return coords[2] < other.coords[2];
    }
    return false; // Точки равны с учетом Tolerance_Constants::POINT_TOLERANCE
}
} // namespace bielcc