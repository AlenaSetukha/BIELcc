#ifndef _POINT_H_
#define _POINT_H_

#include <iostream>

namespace bielcc {

//===========================================================
//----------------------Point structure----------------------
//===========================================================
/**
 * Point structure consist of three coordinates {x, y, z}.
 * @note Two points considered identical with an accuracy POINT_TOLERANCE
*/

/**
    * @brief Point structure
*/

struct Point {
    double coords[3]{};
    
    Point() = default;
    ~Point() = default;
    Point(const Point&) = default;
    Point& operator=(const Point&) = default;
    // Move-constructors
    Point(Point&&) noexcept = default;
    Point& operator=(Point&&) noexcept = default;
    



    Point(double _x, double _y, double _z) {
        coords[0] = _x;
        coords[1] = _y;
        coords[2] =  _z;}
    Point(const double* pnt) {
        coords[0] = pnt[0];
        coords[1] = pnt[1];
        coords[2] = pnt[2];
    }

    bool operator==(const Point& other) const;
    bool operator<(const Point& other) const;
    double operator[](int i) const;
    Point operator+(const Point& other) const;
    Point operator-(const Point& other) const;
    Point operator*(double a) const;

    const double* GetCoord() const {return coords;}
};





/*
struct PointHash {
    std::size_t operator()(const Point& p) const {
        const double inv_tolerance = 1.0 / POINT_TOLERANCE;

        long long qx = static_cast<long long>(std::floor(p.x * inv_tolerance));
        long long qy = static_cast<long long>(std::floor(p.y * inv_tolerance));
        long long qz = static_cast<long long>(std::floor(p.z * inv_tolerance));

        std::size_t seed = 0;
        std::hash<long long> hasher;

        // Магическое число 0x9e3779b9 - это примерно золотое сечение,
        // используется для хорошего распределения хешей.
        seed ^= hasher(qx) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(qy) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(qz) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};*/

}       // namespace bielcc
#endif  // _POINT_H