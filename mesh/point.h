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
    double x, y, z;
    Point(): x(0.0), y(0.0), z(0.0) {}
    Point(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    Point(const double* pnt): x(pnt[0]), y(pnt[1]), z(pnt[2]) {}
    Point(const Point& obj) {x = obj.x;
                             y = obj.y;
                             z = obj.z; };

    ~Point() = default;


    bool operator==(const Point& other) const;
    bool operator<(const Point& other) const;
    double operator[](int i) const;
    Point operator+(const Point& other) const;
    Point operator-(const Point& other) const;
    Point operator*(double a) const;


    void GetCoord(double (&coord)[3]) const {coord[0] = x;
                                             coord[1] = y;
                                             coord[2] = z;}
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