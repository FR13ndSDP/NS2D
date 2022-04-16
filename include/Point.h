#ifndef _POINT_
#define _POINT_
#include <iostream>
#include <cmath>

class Point {
public:
    double coord_x;
    double coord_y;
    Point()
    {
        coord_x = 0;
        coord_y = 0;
    }

    Point(double x, double y)
    {
        coord_x = x;
        coord_y = y;
    }

    ~Point() {};

    friend std::ostream &operator<<(std::ostream &out, Point &p) {
        out << p.coord_x << ", " << p.coord_y;
        return out;
    }

    // distance to (0,0)
    double D() {
        return sqrt(coord_x*coord_x + coord_y*coord_y);
    }

    // reverse direction
    Point operator-() {
        return Point(-1 * this->coord_x, -1 * this->coord_y);
    }

    Point operator-(Point &p) {
        return Point(this->coord_x - p.coord_x, this->coord_y - p.coord_y);
    }

    Point operator-=(Point &p) {
        coord_x -= p.coord_x;
        coord_y -= p.coord_y;
        return *this;
    }

    Point operator+(Point &p) {
        return Point(this->coord_x + p.coord_x, this->coord_y + p.coord_y);
    }

    Point operator+=(Point &p) {
        coord_x += p.coord_x;
        coord_y += p.coord_y;
        return *this;
    }

    // scalar division
    Point operator/(double d) {
        return Point(this->coord_y/d, this->coord_y/d);
    }

    Point operator/=(double d) {
        coord_x /= d;
        coord_y /= d;
        return *this;
    }

    // scalar product
    friend Point operator*(double d, Point &p) {
        return Point(p.coord_x * d, p.coord_y * d);
    }

    friend Point operator*(Point &p, double d) {
        return Point(p.coord_x * d, p.coord_y * d);
    }

    Point operator*=(double d) {
        coord_x *= d;
        coord_y *= d;
        return *this;
    }
    
    // vector dot product
    double operator*(Point &p) {
        return this->coord_x*p.coord_x + this->coord_y*p.coord_y;
    }

    double operator[](int index) {
        if (index == 0)
            return coord_x;
        else
            return coord_y;
    }
};

#endif