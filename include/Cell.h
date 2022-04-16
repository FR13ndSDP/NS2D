#ifndef _CELL_
#define _CELL_

#include <vector>
#include "Point.h"

class Cell {
public:
    // center point
    Point center;
    
    // neighbor points list
    /* 1 <- 0
       |    |
       2 -> 3
    */
    std::vector<Point> points;

    /* 
      x_lo, y_lo, x_hi, y_hi
      unit vector
    */
    std::vector<Point> faceNorm;

    // face area list
    std::vector<double> area;

    double volume;

    Cell() {
        volume = 0;
    };

    // construct with neighbor points
    Cell(Point p0, Point p1, Point p2, Point p3) {
        if (!points.empty()) {
            throw std::runtime_error("Cell is not empty");
        } else {
            points.push_back(p0);
            points.push_back(p1);
            points.push_back(p2);
            points.push_back(p3);
            getCenter();
            getFaceNorm();
            getVolume();
        }
    }
    ~Cell() {};

    void getCenter();

    void getFaceNorm();

    void getVolume();

    Point norm(int i, int j);

    friend std::ostream &operator<<(std::ostream &out, Cell &c) {
        out << "(" << c.center.coord_x << ", " << c.center.coord_y << ")\n"\
            << "volume = " << c.volume;
        return out;
    }
};

#endif