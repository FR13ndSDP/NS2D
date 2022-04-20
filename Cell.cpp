#include "./include/Cell.h"

void Cell::getVolume() {
    Point vec1 = points[2] - points[0];
    Point vec2 = points[3] - points[1];
    // cross product
    volume = (vec2[0]*vec1[1] - vec2[1]*vec1[0]) * 0.5;
}

void Cell::getCenter() {
    double x = 0, y = 0;
    for (int i=0; i<4; i++) {
        x += points[i][0];
        y += points[i][1];
    }
    center.coord_x = x/4.0;
    center.coord_y = y/4.0;
}

Point Cell::norm(int i, int j) {
    Point vec(points[i][1]-points[j][1], points[j][0] - points[i][0]);
    double A = vec.D();
    vec /= A;
    area.push_back(A);
    return vec;
}

void Cell::getFaceNorm() {
    faceNorm.push_back(norm(0,1));
    faceNorm.push_back(norm(1,2));
    faceNorm.push_back(norm(2,3));
    faceNorm.push_back(norm(3,0));
}
