#ifndef _SOLVER_
#define _SOLVER_
#include "Field.h"
#include "Mesh.h"
#include "Flux.h"
#include "params.h"
#include <string>

class Solver {
public:
    Field field;

    Mesh mesh;

    Flux f;
    
    double ***rhs;

    Solver();

    ~Solver();

    //void getInput();
    
    void init(std::string name);

    void setBoundary();

    void computeRHS();

    void computeFlux();

    void computeFlux_x();
    
    void computeFlux_y();

    void timeAdvance();

    void writeFile(std::string name);

    void writeDebug_F();
    
    void writeDebug_U();

    void writeDebug_G();

    double jacx(int i, int j);

    double jacy(int i, int j);

    double* vis_x(int i, int j);

    double* vis_y(int i, int j);
};

#endif