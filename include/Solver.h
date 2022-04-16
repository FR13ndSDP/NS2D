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
    
    void init();

    void setBoundary();

    void computeRHS();

    void computeFlux();

    void timeAdvance();

    void writeFile(std::string name);

    void writeDebug();
};

#endif