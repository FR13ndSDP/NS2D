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

    double dt;

    Solver();

    ~Solver();

    //void getInput();
    
    void init(std::string name);

    void setBoundary_cylinder();
    
    void setBoundary_sod();

    void computeRHS();

    void computeDt(double t_now, double t_end);

    void showRHS();

    void computeFlux();

    void computeFlux_x();
    
    void computeFlux_y();

    void timeAdvance(double t_now, double t_end);

    void writeFile(std::string name);
};

#endif
