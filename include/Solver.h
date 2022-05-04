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

    double **dt_local;

    Solver();

    ~Solver();

    //void getInput();
    
    void init(std::string name);

    void setBoundary();

    void computeRHS();

    void computeDt();

    void showRHS();

    void computeFlux();

    void computeFlux_x();
    
    void computeFlux_y();

    void timeAdvance();

    void writeFile(std::string name);

    double jacx(int i, int j);

    double jacy(int i, int j);

    double* vis_x(int i, int j);

    double* vis_y(int i, int j);

    void fill_corner();

    double mu_t(double T);

    double k_t(double T);
};

#endif
