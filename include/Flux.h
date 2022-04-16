#ifndef _FLUX_
#define _FLUX_
#include "params.h"

class Flux {
public:
    int nfx = nx;
    int nfy = ny;

    // face centered flux
    double ***Fx;
    double ***Fy;

    Flux();

    ~Flux();
    
    void fluxSplit(double *primL, double *primR, double *Flux);
};

#endif