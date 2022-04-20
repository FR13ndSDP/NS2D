#ifndef _FIELD_
#define _FIELD_
#include "params.h"
#include "Mesh.h"

class Field {
public:
    // conservative variables with ghost cell
    double ***U;

    // old time step value
    double ***Un;

    // primitive with ghost cell
    // rho, u, v, p
    double ***prim;

    int ngx = ncx+ 2*nGhost;
    int ngy = ncy + 2*nGhost;

    // inlet for cons
    double inlet[4] = {
        rho_infty,
        rho_infty * Ma * soundSpeed,
        0.0,
        rho_infty *Temprature* R/(G-1) + \
        0.5 * rho_infty * (Ma * soundSpeed) * (Ma * soundSpeed)
    };

    Field();

    ~Field();

    void cons2prim();

    void prim2cons();

    void init();

    void initFromFile(std::string name);

    void U2Un();

    double Temp(int i, int j);
};

#endif