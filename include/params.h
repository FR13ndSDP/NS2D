#ifndef _PARAMS_
#define _PARAMS_
#include <cmath>

/*
y 80
|
|outlet
|
---------x
        240
*/
#define nx 240
#define ny 80
#define ncx (nx-1)
#define ncy (ny-1)
#define nc (ncx*ncy)
#define nGhost 2

#define R 287
#define G 1.4
#define Ma 6.0
#define Re 10000

#define rho_infty 1.29
#define dt 1e-6
#define Temprature 226.5
#define T_wall (4.2*Temprature)
#define soundSpeed sqrt(G * R * Temprature)

#endif