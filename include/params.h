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
#define Ma 8.0
#define Pr 0.77

#define rho_infty 8.148e-5
#define dt 1e-7
#define Temprature 226.5
#define T_wall (4.2*Temprature)
#define soundSpeed sqrt(G * R * Temprature)

#define plot_interval 1000
#define show_interval 10

#endif