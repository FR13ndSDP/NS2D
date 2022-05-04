#include "./include/Solver.h"

double Solver::mu_t(double T) {
    double mu_0 = 1.716e-5, T_ref = 273.15, T_0 = 111;
    double mu = mu_0*pow(T/T_ref, 1.5)*(T_ref + T_0)/(T + T_0);
    return mu;
}

double Solver::k_t(double T) {
    double cp = G * R/(G-1);
    double k = mu_t(T) * cp / Pr;
    return k;
}

// in I-1/2, x direction
double Solver::jacx(int i, int j)
{
    double Dix = mesh.cell(i,j).center[0] - mesh.cell(i-1,j).center[0];
    double Diy = mesh.cell(i,j).center[1] - mesh.cell(i-1,j).center[1];

    double Djx = 0.25 *(
                 mesh.cell(i-1,j+1).center[0] \
                 + mesh.cell(i,j+1).center[0]\
                 - mesh.cell(i-1, j-1).center[0]\
                 - mesh.cell(i, j-1).center[0]
                 );
    double Djy = 0.25 *(
                 mesh.cell(i-1,j+1).center[1] \
                 + mesh.cell(i,j+1).center[1]\
                 - mesh.cell(i-1, j-1).center[1]\
                 - mesh.cell(i, j-1).center[1]
                 );
    return 1.0/(Dix*Djy - Djx*Diy);
}

// in I-1/2, y direction
double Solver::jacy(int i, int j)
{
    double Djx = mesh.cell(i,j).center[0] - mesh.cell(i,j-1).center[0];
    double Djy = mesh.cell(i,j).center[1] - mesh.cell(i,j-1).center[1];

    double Dix = 0.25 *(
                 mesh.cell(i+1,j-1).center[0] \
                 + mesh.cell(i+1,j).center[0]\
                 - mesh.cell(i-1, j-1).center[0]\
                 - mesh.cell(i-1, j).center[0]
                 );
    double Diy = 0.25 *(
                 mesh.cell(i+1,j-1).center[1] \
                 + mesh.cell(i+1,j).center[1]\
                 - mesh.cell(i-1, j-1).center[1]\
                 - mesh.cell(i-1, j).center[1]
                 );
    return 1.0/(Dix*Djy - Djx*Diy);
}

double* Solver::vis_x(int i, int j) {
    double T_avg = 0;
    int iprim = i + nGhost;
    int jprim = j+ nGhost;
    if (i == 0) {
        T_avg = field.Temp(iprim, jprim);
    } else if (i == ncx) {
        T_avg = field.Temp(iprim-1, jprim);
    } else {
        T_avg = (field.Temp(iprim, jprim) + field.Temp(iprim-1, jprim))/2.0;
    }

    double mu = mu_t(T_avg);
    double k = k_t(T_avg);

    double *viscous_x = new double [5];

    double Dix = mesh.cell(i,j).center[0] - mesh.cell(i-1,j).center[0];
    double Diy = mesh.cell(i,j).center[1] - mesh.cell(i-1,j).center[1];

    double Djx = 0.25 *(
                 mesh.cell(i-1,j+1).center[0] \
                 + mesh.cell(i,j+1).center[0]\
                 - mesh.cell(i-1, j-1).center[0]\
                 - mesh.cell(i, j-1).center[0]
                 );
    double Djy = 0.25 *(
                 mesh.cell(i-1,j+1).center[1] \
                 + mesh.cell(i,j+1).center[1]\
                 - mesh.cell(i-1, j-1).center[1]\
                 - mesh.cell(i, j-1).center[1]
                 );
    double Ds = 1.0/(Dix*Djy - Djx*Diy);

    double Diu = field.prim[iprim][jprim][1] - field.prim[iprim-1][jprim][1];
    double Div = field.prim[iprim][jprim][2] - field.prim[iprim-1][jprim][2];
    double DiT = field.Temp(iprim, jprim) - field.Temp(iprim-1, jprim);
    double Dju = 0.25 *(
        field.prim[iprim-1][jprim+1][1] \
        + field.prim[iprim][jprim+1][1] \
        - field.prim[iprim-1][jprim-1][1] \
        - field.prim[iprim][jprim-1][1]
    );
    double Djv = 0.25 *(
        field.prim[iprim-1][jprim+1][2] \
        + field.prim[iprim][jprim+1][2] \
        - field.prim[iprim-1][jprim-1][2] \
        - field.prim[iprim][jprim-1][2]
    );
    double DjT = 0.25 *(
        field.Temp(iprim-1, jprim+1) \
        + field.Temp(iprim, jprim+1) \
        - field.Temp(iprim-1, jprim-1) \
        - field.Temp(iprim, jprim-1)
    );
    double ux = (Diu*Djy-Dju*Diy)*Ds;
    double vx = (Div*Djy-Djv*Diy)*Ds;
    double Tx = (DiT*Djy-DjT*Diy)*Ds;
    double uy = (-Diu*Djx+Dju*Dix)*Ds;
    double vy = (-Div*Djx+Djv*Dix)*Ds;
    double Ty = (-DiT*Djx+DjT*Dix)*Ds;

    viscous_x[0] = ((4.0/3.0)*ux-(2.0/3.0)*vy) * mu;
    viscous_x[1] = ((4.0/3.0)*vy-(2.0/3.0)*ux) * mu;
    viscous_x[2] = (uy+vx)*mu;
    double u1= (field.prim[iprim][jprim][1] + field.prim[iprim-1][jprim][1])*0.50;
    double v1= (field.prim[iprim][jprim][2] + field.prim[iprim-1][jprim][2])*0.50;

    viscous_x[3] = u1*viscous_x[0] + v1*viscous_x[2] + k*Tx;
    viscous_x[4] = u1*viscous_x[2] + v1*viscous_x[1] + k*Ty;

    return viscous_x;
}


double* Solver::vis_y(int i, int j) {
    double T_avg = 0;
    int iprim = i + nGhost;
    int jprim = j + nGhost;
    if (j == 0) {
        T_avg = field.Temp(iprim, jprim);
    } else if (j == ncy){
        T_avg = field.Temp(iprim, jprim-1);
    } else {
        T_avg = (field.Temp(iprim, jprim) + field.Temp(iprim, jprim-1))/2.0;
    }
    double mu = mu_t(T_avg);
    double k = k_t(T_avg);
    double *viscous_y = new double [5];
    double Djx = mesh.cell(i,j).center[0] - mesh.cell(i,j-1).center[0];
    double Djy = mesh.cell(i,j).center[1] - mesh.cell(i,j-1).center[1];

    double Dix = 0.25 *(
                 mesh.cell(i+1,j-1).center[0] \
                 + mesh.cell(i+1,j).center[0]\
                 - mesh.cell(i-1, j-1).center[0]\
                 - mesh.cell(i-1, j).center[0]
                 );
    double Diy = 0.25 *(
                 mesh.cell(i+1,j-1).center[1] \
                 + mesh.cell(i+1,j).center[1]\
                 - mesh.cell(i-1, j-1).center[1]\
                 - mesh.cell(i-1, j).center[1]
                 );
    double Ds = 1.0/(Dix*Djy - Djx*Diy);

    double Dju = field.prim[iprim][jprim][1] - field.prim[iprim][jprim-1][1];
    double Djv = field.prim[iprim][jprim][2] - field.prim[iprim][jprim-1][2];
    double DjT = field.Temp(iprim, jprim) - field.Temp(iprim, jprim-1);
    double Diu = 0.25 *(
        field.prim[iprim+1][jprim-1][1] \
        + field.prim[iprim+1][jprim][1] \
        - field.prim[iprim-1][jprim-1][1] \
        - field.prim[iprim-1][jprim][1]
    );
    double Div = 0.25 *(
        field.prim[iprim+1][jprim-1][2] \
        + field.prim[iprim+1][jprim][2] \
        - field.prim[iprim-1][jprim-1][2] \
        - field.prim[iprim-1][jprim][2]
    );
    double DiT = 0.25 *(
        field.Temp(iprim+1, jprim-1) \
        + field.Temp(iprim+1, jprim) \
        - field.Temp(iprim-1, jprim-1) \
        - field.Temp(iprim-1, jprim)
    );
    double ux = (Diu*Djy - Dju*Diy) * Ds;
    double vx = (Div*Djy - Djv*Diy) * Ds;
    double Tx = (DiT*Djy - DjT*Diy) * Ds;

    double uy=(-Diu*Djx+Dju*Dix)*Ds;
    double vy=(-Div*Djx+Djv*Dix)*Ds;
    double Ty=(-DiT*Djx+DjT*Dix)*Ds;

    viscous_y[0] = ((4.0/3.0)*ux-(2.0/3.0)*vy) * mu;
    viscous_y[1] = ((4.0/3.0)*vy-(2.0/3.0)*ux) * mu;
    viscous_y[2] = (uy+vx)*mu;
    double u1= (field.prim[iprim][jprim][1] + field.prim[iprim][jprim-1][1])*0.50;
    double v1= (field.prim[iprim][jprim][2] + field.prim[iprim][jprim-1][2])*0.50;

    viscous_y[3] = u1*viscous_y[0] + v1*viscous_y[2] + k*Tx;
    viscous_y[4] = u1*viscous_y[2] + v1*viscous_y[1] + k*Ty;

    return viscous_y;
}

void Solver::fill_corner() {
    for (int i=0; i<4; i++) {
        field.prim[1][1][i] = field.prim[1][2][i] + field.prim[2][1][i] - field.prim[2][2][i];
        field.prim[1][ncy+nGhost][i] = field.prim[1][ncy+nGhost-1][i] + field.prim[2][ncy+nGhost][i]\
                                      - field.prim[2][ncy+nGhost-1][i];
        field.prim[ncx+nGhost][1][i] = field.prim[ncx+nGhost][2][i] + field.prim[ncx+nGhost-1][1][i]\
                                     - field.prim[ncx+nGhost-1][2][i];
        field.prim[ncx+nGhost][ncy+nGhost][i] = field.prim[ncx+nGhost][ncy+nGhost-1][i] + field.prim[ncx+nGhost-1][ncy+nGhost][i]\
                                     - field.prim[ncx+nGhost-1][ncy+nGhost-1][i];
    }
}