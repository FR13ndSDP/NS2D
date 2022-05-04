#include "./include/Solver.h"
#include "./include/params.h"

void Solver::setBoundary()
{
    field.cons2prim();
    // intlet & wall
    for (int i = 2; i < ncx + nGhost; i++)
    {
        // TODO: constant wall temperature
        // wall
        double rho = field.prim[i][2][0];
        double u = field.prim[i][2][1];
        double v = field.prim[i][2][2];
        double p = (G-1) * (field.U[i][2][3] - 0.5 * rho * (u*u + v*v));
        double T_cell = field.Temp(i, 2);
        double T_ghost = 2 * T_wall - T_cell;
        if (T_ghost < 0.5 * T_cell) T_ghost = 0.5 * T_cell;

        field.U[i][1][0] = p/(R * T_ghost);
        field.U[i][1][1] = -u * field.U[i][1][0];
        field.U[i][1][2] = -v * field.U[i][1][0];
        field.U[i][1][3] = p / (G-1) + 0.5 * field.U[i][1][0] * (u*u + v*v);

        rho = field.prim[i][3][0];
        u = field.prim[i][3][1];
        v = field.prim[i][3][2];
        p = (G-1) * (field.U[i][3][3] - 0.5 * rho * (u*u + v*v));
        T_cell = field.Temp(i, 3);
        T_ghost = 2 * T_wall - T_cell;
        if (T_ghost < 0.5 * T_cell) T_ghost = 0.5 * T_cell;

        field.U[i][0][0] = p/(R * T_ghost);
        field.U[i][0][1] = -u * field.U[i][0][0];
        field.U[i][0][2] = -v * field.U[i][0][0];
        field.U[i][0][3] = p / (G-1) + 0.5 * field.U[i][0][0] * (u*u + v*v);
        // inlet
        field.U[i][ncy+nGhost][0] = field.inlet[0];
        field.U[i][ncy+nGhost][1] = field.inlet[1];
        field.U[i][ncy+nGhost][2] = field.inlet[2];
        field.U[i][ncy+nGhost][3] = field.inlet[3];

        field.U[i][ncy+nGhost+1][0] = field.inlet[0];
        field.U[i][ncy+nGhost+1][1] = field.inlet[1];
        field.U[i][ncy+nGhost+1][2] = field.inlet[2];
        field.U[i][ncy+nGhost+1][3] = field.inlet[3];
    }

    for (int j = 2;  j < ncy + nGhost; j++)
    {
        // down
        field.U[1][j][0] = field.U[2][j][0];
        field.U[1][j][1] = field.U[2][j][1];
        field.U[1][j][2] = field.U[2][j][2];
        field.U[1][j][3] = field.U[2][j][3];

        field.U[0][j][0] = field.U[3][j][0];
        field.U[0][j][1] = field.U[3][j][1];
        field.U[0][j][2] = field.U[3][j][2];
        field.U[0][j][3] = field.U[3][j][3];

        // up
        field.U[ncx + nGhost][j][0] = field.U[ncx + nGhost - 1][j][0];
        field.U[ncx + nGhost][j][1] = field.U[ncx + nGhost - 1][j][1];
        field.U[ncx + nGhost][j][2] = field.U[ncx + nGhost - 1][j][2];
        field.U[ncx + nGhost][j][3] = field.U[ncx + nGhost - 1][j][3];

        field.U[ncx + nGhost + 1][j][0] = field.U[ncx + nGhost - 2][j][0];
        field.U[ncx + nGhost + 1][j][1] = field.U[ncx + nGhost - 2][j][1];
        field.U[ncx + nGhost + 1][j][2] = field.U[ncx + nGhost - 2][j][2];
        field.U[ncx + nGhost + 1][j][3] = field.U[ncx + nGhost - 2][j][3];
    }
}