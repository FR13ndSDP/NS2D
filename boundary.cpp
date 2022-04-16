#include "./include/Solver.h"
#include "./include/params.h"

void Solver::setBoundary()
{
    double u,v;
    // intlet & wall
    for (int i = 2; i < ncx + nGhost; i++)
    {
        // wall
        u = field.U[i][2][1];
        v = field.U[i][2][2];
        Point n1 = mesh.cell(i-2, 0).faceNorm[3];
        Point vel1(u,v);
        Point vn1 = 2 * (vel1 * n1) * n1;
        vel1 -= vn1;
        field.U[i][1][0] = field.U[i][2][0];
        field.U[i][1][1] = vel1[0];
        field.U[i][1][2] = vel1[1];
        field.U[i][1][3] = field.U[i][2][3];

        u = field.U[i][3][1];
        v = field.U[i][3][2];
        Point n2 = mesh.cell(i-2, 0).faceNorm[3];
        Point vel2(u,v);
        Point vn2 = 2 * (vel2 * n2) * n2;
        vel2 -= vn2;
        field.U[i][0][0] = field.U[i][3][0];
        field.U[i][0][1] = vel2[0];
        field.U[i][0][2] = vel2[1];
        field.U[i][0][3] = field.U[i][3][3];

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