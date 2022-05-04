#include "./include/Solver.h"

double minMod(double var1, double var2) {
    if (var1*var2 > 0) {
        return fabs(var1) > fabs(var2)? var2: var1;
    } else {
        return 0;
    }
}

void Solver::computeFlux() {
    field.cons2prim();
    fill_corner();
    computeFlux_x();
    computeFlux_y();
}

// TODO: nnd bugs, wall condition
// compute flux on face
void Solver::computeFlux_y()
{
    double *UL = new double[4];
    double *UR = new double[4];
    // rho u v p
    double *primL = new double[4];
    double *primR = new double[4];

    double *primL_rotate = new double[4];
    double *primR_rotate = new double[4];

    double *fy = new double[4];
    double ni, nj;
    double A;

    double *vy = new double [5];

    // y direction (0,I)->(0,I+1)
    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < f.nfy; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                // make wall boundary more robust
                UL[k] = field.U[i + 2][j + 1][k] + 0.5 *minMod(field.U[i + 2][j + 2][k]-field.U[i + 2][j + 1][k], field.U[i + 2][j + 1][k]-field.U[i + 2][j][k]);
                UR[k] = field.U[i + 2][j + 2][k] - 0.5 *minMod(field.U[i + 2][j + 2][k]-field.U[i + 2][j + 1][k], field.U[i + 2][j + 3][k]-field.U[i + 2][j+2][k]);
            }

            primL[0] = UL[0];
            primL[1] = UL[1] / UL[0];
            primL[2] = UL[2] / UL[0];
            primL[3] = (G - 1.0) * (UL[3] - 0.5 * primL[0] * (primL[1] * primL[1] + primL[2] * primL[2]));

            primR[0] = UR[0];
            primR[1] = UR[1] / UR[0];
            primR[2] = UR[2] / UR[0];
            primR[3] = (G - 1.0) * (UR[3] - 0.5 * primR[0] * (primR[1] * primR[1] + primR[2] * primR[2]));
            
            if (j != f.nfy - 1)
            {
                ni = -mesh.cell(i, j).faceNorm[3][0];
                nj = -mesh.cell(i, j).faceNorm[3][1];
                A = mesh.cell(i, j).area[3];
            }
            else
            {
                ni = mesh.cell(i, j - 1).faceNorm[1][0];
                nj = mesh.cell(i, j - 1).faceNorm[1][1];
                A = mesh.cell(i, j - 1).area[1];
            }

            primL_rotate[0] = primL[0];
            primL_rotate[1] = primL[1] * ni + primL[2] * nj;
            primL_rotate[2] = -primL[1] * nj + primL[2] * ni;
            primL_rotate[3] = primL[3];

            primR_rotate[0] = primR[0];
            primR_rotate[1] = primR[1] * ni + primR[2] * nj;
            primR_rotate[2] = -primR[1] * nj + primR[2] * ni;
            primR_rotate[3] = primR[3];
            
            f.fluxSplit(primL_rotate, primR_rotate, fy);

            f.Fy[i][j][1] = (fy[1] * ni - fy[2] * nj) * A;
            f.Fy[i][j][2] = (fy[1] * nj + fy[2] * ni) * A;
            f.Fy[i][j][0] = fy[0] * A;
            f.Fy[i][j][3] = fy[3] * A;

            // viscous part
            vy = vis_y(i, j);
            f.Fy[i][j][1] += -(vy[0] * ni + vy[2] * nj) * A;
            f.Fy[i][j][2] += -(vy[2] * ni + vy[1] * nj) * A;
            f.Fy[i][j][3] += -(vy[3] * ni + vy[4] * nj) * A;
        }
    }

    delete[] UL;
    delete[] UR;
    delete[] primL;
    delete[] primR;
    delete[] primL_rotate;
    delete[] primR_rotate;
    delete[] fy;
    delete[] vy;
}

void Solver::computeFlux_x() {
    double *UL = new double[4];
    double *UR = new double[4];
    // rho u v p
    double *primL = new double[4];
    double *primR = new double[4];

    double *primL_rotate = new double[4];
    double *primR_rotate = new double[4];

    double *fx = new double[4];
    double ni, nj;
    double A;

    double *vx = new double [5];
    // x direction (I,0)->(I+1,0)
    for (int i = 0; i < f.nfx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                UL[k] = field.U[i + 1][j + 2][k] + 0.5 *minMod(field.U[i + 2][j + 2][k]-field.U[i + 1][j + 2][k], field.U[i + 1][j + 2][k]-field.U[i][j+2][k]);
                UR[k] = field.U[i + 2][j + 2][k] - 0.5 *minMod(field.U[i + 2][j + 2][k]-field.U[i + 1][j + 2][k], field.U[i + 3][j + 2][k]-field.U[i + 2][j+2][k]);
            }
            primL[0] = UL[0];
            primL[1] = UL[1] / UL[0];
            primL[2] = UL[2] / UL[0];
            primL[3] = (G - 1.0) * (UL[3] - 0.5 * primL[0] * (primL[1] * primL[1] + primL[2] * primL[2]));

            primR[0] = UR[0];
            primR[1] = UR[1] / UR[0];
            primR[2] = UR[2] / UR[0];
            primR[3] = (G - 1.0) * (UR[3] - 0.5 * primR[0] * (primR[1] * primR[1] + primR[2] * primR[2]));

            if (i != f.nfx - 1)
            {
                ni = -mesh.cell(i, j).faceNorm[0][0];
                nj = -mesh.cell(i, j).faceNorm[0][1];
                A = mesh.cell(i, j).area[0];
            }
            else
            {
                ni = mesh.cell(i - 1, j).faceNorm[2][0];
                nj = mesh.cell(i - 1, j).faceNorm[2][1];
                A = mesh.cell(i - 1, j).area[2];
            }
            primL_rotate[0] = primL[0];
            primL_rotate[1] = primL[1] * ni + primL[2] * nj;
            primL_rotate[2] = -primL[1] * nj + primL[2] * ni;
            primL_rotate[3] = primL[3];

            primR_rotate[0] = primR[0];
            primR_rotate[1] = primR[1] * ni + primR[2] * nj;
            primR_rotate[2] = -primR[1] * nj + primR[2] * ni;
            primR_rotate[3] = primR[3];

            f.fluxSplit(primL_rotate, primR_rotate, fx);

            f.Fx[i][j][1] = (fx[1] * ni - fx[2] * nj) * A;
            f.Fx[i][j][2] = (fx[1] * nj + fx[2] * ni) * A;
            f.Fx[i][j][0] = fx[0] * A;
            f.Fx[i][j][3] = fx[3] * A;

            // viscous part
            f.Fx[i][j][1] += -(vx[0] * ni + vx[2] * nj) * A;
            f.Fx[i][j][2] += -(vx[2] * ni + vx[1] * nj) * A;
            f.Fx[i][j][3] += -(vx[3] * ni + vx[4] * nj) * A;
        }
    }

    delete[] UL;
    delete[] UR;
    delete[] primL;
    delete[] primR;
    delete[] primL_rotate;
    delete[] primR_rotate;
    delete[] fx;
    delete[] vx;
}
