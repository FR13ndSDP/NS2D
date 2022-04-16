#include "./include/Solver.h"

void Solver::timeAdvance()
{
#ifdef EULER
    computeFlux();
    computeRHS();
    for (int i = nGhost; i < ncx + nGhost; i++)
    {
        for (int j = nGhost; j < ncy + nGhost; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i][j][k] += dt * rhs[i - 2][j - 2][k];
            }
        }
    }
#endif

#ifdef RK2
// stage 1
    field.U2Un();
    computeFlux();
    computeRHS();
    for (int i = nGhost; i < ncx + nGhost; i++)
    {
        for (int j = nGhost; j < ncy + nGhost; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i][j][k] += dt * rhs[i - 2][j - 2][k];
            }
        }
    }
// stage 2
    computeFlux();
    computeRHS();
    for (int i = nGhost; i < ncx + nGhost; i++)
    {
        for (int j = nGhost; j < ncy + nGhost; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i][j][k] = 0.5 * field.Un[i][j][k] + \
                0.5*(field.U[i][j][k] + dt * rhs[i - 2][j - 2][k]);
            }
        }
    }
#endif

#ifdef RK3
// stage 1
    field.U2Un();
    computeFlux();
    computeRHS();
    for (int i = nGhost; i < ncx + nGhost; i++)
    {
        for (int j = nGhost; j < ncy + nGhost; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i][j][k] = field.Un[i][j][k] + dt * rhs[i - 2][j - 2][k];
            }
        }
    }
// stage 2
    computeFlux();
    computeRHS();
    for (int i = nGhost; i < ncx + nGhost; i++)
    {
        for (int j = nGhost; j < ncy + nGhost; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i][j][k] = 0.75 * field.Un[i][j][k] + \
                                   0.25*(field.U[i][j][k] + dt * rhs[i - 2][j - 2][k]);
            }
        }
    }
// stage 3
    field.U2Un();
    computeFlux();
    computeRHS();
    for (int i = nGhost; i < ncx + nGhost; i++)
    {
        for (int j = nGhost; j < ncy + nGhost; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i][j][k] = 1.0/3.0 * field.Un[i][j][k] + \
                                   2.0/3.0 * (field.U[i][j][k] + dt * rhs[i - 2][j - 2][k]);
            }
        }
    }
#endif
}