#include "./include/Solver.h"

void Solver::timeAdvance()
{
#ifdef EULER
    setBoundary();
    computeFlux();
    computeRHS();
    computeDt();
    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i+nGhost][j+nGhost][k] += -dt_local[i][j] * rhs[i][j][k];
            }
        }
    }
#endif

#ifdef RK2
// stage 1
    setBoundary();
    field.U2Un();
    computeFlux();
    computeRHS();
    computeDt();
    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i+nGhost][j+nGhost][k] += -dt_local[i][j] * rhs[i][j][k];
            }
        }
    }
// stage 2
    setBoundary();
    computeFlux();
    computeRHS();
    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i+nGhost][j+nGhost][k] = 0.5 * field.Un[i+nGhost][j+nGhost][k] + \
                0.5*(field.U[i+nGhost][j+nGhost][k] - dt_local[i][j] * rhs[i][j][k]);
            }
        }
    }
#endif
}