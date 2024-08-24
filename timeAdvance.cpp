#include "./include/Solver.h"

void Solver::timeAdvance(double time, double stop_time)
{
#ifdef EULER
#ifdef SOD
    setBoundary_sod();
#else
    setBoundary_cylinder();
#endif
    computeFlux();
    computeRHS();
    computeDt(time, stop_time);
    
    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i+nGhost][j+nGhost][k] += -rhs[i][j][k];
            }
        }
    }
#endif

#ifdef RK2
// stage 1
#ifdef SOD
    setBoundary_sod();
#else
    setBoundary_cylinder();
#endif
    field.U2Un();
    computeFlux();
    computeRHS();
    computeDt(time, stop_time);
    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i+nGhost][j+nGhost][k] += -rhs[i][j][k];
            }
        }
    }
// stage 2
#ifdef SOD
    setBoundary_sod();
#else
    setBoundary_cylinder();
#endif
    computeFlux();
    computeRHS();
    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                field.U[i+nGhost][j+nGhost][k] = 0.5 * field.Un[i+nGhost][j+nGhost][k] + \
                0.5*(field.U[i+nGhost][j+nGhost][k] - rhs[i][j][k]);
            }
        }
    }
#endif
}