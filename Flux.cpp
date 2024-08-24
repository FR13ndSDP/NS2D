#include "./include/Flux.h"
#include <valarray>

using label   = long long int;
using scalar  = double;
using tensor1 = std::valarray<double>;
using tensor2 = std::valarray<std::valarray<double>>;

Flux::Flux()
{
    Fx = new double **[nfx];
    for (int i = 0; i < nfx; i++)
    {
        Fx[i] = new double *[ncy];
        for (int j = 0; j < ncy; j++)
        {
            Fx[i][j] = new double[4];
        }
    }

    Fy = new double **[ncx];
    for (int i = 0; i < ncx; i++)
    {
        Fy[i] = new double *[nfy];
        for (int j = 0; j < nfy; j++)
        {
            Fy[i][j] = new double[4];
        }
    }
}

Flux::~Flux()
{
    for (int i = 0; i < nfx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            delete[] Fx[i][j];
        }
        delete[] Fx[i];
    }

    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < nfy; j++)
        {
            delete[] Fy[i][j];
        }
        delete[] Fy[i];
    }

    delete[] Fx;
    delete[] Fy;
}

void cal_Moment_xi(const double &lambda, tensor1 &Mxi)
{
    Mxi[0] = 1.0;                                     // <\xi^0>
    Mxi[1] = K/(2.0*lambda);                          // <\xi^2>
    Mxi[2] = (3.0*K + K*(K-1.0))/(4*lambda*lambda);   // <\xi^4>
}

void cal_Moment_uPositive(const double *prim, tensor1 &MuL, tensor1 &Mxi)// x<0, left
{
    
    const double U = prim[1], lambda = prim[3];

    MuL[0] = 0.5 * erfc( -sqrt(lambda)* U ) ;
    MuL[1] = U * MuL[0] + 0.5*exp( -lambda*U*U ) / sqrt( M_PI*lambda );

    for (std::size_t i = 2; i < (MuL).size(); ++i)
    {
        MuL[i] = U * MuL[i - 1] + 0.5 * (i - 1) * MuL[i - 2] / lambda;
    }

    cal_Moment_xi(lambda, Mxi);
}

void cal_Moment_uNegative(const double *prim, tensor1 &MuR, tensor1 &Mxi)// x>0, right
{

    const double U = prim[1], lambda = prim[3];

    MuR[0] = 0.5 * erfc( sqrt(lambda)* U ) ;
    MuR[1] = U * MuR[0] - 0.5*exp( -lambda*U*U ) / sqrt( M_PI*lambda );

    for (std::size_t i = 2; i < (MuR).size(); ++i)
    {
        MuR[i] = U * MuR[i - 1] + 0.5 * (i - 1) * MuR[i - 2] / lambda;
    }
   
    cal_Moment_xi(lambda, Mxi);
}

void cal_Moment_v(const double *prim, tensor1 &Mu, tensor1 &Mxi)
{
    
    const double U = prim[2], lambda = prim[3];
    
    Mu[0] = 1.0;                                     // <u^0>
    Mu[1] = U;                                       // <u^1>

    for (std::size_t i = 2; i < (Mu).size(); ++i)
    {
        Mu[i] = U * Mu[i - 1] + 0.5 * (i - 1) * Mu[i - 2] / lambda;
    }

    cal_Moment_xi(lambda, Mxi);
}

tensor1 Moment_half(const tensor1 &slope, label m, label n, label k, const tensor1 &Mu, const tensor1 &Mv, const tensor1 &Mxi)
{
    tensor1 moment(4);

    moment[0] = slope[0] * Mu[m] * Mv[k] * Mxi[n] \
                + slope[1] * Mu[m+1] * Mv[k] * Mxi[n] \
                + slope[2] * Mu[m] * Mv[k+1] * Mxi[n] \
                + slope[3] * 0.5 * (Mu[m+2] * Mv[k] * Mxi[n] + Mu[m] * Mv[k+2] * Mxi[n] + Mu[m] * Mv[k] * Mxi[n+1]);

    moment[1] = slope[0] * Mu[m+1] * Mv[k] * Mxi[n] \
                + slope[1] * Mu[m+2] * Mv[k] * Mxi[n] \
                + slope[2] * Mu[m+1] * Mv[k+1] * Mxi[n] \
                + slope[3] * 0.5 * (Mu[m+3] * Mv[k] * Mxi[n] + Mu[m+1] * Mv[k+2] * Mxi[n] + Mu[m+1] * Mv[k] * Mxi[n+1]);

    moment[2] = slope[0] * Mu[m] * Mv[k+1] * Mxi[n] \
                + slope[1] * Mu[m+1] * Mv[k+1] * Mxi[n] \
                + slope[2] * Mu[m] * Mv[k+2] * Mxi[n] \
                + slope[3] * 0.5 * (Mu[m+2] * Mv[k+1] * Mxi[n] + Mu[m] * Mv[k+3] * Mxi[n] + Mu[m] * Mv[k+1] * Mxi[n+1]);

    moment[3] = 0.5 * (slope[0] * (Mu[m+2] * Mv[k] * Mxi[n] + Mu[m] * Mv[k+2] * Mxi[n] + Mu[m] * Mv[k] * Mxi[n+1]) \
                + slope[1] * (Mu[m+3] * Mv[k] * Mxi[n] + Mu[m+1] * Mv[k+2] * Mxi[n] + Mu[m+1] * Mv[k] * Mxi[n+1]) \
                + slope[2] * (Mu[m+2] * Mv[k+1] * Mxi[n] + Mu[m] * Mv[k+3] * Mxi[n] + Mu[m] * Mv[k+1] * Mxi[n+1]) \
                + slope[3] * 0.5 * (Mu[m+4] * Mv[k] * Mxi[n] + Mu[m] * Mv[k+4] * Mxi[n] + Mu[m] * Mv[k] * Mxi[n+2] \
                + 2.0 * (Mu[m+2] * Mv[k+2] * Mxi[n] + Mu[m+2] * Mv[k] * Mxi[n+1] + Mu[m] * Mv[k+2] * Mxi[n+1])));

    return moment;
}

void Flux::fluxSplit(double *primL, double *primR, double *Flux, double dt)
{
    tensor1 slope = {1, 0, 0, 0};
    
    double rhoL = primL[0];
    double rhoR = primR[0];

    tensor1 MuL (7), MuR (7);
    tensor1 MvL (7), MvR (7);
    tensor1 MxiL(3), MxiR(3);

    // first order part
    cal_Moment_uPositive(primL, MuL, MxiL);
    cal_Moment_uNegative(primR, MuR, MxiR);
    cal_Moment_v(primL, MvL, MxiL);
    cal_Moment_v(primR, MvR, MxiR);

    tensor1 t(2), flux(4);
    t[0] = dt;
    t[1] = - 0.5 * dt * dt;

    flux = dt * ( rhoL * Moment_half(slope, 1, 0, 0, MuL, MvL, MxiL) + rhoR * Moment_half(slope, 1, 0, 0, MuR, MvR, MxiR));

    Flux[0] = flux[0];
    Flux[1] = flux[1];
    Flux[2] = flux[2];
    Flux[3] = flux[3];
}
