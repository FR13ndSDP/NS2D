#include "./include/Flux.h"

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

// use primitive variable

void Flux::fluxSplit(double *QL, double *QR, double *Flux)
{
    double dl, uul, vvl, pl, al, hl, dr, uur, vvr, pr, ar, hr;
    double a, Ml, Mr, Mp4, Mm4, Pp5, Pm5, M, mp, mm, p, dm, xm;
    double *fp = new double[4];
    double *fm = new double[4];

    dl = QL[0];
    uul = QL[1];
    vvl = QL[2];
    pl = QL[3];
    dr = QR[0];
    uur = QR[1];
    vvr = QR[2];
    pr = QR[3];
    hl = G * pl / ((G - 1.0) * dl) + 0.5 * (uul * uul + vvl * vvl);
    hr = G * pr / ((G - 1.0) * dr) + 0.5 * (uur * uur + vvr * vvr);

    al = sqrt(G * pl / dl);
    ar = sqrt(G * pr / dr);
    a = 0.5 * (al + ar);
    Ml = uul / a;
    Mr = uur / a;

    if (fabs(Ml) >= 1.0)
    {
        Mp4 = 0.5 * (Ml + fabs(Ml));
    }
    else
    {
        Mp4 = 0.25 * (Ml + 1.0) * (Ml + 1.0) + 0.125 * (Ml * Ml - 1.0) * (Ml * Ml - 1.0);
    }

    if (fabs(Mr) >= 1.0)
    {
        Mm4 = 0.5 * (Mr - fabs(Mr));
    }
    else
    {
        Mm4 = -0.25 * (Mr - 1.0) * (Mr - 1.0) - 0.125 * (Mr * Mr - 1.0) * (Mr * Mr - 1.0);
    }

    M = Mp4 + Mm4;
    mp = dl * a * fmax(0.0, M);
    mm = dr * a * fmin(0.0, M);
    xm = mp + mm;
    if (M >= 0)
    {
        dm = a * fabs(M) * dl;
    }
    else
    {
        dm = a * fabs(M) * dr;
    }

    if (fabs(Ml) >= 1.0)
        Pp5 = 0.5 * (Ml + fabs(Ml)) / Ml;
    else
        Pp5 = 0.25 * (Ml + 1.0) * (Ml + 1.0) * (2.0 - Ml) + 3.0 * Ml * (Ml * Ml - 1.0) * (Ml * Ml - 1.0) / 16.0;
    if (fabs(Mr) >= 1.0)
        Pm5 = 0.5 * (Mr - fabs(Mr)) / Mr;
    else
        Pm5 = 0.25 * (Mr - 1.0) * (Mr - 1.0) * (2.0 + Mr) - 3.0 * Mr * (Mr * Mr - 1.0) * (Mr * Mr - 1.0) / 16.0;

    p = pl * Pp5 + pr * Pm5;

    fp[0] = 1.0;
    fp[1] = uul;
    fp[2] = vvl;
    fp[3] = hl;

    fm[0] = 1.0;
    fm[1] = uur;
    fm[2] = vvr;
    fm[3] = hr;

    Flux[0] = 0.5 * xm * (fp[0] + fm[0]) - 0.5 * dm * (fm[0] - fp[0]);
    Flux[1] = 0.5 * xm * (fp[1] + fm[1]) - 0.5 * dm * (fm[1] - fp[1]) + p;
    Flux[2] = 0.5 * xm * (fp[2] + fm[2]) - 0.5 * dm * (fm[2] - fp[2]);
    Flux[3] = 0.5 * xm * (fp[3] + fm[3]) - 0.5 * dm * (fm[3] - fp[3]);

    delete [] fp;
    delete [] fm;
}
