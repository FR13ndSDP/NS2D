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

void Flux::fluxSplit(double *primL, double *primR, double *Flux)
{
    double *fp = new double[4];
    double *fn = new double[4];

    double rhoL = primL[0], uL = primL[1], vL = primL[2], pL = primL[3];
    double rhoR = primR[0], uR = primR[1], vR = primR[2], pR = primR[3];

    double cL = sqrt(G * pL / rhoL);
    double cR = sqrt(G * pR / rhoR);

    double ML = uL / cL, MR = uR / cR;

    if (ML >= 1.0)
    {
        fp[0] = rhoL * uL;
        fp[1] = rhoL * uL * uL + pL;
        fp[2] = rhoL * uL * vL;
        fp[3] = uL * (G * pL / (G - 1.0) + 0.5 * rhoL * (uL * uL + vL * vL));
    }
    else if (fabs(ML) < 1.0)
    {
        double Mp = 0.25 * (1 + ML) * (1 + ML);
        double tmp = rhoL * cL * Mp;
        fp[0] = tmp;
        fp[1] = tmp * ((G - 1.0) * uL + 2.0 * cL) / G;
        fp[2] = tmp * vL;
        fp[3] = tmp * (pow((G - 1.0) * uL + 2.0 * cL, 2) * 0.5 / (G * G - 1.0) + 0.5 * vL * vL);
    }
    else if (ML <= -1.0)
    {
        fp[0] = 0.0;
        fp[1] = 0.0;
        fp[2] = 0.0;
        fp[3] = 0.0;
    }

    if (MR >= 1.0)
    {
        fn[0] = 0.0;
        fn[1] = 0.0;
        fn[2] = 0.0;
        fn[3] = 0.0;
    }
    else if (fabs(MR) < 1.0)
    {
        double Mn = -0.25 * (MR - 1) * (MR - 1);
        double tmpn = rhoR * cR * Mn;
        fn[0] = tmpn;
        fn[1] = tmpn * ((G - 1.0) * uR - 2.0 * cR) / G;
        fn[2] = tmpn * vR;
        fn[3] = tmpn * (pow((G - 1.0) * uR - 2.0 * cR, 2) * 0.5 / (G * G - 1.0) + 0.5 * vR * vR);
    }
    else if (MR <= -1.0)
    {
        fn[0] = rhoR * uR;
        fn[1] = rhoR * uR * uR + pR;
        fn[2] = rhoR * uR * vR;
        fn[3] = uR * (G * pR / (G - 1.0) + 0.5 * rhoR * (uR * uR + vR * vR));
    }

    for (int i = 0; i < 4; i++)
    {
        Flux[i] = fp[i] + fn[i];
    }

    delete[] fp;
    delete[] fn;
}

// use primitive variable
/*
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
*/

/*
void Flux::fluxSplit(double *QL, double *QR, double *flux)
{
    double dl, uul, vvl, pl, al, dr, uur, vvr, pr, ar;
    double tmp0, tmp1, tmp3, E1P, E2P, E3P, E1M, E2M, E3M;
    double *fp = new double[4], *fm = new double[4];

    dl = QL[0];
    uul = QL[1];
    vvl = QL[2];
    pl = QL[3];
    dr = QR[0];
    uur = QR[1];
    vvr = QR[2];
    pr = QR[3];

    al = sqrt(G * pl / dl);
    ar = sqrt(G * pr / dr);

    tmp1 = 2 * (G - 1.0);
    tmp3 = (3.0 - G) / (2 * (G - 1.0));

    E1P = (uul + fabs(uul)) * 0.5;
    E2P = (uul - al + fabs(uul - al)) * 0.5;
    E3P = (uul + al + fabs(uul + al)) * 0.5;
    tmp0 = dl / (2 * G);
    fp[0] = tmp0 * (tmp1 * E1P + E2P + E3P);
    fp[1] = tmp0 * (tmp1 * E1P * uul + E2P * (uul - al) + E3P * (uul + al));
    fp[2] = tmp0 * (tmp1 * E1P * vvl + E2P * vvl + E3P * vvl);
    fp[3] = tmp0 * (E1P * (G - 1.0) * (uul * uul + vvl * vvl) + E2P * ((uul - al) * (uul - al) + vvl * vvl) * 0.5 + E3P * ((uul + al) * (uul + al) + vvl * vvl) * 0.5 + tmp3 * al * al * (E2P + E3P));

    E1M = (uur - fabs(uur)) * 0.5;
    E2M = (uur - ar - fabs(uur - ar)) * 0.5;
    E3M = (uur + ar - fabs(uur + ar)) * 0.5;
    tmp0 = dr / (2 * G);
    fm[0] = tmp0 * (tmp1 * E1M + E2M + E3M);
    fm[1] = tmp0 * (tmp1 * E1M * uur + E2M * (uur - ar) + E3M * (uur + ar));
    fm[2] = tmp0 * (tmp1 * E1M * vvr + E2M * vvr + E3M * vvr);
    fm[3] = tmp0 * (E1M * (G - 1.0) * (uur * uur + vvr * vvr) + E2M * ((uur - ar) * (uur - ar) + vvr * vvr) * 0.5 + E3M * ((uur + ar) * (uur + ar) + vvr * vvr) * 0.5 + tmp3 * ar * ar * (E2M + E3M));

    for (int i = 0; i < 4; i++)
        flux[i] = fp[i] + fm[i];
        
    delete[] fp;
    delete[] fm;
}
*/