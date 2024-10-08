#include "./include/Solver.h"

#include <string>

Solver::Solver() {
    rhs = new double **[ncx];
    for (int i = 0; i < ncx; i++) {
        rhs[i] = new double *[ncy];
        for (int j = 0; j < ncy; j++) {
            rhs[i][j] = new double[4];
        }
    }

    dt = 1e-8; // initial dt;
};

Solver::~Solver() {
    for (int i = 0; i < ncx; i++) {
        for (int j = 0; j < ncy; j++) {
        delete[] rhs[i][j];
        }
        delete[] rhs[i];
    }
    delete[] rhs;
};

void Solver::init(std::string name) {
    if (name.empty()) {
#ifdef SOD
        field.init_sod();
#else 
        field.init_cylinder();
#endif
    } else {
        field.initFromFile(name);
    }
    mesh.readMesh("mesh2d.dat");
}

void Solver::writeFile(std::string name) {
    field.cons2prim();
    std::ofstream file(name, std::ios::out);
    file << "TITLE = \"CONE\"" << std::endl;
    file << "VARIABLES = \"x\", \"y\", \"rho\", \"u\", \"v\", \"lambda\", \"T\""
         << std::endl;
    file << "ZONE I=" << ncy << " J=" << ncx << std::endl;
    for (int i = 0; i < ncx; i++) {
        for (int j = 0; j < ncy; j++) {
            file << mesh.cell(i, j).center.coord_x << " "
                 << mesh.cell(i, j).center.coord_y << " ";
        for (int k = 0; k < 4; k++) {
            file << field.prim[i + 2][j + 2][k] << " ";
        }
        file << field.Temp(i+2, j+2) << std::endl;
        }
    }
    file.close();
}

void Solver::computeRHS() {
    for (int i = 0; i < ncx; i++) {
        for (int j = 0; j < ncy; j++) {
            for (int k = 0; k < 4; k++) {
                rhs[i][j][k] = (f.Fx[i + 1][j][k] - f.Fx[i][j][k] + f.Fy[i][j + 1][k] -
                        f.Fy[i][j][k])/mesh.cell(i,j).volume;
            }
      }
    }
}

void Solver::showRHS() {
    double rhs_max = 0;
    double rhs_rms = 0;
    
    for (int i = 1; i < ncx-1; i++) {
        for (int j = 1; j < ncy-1; j++) {
            rhs_rms += rhs[i][j][0] * rhs[i][j][0];

            if (fabs(rhs[i][j][0]) > rhs_max) rhs_max = fabs(rhs[i][j][0]);
        }
    }

    rhs_rms = sqrt(rhs_rms) / (nc-2);

    std::cout << "max_rhs: " << rhs_max << std::endl;
    std::cout << "rms_rhs: " << rhs_rms << std::endl;
}


void Solver::computeDt(double time, double stop_time) {
    double mindt = 1e16;

    for (int i = 0; i < ncx; i++) {
        for (int j = 0; j < ncy; j++) {
        double s0 = mesh.cell(i,j).volume;
        double si=0.5 * (mesh.cell(i,j).area[0] + mesh.cell(i,j).area[2]);
        double sj=0.5 * (mesh.cell(i,j).area[1] + mesh.cell(i,j).area[3]);
        Point ni = (mesh.cell(i,j).faceNorm[0] - mesh.cell(i,j).faceNorm[2]);
        Point nj = (mesh.cell(i,j).faceNorm[3] - mesh.cell(i,j).faceNorm[1]);
        double uni = 0.5 *(field.prim[i+nGhost][j+nGhost][1]*ni[0] + field.prim[i+nGhost][j+nGhost][2]*ni[1]);
        double unj = 0.5 *(field.prim[i+nGhost][j+nGhost][1]*nj[0] + field.prim[i+nGhost][j+nGhost][2]*nj[1]);
        double cc = sqrt(G*field.prim[i+nGhost][j+nGhost][3]/field.prim[i+nGhost][j+nGhost][0]);
        double Lci=(fabs(uni)+cc)*si;
        double Lcj=(fabs(unj)+cc)*sj;

        double dt_l=s0/std::max(Lci, Lcj);

        if (dt_l < mindt)
            mindt = dt_l;
        }
    }

    dt = mindt * CFL;

    if (time + dt > stop_time) 
        dt = stop_time - time;
}