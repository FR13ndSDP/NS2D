#include "./include/Field.h"
#include "./include/params.h"
#include <fstream>

Field::Field() {
    // field with ghost cell
    U = new double **[ngx];
    Un = new double **[ngx];
    prim = new double **[(ngx)];
    for (int i=0; i< ngx; i++ ) {
        U[i] = new double *[ngy];
        Un[i] = new double *[ngy];
        prim[i] = new double *[ngy];
        for (int j=0; j< ngy; j++) {
            U[i][j] = new double [4];
            Un[i][j] = new double [4];
            prim[i][j] = new double [4];
        }
    }
}

Field::~Field() {
    for (int i=0; i<ngx; i++ ) {
        for (int j=0; j<ngy; j++) {
            delete [] U[i][j];
            delete [] Un[i][j];
            delete [] prim[i][j];
        }
        delete [] U[i];
        delete [] Un[i];
        delete [] prim[i];
    }
    delete [] U;
    delete [] Un;
    delete [] prim;
}

void Field::U2Un() {
    for (int i=0; i<ngx; i++) {
        for (int j=0; j<ngy; j++) {
            for (int k=0; k<4; k++) {
                Un[i][j][k] = U[i][j][k];
            }
        }
    }
}

void Field::init() {
    for (int i=0; i<ngx; i++) {
        for (int j=0; j<ngy; j++) {
            // ic for the cone here
            U[i][j][0] = inlet[0];
            U[i][j][1] = inlet[1];
            U[i][j][2] = inlet[2];
            U[i][j][3] = inlet[3];
        }
    }
    cons2prim();
}

void Field::initFromFile(std::string name) {
    std::ifstream dataFile(name);
    if (!dataFile.is_open()) throw std::runtime_error("could not open file");
    std::string nonesense, line, var;
    int i=0, j=0;

    if (dataFile.good()) {
        // skip first three lines
        std::getline(dataFile, nonesense);
        std::getline(dataFile, nonesense);
        std::getline(dataFile, nonesense);
    }

    // read all the points
    while (!dataFile.eof()) {
        dataFile >> nonesense >> nonesense\
        >> prim[i+2][j+2][0] >> prim[i+2][j+2][1] >> prim[i+2][j+2][2] >> prim[i+2][j+2][3] >> nonesense;
        if (j == ncy-1) {
            j = 0; i++;
        } else {
            j++;
        }
    }
    dataFile.close();
    prim2cons();
}

void Field::cons2prim() {
    for (int i=0; i<ngx; i++) {
        for (int j=0; j<ngy; j++) {
            prim[i][j][0] = U[i][j][0];
            prim[i][j][1] = U[i][j][1]/U[i][j][0];
            prim[i][j][2] = U[i][j][2]/U[i][j][0];
            prim[i][j][3] = (U[i][j][3] - 0.5 * prim[i][j][0] *(\
                             prim[i][j][1]*prim[i][j][1] + prim[i][j][2]*prim[i][j][2])) \
                             * (G-1.0);
        }
    }
}

void Field::prim2cons() {
    for (int i=0; i<ngx; i++) {
        for (int j=0; j<ngy; j++) {
            U[i][j][0] = prim[i][j][0];
            U[i][j][1] = prim[i][j][0] * prim[i][j][1];
            U[i][j][2] = prim[i][j][0] * prim[i][j][2];
            U[i][j][3] = prim[i][j][3]/(G-1.0) \
                         + 0.5*prim[i][j][0]*\
                         (prim[i][j][1]*prim[i][j][1] + prim[i][j][2]*prim[i][j][2]);
        }
    }
}

double Field::Temp(int i, int j) {
    return prim[i][j][3]/(R * prim[i][j][0]);
}