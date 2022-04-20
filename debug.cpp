#include "./include/Solver.h"

void Solver::writeDebug_F() {
    std::ofstream file("debug_F.dat", std::ios::out);
    file << "TITLE = \"CONE\"" << std::endl;
    //file << "VARIABLES = \"x\", \"y\", \"varx1\", \"varx2\", \"varx3\", \"varx4\", "
    //     << "\"vary1\", \"vary2\", \"vary3\", \"vary4\""<< std::endl;
    file << "VARIABLES =\"x\", \"y\", \"x1\", \"x2\", \"x3\", \"x4\", \"y1\", \"y2\", \"y3\", \"y4\"" << std::endl;
    file << "ZONE I=" << ncy << " J=" << ncx << std::endl;
    for (int i = 0; i < ncx; i++) {
        for (int j=0; j<ncy; j++) {
            file << mesh.cell(i,j).center.coord_x << " "\
                 << mesh.cell(i,j).center.coord_y << " "; 

            file << f.Fx[i][j][0] << " " << f.Fx[i][j][1] << " "<<  f.Fx[i][j][2]
            << " " << f.Fx[i][j][3] << " ";

            file << f.Fy[i][j][0] << " " << f.Fy[i][j][1] << " "<<  f.Fy[i][j][2]
            << " " << f.Fy[i][j][3] << std::endl;
        }
    }


    file.close();
}

void Solver::writeDebug_U() {
    std::ofstream file("debug_U.dat", std::ios::out);
    file << "TITLE = \"CONE\"" << std::endl;
    //file << "VARIABLES = \"x\", \"y\", \"varx1\", \"varx2\", \"varx3\", \"varx4\", "
    //     << "\"vary1\", \"vary2\", \"vary3\", \"vary4\""<< std::endl;
    file << "VARIABLES =\"x\", \"y\", \"x1\", \"x2\", \"x3\", \"x4\"" << std::endl;
    file << "ZONE I=" << ncy << " J=" << ncx << std::endl;
    for (int i = 0; i < ncx; i++) {
        for (int j=0; j<ncy; j++) {
            file << mesh.cell(i,j).center.coord_x << " "\
                 << mesh.cell(i,j).center.coord_y << " "; 

            file << field.U[i][j][0] << " " << field.U[i][j][1] << " "<< field.U[i][j][2]
            << " " << field.U[i][j][3] << std::endl;
        }
    }


    file.close();
}

void Solver::writeDebug_G() {
    std::ofstream file("debug_G.dat", std::ios::out);
    file << "TITLE = \"CONE\"" << std::endl;
    //file << "VARIABLES = \"x\", \"y\", \"varx1\", \"varx2\", \"varx3\", \"varx4\", "
    //     << "\"vary1\", \"vary2\", \"vary3\", \"vary4\""<< std::endl;
    file << "VARIABLES =\"x\", \"y\", \"volume\"" << std::endl;
    file << "ZONE I=" << ncy << " J=" << ncx << std::endl;
    for (int i = 0; i < ncx; i++) {
        for (int j=0; j<ncy; j++) {
            file << mesh.cell(i,j).center.coord_x << " "\
                 << mesh.cell(i,j).center.coord_y << " "; 

            file << mesh.cell(i,j).area[0] << std::endl;
        }
    }


    file.close();
}

