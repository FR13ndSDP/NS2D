#include "./include/Solver.h"

void Solver::writeDebug() {
    std::ofstream file("debug.dat", std::ios::out);
    file << "TITLE = \"CONE\"" << std::endl;
    //file << "VARIABLES = \"x\", \"y\", \"varx1\", \"varx2\", \"varx3\", \"varx4\", "
    //     << "\"vary1\", \"vary2\", \"vary3\", \"vary4\""<< std::endl;
    file << "VARIABLES =\"x\", \"y\", \"rhs0\", \"rhs1\", \"rhs2\", \"rhs3\"" << std::endl;
    file << "ZONE I=" << ncy << " J=" << ncx << std::endl;
    for (int i = 0; i < ncx; i++) {
        for (int j=0; j<ncy; j++) {
            file << mesh.cell(i,j).center.coord_x << " "\
                 << mesh.cell(i,j).center.coord_y << " "; 

            //file << f.Fx[i][j][0] << " " << f.Fx[i][j][1] << " "<<  f.Fx[i][j][2]
            //<< " " << f.Fx[i][j][3] << " ";

            //file << f.Fy[i][j][0] << " " << f.Fy[i][j][1] << " "<<  f.Fy[i][j][2]
            //<< " " << f.Fy[i][j][3] << std::endl;
            file << rhs[i][j][0] << " " << rhs[i][j][1] << " "<<  rhs[i][j][2] \
                 << " " << rhs[i][j][3] << std::endl;
        }
    }


    file.close();
}
