#include "./include/Solver.h"
#include <string>

Solver::Solver()
{
    rhs = new double **[ncx];
    for (int i = 0; i < ncx; i++)
    {
        rhs[i] = new double *[ncy];
        for (int j = 0; j < ncy; j++)
        {
            rhs[i][j] = new double[4];
        }
    }
};

Solver::~Solver()
{
    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            delete[] rhs[i][j];
        }
        delete[] rhs[i];
    }
    delete[] rhs;
};

void Solver::init(std::string name)
{
    if (name.empty())
    {
        field.init();
    } else {
        field.initFromFile(name);
    }
    mesh.readMesh("mesh2d.dat");
}

void Solver::writeFile(std::string name)
{
    field.cons2prim();
    std::ofstream file(name, std::ios::out);
    file << "TITLE = \"CONE\"" << std::endl;
    file << "VARIABLES = \"x\", \"y\", \"rho\", \"u\", \"v\", \"p\", \"T\"" << std::endl;
    file << "ZONE I=" << ncy << " J=" << ncx << std::endl;
    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            file << mesh.cell(i, j).center.coord_x << " "
                 << mesh.cell(i, j).center.coord_y << " ";
            for (int k = 0; k < 4; k++)
            {
                file << field.prim[i + 2][j + 2][k] << " ";
            }
            file << field.Temp(i+2, j+2) << std::endl;
        }
    }
    file.close();
}

void Solver::computeRHS()
{
    double rhs_max = 0;
    for (int i = 0; i < ncx; i++)
    {
        for (int j = 0; j < ncy; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                rhs[i][j][k] = (f.Fx[i + 1][j][k] - f.Fx[i][j][k] + f.Fy[i][j+1][k] - f.Fy[i][j][k]);
            }
            if (fabs(rhs[i][j][0]) > rhs_max)
                    rhs_max = fabs(rhs[i][j][0]);
        }
    }
    std::cout << "max rhs : " << rhs_max << std:: endl;
}
