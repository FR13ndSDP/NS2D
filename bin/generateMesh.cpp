#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

int main()
{
    double L = 0.3;
    double H = 1.0;
    std::ofstream file("mesh2d.dat", std::ios::out);
    file << "VARIABLES = \"x\", \"y\"" << std::endl;
    file << "ZONE I= 240"<< " J=80" << std::endl;
    for (int j = 0; j < 80; j++)
    {
        for (int i = 0; i < 240; i++)
        {
            file << i*(L/(240-1)) << " "
                 << j*(H/(80-1)) << " " << std::endl;
        }
    }
    file.close();
    return 0;
}