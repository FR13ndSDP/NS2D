#include "./include/Solver.h"
#include <iostream>
#ifdef TIMING
#include <time.h>
#endif

int main(int argc, char *argv[])
{
    Solver s;
    std::string filename = "step_5000.dat";
    int step = 5000;
    s.init(filename);

    for (int i = step; i < step + std::stoi(argv[1]); i++)
    {   
        if (i%plot_interval == 0) {
            std::cout << "step : " << i << std::endl;
            s.writeFile("step_"+std::to_string(i)+".dat");
        }
        s.timeAdvance();
    }
    s.writeFile("final.dat");
    return 0;
}
