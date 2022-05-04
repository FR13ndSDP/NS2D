#include "./include/Solver.h"
#include <iostream>
#ifdef TIMING
#include <time.h>
#endif

int main(int argc, char *argv[])
{
    Solver s;
    std::string filename = "";
    int step = 1;
    double time = step * dt;
    s.init(filename);
    for (int i = step; i < step + std::stoi(argv[1]); i++)
    {   
        s.timeAdvance();
        if (i%show_interval == 0) {
            std::cout << "step : " << i << " time :" << time << std::endl;
            s.showRHS();
        }
        if (i%plot_interval == 0) {
            s.writeFile("step_"+std::to_string(i)+".dat");
        }
        time += dt;
    }
    s.writeFile("final.dat");
    return 0;
}
