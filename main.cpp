#include "./include/Solver.h"
#include <iostream>
#ifdef TIMING
#include <time.h>
#endif

int main(int argc, char *argv[])
{
    int stop_step = std::stoi(argv[1]);
    double stop_time = std::stod(argv[2]);

    Solver s;
    std::string filename = "";
    int step = 1;
    double time = 0.0;
    s.init(filename);
    for (int i = step; i <= stop_step; i++)
    {   
        time += s.dt;

        s.timeAdvance(time, stop_time);
        if (i%show_interval == 0) {
            std::cout << "step : " << i << " time :" << time << " dt :"<< s.dt << std::endl;
            s.showRHS();
        }
        if (i%plot_interval == 0) {
            s.writeFile("step_"+std::to_string(i)+".dat");
        }

        if (time >= stop_time) 
            break;
    }
    std::cout << "end time : " << time << std::endl;
    s.writeFile("final.dat");
    return 0;
}
