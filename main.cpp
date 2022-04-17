#include "./include/Solver.h"
#include <iostream>
#ifdef TIMING
#include <time.h>
#endif

int main(int argc, char *argv[])
{
#ifdef TIMING
    clock_t start_t, end_t;
    start_t = clock();
#endif
    Solver s;
    std::string filename = "step_10000.dat";
    int step = 10000;
    s.init(filename);
    s.writeFile("initial.dat");
    for (int i = step; i <= step + std::stoi(argv[1]); i++)
    {
        s.timeAdvance();
        if (i%plot_interval == 0) {
            std::cout << "step : " << i << std::endl;
            s.writeFile("step_"+std::to_string(i)+".dat");
        }
    }
    s.writeFile("final.dat");

#ifdef TIMING
    end_t = clock();
    std ::cout << double(end_t - start_t) / CLOCKS_PER_SEC << " s elapsed" << std::endl;
#endif
    return 0;
}
