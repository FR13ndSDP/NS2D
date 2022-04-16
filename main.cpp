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
    s.init();
    s.setBoundary();
    s.writeFile("initial.dat");
    for (int i = 1; i <= std::stoi(argv[1]); i++)
    {
        s.timeAdvance();
        s.setBoundary();
        if (i%100 == 0) {
            std::cout << "step : " << i << std::endl;
            s.writeFile("step_"+std::to_string(i)+".dat");
        }
    }
    s.writeFile("final.dat");
    s.writeDebug();

#ifdef TIMING
    end_t = clock();
    std ::cout << double(end_t - start_t) / CLOCKS_PER_SEC << " s elapsed" << std::endl;
#endif
    return 0;
}
