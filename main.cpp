#include "WignerSolver.hpp"
#include "Poisson1D.hpp"
#include <iostream>
#include <chrono>


int main(){
    auto t_start = std::chrono::steady_clock::now();

    // ----------------------
    // Wigner/Boltzmann test
	// ----------------------
    // testBTE();
    testBoltzmannPoisson();

    // ---------------------------
    // Evaluating calculation time
    // ---------------------------
    auto t_end = std::chrono::steady_clock::now();
    auto t_elapsed =  std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
    std::cout << "# RUN TIME: " << t_elapsed.count() << " ms" << " (" 
              << int(t_elapsed.count()/1000./60.) << " min "
              << int(t_elapsed.count()/1000.)%60 << " s "
              << int(t_elapsed.count())%1000 << " ms)" << std::endl;

    return 0;
}
