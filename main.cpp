#include "WignerSolver.hpp"
#include "Poisson1D.hpp"
#include <iostream>
#include <chrono>


/**
 * Main function for running simulations (input if you prefer)
 * @todo create separate functions and clean up
 */
int main(){

    // -------------------------------------------------------------
    // Set up grid, Wigner/Boltzmann equation, and system parameters
    // -------------------------------------------------------------
    // size_t nx = 200, nk = 200;
    // double lD = 1000, lC = 0;
    // double k_max = 0.1;
    // WignerSolver f(nx, lD, lC, nk, k_max);    

	// ---------------------------------------
	// System/simulation parameters
    // Constants are given in src/lib.hpp file
	// ---------------------------------------
    // f.set_m(M_GaAs);
    // f.set_temp(TEMP);
    // f.set_epsilonR(EPS_GaAs);

    // ---------------------------------
	// Fermi level in left/right contact
    // ---------------------------------
    // f.set_uL( calcFermiEn(ND*AU::cm3, f.get_m(), f.get_temp()) );
    // f.set_uR( calcFermiEn(ND*AU::cm3, f.get_m(), f.get_temp()) );
    // f.set_uL(0.1/AU::eV);
    // f.set_uR(0.1/AU::eV);

    // -----------------------------------
	// Miscellaneous simulation parameters
    // -----------------------------------
    // f.set_dt(0.1E-15/AU::s);
    // f.set_useQC(false);      // Quantum correction term (third 'p' derivative)?

	// -----------------
    // Dissipation terms
	// -----------------
    // f.set_scR(0), f.set_scM(0); // 1./(1e-12/AU::s)
    // f.set_scG(0), f.set_scF(0), f.set_lambda(0);

    // --------------------------------------
    // Setting up potential bias and barriers
	// --------------------------------------
    // std::cout<<"# Setting up potential"<<std::endl;
    // f.set_uBias(0.2/AU::eV);
    // Poisson1D p(f.get_nx(), f.get_dx());
    // p.set_boundary_conditions(f.get_uL(), f.get_uR());
    // p.set_epsilonR(EPS_GaAs);
    // p.solve();
    // f.set_uC(p.get_uNew());
    // f.addRectBarr(0.1/AU::eV, 400, 100, 10);
    // f.addRectBarr(0.3/AU::eV, 700, 100, 10);
    // f.addRectBarr(0.3/AU::eV, 1750/AU::nm, 200/AU::nm, 10);
    // f.addRectBarr(0.3/AU::eV, 2250/AU::nm, 200/AU::nm, 10);
    // f.addGaussBarr(0.3/AU::eV, 500/AU::nm, 100/AU::nm);
	// load_poisson_pot: loads potential from binary file to uStart variable
    // f.load_poisson_pot("poisson_pot_100meV_4e4it.bin");
    // f.set_uC(f.get_uStart());

	// ------------------------------------------------
    // Boundary conditions
    // 0 -> 0 (closed system, no carrier inflow)
	// 1 -> Supply function, 2:4 -> Supply function convolution
	// -1 -> Gaussian, -2:-4 -> Gauss function convolution
    // where convolution is with Lorentzian (2, -2), Gaussian (3, -3) or Voigt (4, -4) profile
	// ------------------------------------------------
    // std::cout<<"# Setting up BC"<<std::endl;
    // f.setBoundCond(1);

    // --------------------------------------
    // Setting equilibrium function from file
	// --------------------------------------
	// if no file is given, equilibrium function is calculated from boundary conditions
    // f.setEquilibriumFunction(); --- IGNORE ---

    // Print system parameters
    // f.printParam();

    // ----------------------
    // Start calculation time
    // ----------------------
    auto t_start = std::chrono::steady_clock::now();

    // ----------------------
    // Wigner/Boltzmann test
	// ----------------------
    testBTE();

    // ---------------------
    // Poisson equation test
    // ---------------------
    // std::cout<<"# Solving Poisson equation"<<std::endl;
    // Poisson1D::testUniformCharge();
    // Poisson1D::testExponentCharge();
    // Poisson1D::testSineCharge();
    // Poisson1D::testGrid();
    // Poisson1D::testChargedPlane();
    // Poisson1D::testSelfConsistency();

    //
    // Evaluating calculation time
    //
    auto t_end = std::chrono::steady_clock::now();
    auto t_elapsed =  std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
    std::cout << "# RUN TIME: " << t_elapsed.count() << " ms" << " (" 
              << int(t_elapsed.count()/1000./60.) << " min "
              << int(t_elapsed.count()/1000.)%60 << " s "
              << int(t_elapsed.count())%1000 << " ms)" << std::endl;

    return 0;
}
