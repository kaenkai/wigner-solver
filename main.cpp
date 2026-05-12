#include "lib.hpp"
#include "WignerSolver.hpp"
#include "Poisson1D.hpp"

#include <armadillo>
#include <chrono>
#include <iostream>

using namespace AtomicUnits;


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
    // f.set_uL( calcFermiEn(ND*AU_cm3, f.get_m(), f.get_temp()) );
    // f.set_uR( calcFermiEn(ND*AU_cm3, f.get_m(), f.get_temp()) );
    // f.set_uL(0.1/AU_eV);
    // f.set_uR(0.1/AU_eV);

    // -----------------------------------
	// Miscellaneous simulation parameters
    // -----------------------------------
    // f.set_dt(0.1E-15/AU_s);
    // f.set_useQC(false);      // Quantum correction term (third 'p' derivative)?

	// -----------------
    // Dissipation terms
	// -----------------
    // f.set_rR(0), f.set_rM(0); // 1./(1e-12/AU_s)
    // f.set_rG(0), f.set_rF(0), f.set_lambda(0);

    // --------------------------------------
    // Setting up potential bias and barriers
	// --------------------------------------
    // cout<<"# Setting up potential"<<endl;
    // f.set_uBias(0.2/AU_eV);
    // Poisson1D p(f.get_nx(), f.get_dx());
    // p.set_boundary_conditions(f.get_uL(), f.get_uR());
    // p.set_epsilonR(EPS_GaAs);
    // p.solve();
    // f.set_uC(p.get_uNew());
    // f.addRectBarr(0.1/AU_eV, 400, 100, 10);
    // f.addRectBarr(0.3/AU_eV, 700, 100, 10);
    // f.addRectBarr(0.3/AU_eV, 1750/AU_nm, 200/AU_nm, 10);
    // f.addRectBarr(0.3/AU_eV, 2250/AU_nm, 200/AU_nm, 10);
    // f.addGaussBarr(0.3/AU_eV, 500/AU_nm, 100/AU_nm);
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
    // cout<<"# Setting up BC"<<endl;
    // f.setBoundCond();

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

    WignerSolver::testBTE();

	// --------------------------
    // Wave packet time evolution
    // move to a test function
	// --------------------------
	/*
    f.addWavePacket(500/AU_nm, 100/AU_nm, 0.05, 0.005);  // sqrt(2*f.get_m()*f.get_uL())
    f.addWavePacket(3500/AU_nm, 100/AU_nm, -0.05, 0.005);  // sqrt(2*f.get_m()*f.get_uL())
    double t_total = 1000e-15/AU_s, t = 0, dt = f.get_dt();
    while (t <= t_total) {
        t += dt;
        f.solveTimeEv();
        f.saveDistFun();
        cout<<t*AU_s*1e15<<' '<<f.calcEK()<<' '<<sqrt(f.calcEK2())<<' '<<f.calcEX()*AU_nm<<endl;
    }
	*/

    // ---------------------
    // Poisson equation test
    // ---------------------
    // cout<<"# Solving Poisson equation"<<endl;
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
    cout<<"# RUN TIME: "<<t_elapsed.count()<<" ms";
    cout<<" ("<<int(t_elapsed.count()/1000./60.)<<" min "<<int(t_elapsed.count()/1000.)%60<<" s "<<int(t_elapsed.count())%1000<<" ms)"<<endl;

    return 0;
}
