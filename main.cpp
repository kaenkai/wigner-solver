#include "src/lib.hpp"
#include "src/WignerFunction.hpp"
#include "src/poisson1D.hpp"

#include <chrono>

using namespace AtomicUnits;

int main(){
    size_t nx = 200, nk = 200;
    double lD = 1000/AU_nm, lC = 250/AU_nm;
    double k_max = 0.05;  // -1, 0.15

    WignerFunction f(nx, lD, lC, nk, k_max);

    arma::vec x_val = f.get_x_arr(), k_val = f.get_k_arr();

	//
	// System/simulation parameters
    // Constants are given in src/lib.hpp file
	//
    f.set_m(M_GaAs);
    f.set_temp(TEMP);
    f.set_epsilonR(EPS_GaAs);
    f.set_cD(ND*AU_cm3); // Carrier density in [AU]
	// Fermi level in left/right contact
    f.set_uL( calcFermiEn(f.get_cD(), f.get_m(), f.get_temp()) );
    f.set_uR( calcFermiEn(f.get_cD(), f.get_m(), f.get_temp()) );
	// Miscellaneous simulation parameters
    f.set_dt(0.1E-15/AU_s);
    f.set_useQC(false);  // Quantum correction term (third 'p' derivative)?
    f.set_useNLP(false);  // Calculations with non-local potential?
    f.set_uBias_BC(true);  // Voltage bias given through BC?

	//
    // Dissipation terms
	//
    f.set_rR(0), f.set_rM(0); // 1./(1e-12/AU_s)
    f.set_rG(0), f.set_rF(0), f.set_lambda(0);

    //
	// Setting up differentional scheme
    // Schemes implemented: "CDS1", "UDS1", "UDS2", "UDS3", "HDS22"
	//
    f.set_diffSch_K("UDS2");
    f.set_diffSch_P("UDS2");
    f.set_diffSch_J("UDS2");

	//
    // Boundary conditions
	//
    // 0 -> 0 (closed system, no carrier inflow)
	// 1 -> SF, 2:4 -> SF convolution
	// -1 -> Gauss, -2:-4 -> Gauss function convolution
	//
    cout<<"# Setting up BC"<<endl;
    f.set_bcType(1);

    //
    // Setting up potential bias and barriers
	//
    cout<<"# Setting up potential"<<endl;
    double u_bias = 0.0/AU_eV;
    f.set_uBias(u_bias);
    // f.setPotBias(0.1/AU_eV);
    // f.addRectBarr(0.3/AU_eV, 300, 100, 10);
    // f.addRectBarr(0.3/AU_eV, 700, 100, 10);
    // f.addRectBarr(0.3/AU_eV, 1750/AU_nm, 200/AU_nm, 10);
    // f.addRectBarr(0.3/AU_eV, 2250/AU_nm, 200/AU_nm, 10);
    // f.addGaussBarr(0.3/AU_eV, 500/AU_nm, 100/AU_nm);
    //
	// load_poisson_pot: loads potential from binary file to uStart variable
    // f.load_poisson_pot("poisson_pot_100meV_4e4it.bin");
    // f.set_uC(f.get_uStart());

    //
    // Setting equilibrium function from file
	//
	// if false Wigner/Boltzmann is solved for 0 bias and no dissipation
    // f.setEquilibriumFunction("OutData/wf_feq_BP.bin", true);

    // Print siulation parameters
    f.printParam();

    //
    // Start calculation time
    //
    auto t_start = std::chrono::steady_clock::now();

    //
    // Boltzmann-Poisson test
	//
    /*
    cout<<"## Solving BTE"<<endl;
    f.solveWignerEq();
    f.calcCD_X();
    f.set_doping_profile(0.02);
    //
    Poisson1D p(f.get_nx(), f.get_dx());
    p.set_epsilonR(f.get_epsilonR()), p.set_temp(f.get_temp());  // Permittivity and temperature
    p.set_boundary_conditions(u_bias/2., -u_bias/2.);  // Dirichlet BC
    p.rho_ = f.get_nD() - f.get_cdX();
    p.solve();
    f.set_uC(p.uNew_);
    cout<<"## BTE done"<<endl;
    */

	//
    // Wave packet time evolution
	//
	/*
    f.addWavePacket(500/AU_nm, 100/AU_nm, 0.05, 0.005);  // sqrt(2*f.get_m()*f.get_uL())
    f.addWavePacket(3500/AU_nm, 100/AU_nm, -0.05, 0.005);  // sqrt(2*f.get_m()*f.get_uL())
    double t_total = 1000e-15/AU_s, t = 0, dt = f.get_dt();
    while (t <= t_total) {
        t += dt;
        f.solveTimeEv();
        f.saveWignerFun();
        cout<<t*AU_s*1e15<<' '<<f.calcEK()<<' '<<sqrt(f.calcEK2())<<' '<<f.calcEX()*AU_nm<<endl;
    }
	*/

    //
    // Poisson equation test
    //
    cout<<"# Solving Poisson equation"<<endl;
    Poisson1D::testUniformCharge();

    //
    // Boltzmann-Poisson
	//
	/*
    cout<<"# Solving B-P set of equations"<<endl;
    (uBias, alpha, beta, n_max, timeDependent)
    alpha - density mixing, beta - potential mixing
    f.solveWignerPoisson(0.1/AU_eV, 1E-3, 1, 1, false);
    f.saveWignerFun();
	*/

    //
    // Schrödinger equation
	//
    // cout<<"# Solving Schrödinger equation"<<endl;
    // f.solveSchrEq();

    //
    // Printing and saving results to files located it "OutData" folder
    //
    f.saveWignerFun();
    f.saveTest();
    f.printResults();

    //
    // Evaluating calculation time
    //
    auto t_end = std::chrono::steady_clock::now();
    auto t_elapsed =  std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
    cout<<"# RUN TIME: "<<t_elapsed.count()<<" ms";
    cout<<" ("<<int(t_elapsed.count()/1000./60.)<<" min "<<int(t_elapsed.count()/1000.)%60<<" s "<<int(t_elapsed.count())%1000<<" ms)"<<endl;

    return 0;
}
