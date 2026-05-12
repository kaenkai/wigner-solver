#include "WignerSolver.hpp"
#include "Poisson1D.hpp"


/********************
 * Wigner function tools
 ********************/

 void WignerSolver::testBTE() {
    // ------------------------------------------------------
    // Set up grid, Boltzmann equation, and system parameters
    // ------------------------------------------------------
    size_t nx = 200, nk = 200;
    double lD = 1000, lC = 0;
    double k_max = 0.4;
    WignerSolver f(nx, lD, lC, nk, k_max);

    f.set_m(1);
    f.set_temp(300);
    f.set_epsilonR(1);

    f.set_uL(0.1/AU_eV);
    f.set_uR(0.1/AU_eV);
    f.set_uBias(0.2/AU_eV);

    f.setBoundCond();

    f.printParam();

    // -------------------------------------------
    // Solve Poisson equation for linear potential
    // -------------------------------------------
    Poisson1D p(f.get_nx(), f.get_dx());
    p.set_boundary_conditions(f.get_uL(), f.get_uR());
    p.set_epsilonR(EPS_GaAs);
    p.solve();
    f.set_uC(p.get_uNew());

    // -----------------------
    // Add rectangular barrier
    // -----------------------
    f.addRectBarr(0.1/AU_eV, 500, 100, 10);

    // ----------------------------------
    // Solve Boltzmann transport equation
    // ----------------------------------

    f.solveBTE();

    f.saveDistFun();

    arma::vec cdX = f.calcCD_X();
    cout << "# Current density: " << f.calcCurrentDensity()*AU_A/AU_cm2 << " A/cm^2" << endl;
    cout << "# Carrier density in x space: " << cdX(0)/AU_cm3 << " cm^-3 at x=0 and " << cdX(nx-1)/AU_cm3 << " cm^-3 at x=L" << endl;
    cout << "# Normalization of distribution function: " << f.calcNorm() << endl;
    cout << "# Min and max of distribution function: " << f.f_.min() << ", " << f.f_.max() << endl;

    arma::field<std::string> header(4);
    arma::mat out_data;
    out_data.insert_cols(0, arma::linspace(0, f.get_l(), f.get_nx())), header(0) = "x [au]";
    out_data.insert_cols(1, f.get_u()*AU_eV), header(1) = "U [eV]";
    out_data.insert_cols(2, f.get_currD()/arma::norm(f.get_currD())), header(2) = "J(x)/||J(x)||";
    out_data.insert_cols(3, f.calcCD_X()/AU_cm3), header(3) = "n [cm^{-3}]";
	// out_data.insert_cols(8, f.get_uB()*AU_eV), header(8) = "U^B [eV]";
	// out_data.insert_cols(9, f.get_uC()*AU_eV), header(9) = "U^C [eV]";
    out_data.save( arma::csv_name("output/test.csv", header) );

 }
