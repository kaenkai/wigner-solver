#include "lib.hpp"
#include "WignerSolver.hpp"
#include "Poisson1D.hpp"
#include <armadillo>
#include <cmath>


// ------------
// Declarations 
// ------------


void testZeroPotential();
void testStationaryGWP();
void testNonstationaryGWP();


// -----------
// Definitions 
// -----------


/**
 * Zero potential test for verificating integrals
 * for uBias=0 and no dissipation density function should be independent of x
 * hence int dx f(k) = L*f(k) = L*BC(k)
 */
void testZeroPotential() {
    size_t nx = 150, nk = 150;
    double lD = 1500, lC = 0;
    double k_max = 0.1;
    WignerSolver f(nx, lD, lC, nk, k_max);

    f.set_m(0.067);
    f.set_temp(300);
    f.set_epsilonR(13.1);
    f.set_uL( calcFermiEn(2E18*AU_cm3, f.get_m(), f.get_temp()) );
    f.set_uR( calcFermiEn(2E18*AU_cm3, f.get_m(), f.get_temp()) );
    f.set_uBias(0./AU_eV);
    f.setBoundCond(1);

    f.solveBTE();
    f.saveDistFun();

    arma::vec cdX = f.calcCD_X();
    arma::mat cdX_test;
    cdX_test.insert_cols(0, f.get_x());
    cdX_test.insert_cols(1, cdX/AU_cm3);
    cdX_test.insert_cols(2, f.calcCurrentDensity()*AU_A/AU_cm2);
    cdX_test.save("output/cdX_test.out", arma::raw_ascii);

    arma::vec cdK = f.calcCD_K();
    arma::mat cdK_test;
    cdK_test.insert_cols(0, f.get_k());
    cdK_test.insert_cols(1, f.get_bc());
    cdK_test.insert_cols(2, cdK);
    cdK_test.insert_cols(3, f.get_bc()-cdK/f.get_l());
    cdK_test.save("output/cdK_test.out", arma::raw_ascii);

    cout << "# Distribution function expected value in X: " << f.calcEX() << endl;
    cout << "# Distribution function expected value in K: " << f.calcEK() << endl;
    cout << "# Distribution function normalization: " << f.calcNorm() << endl;
    cout << "# Carrier density in x space: " << cdX(0)/AU_cm3 << " cm^-3 at x=0 and " << cdX(nx-1)/AU_cm3 << " cm^-3 at x=L" << endl;
    cout << "# Min and max of distribution function: " << f.get_f().min() << ", " << f.get_f().max() << endl;
}


void testStationaryGWP() {
    size_t nx = 150, nk = 100;
    double lD = 2137, lC = 0;
    double k_max = 0.1;
    WignerSolver f(nx, lD, lC, nk, k_max);

    f.set_m(0.067);
    f.set_temp(300);
    f.set_epsilonR(13.1);

    f.addWavePacket(666.666, 67., 0.049, 0.00411);
    f.saveDistFun();

    arma::mat test;
    test.insert_cols(0, f.get_k());
    test.insert_cols(1, f.calcCD_K());
    test.save("output/test.out", arma::raw_ascii);

    cout << "# Gaussian wave packet localization, blur" << endl;
    cout << "# x-space: " << f.calcEX() << ", " << f.calcSDX() << endl;
    cout << "# k-space: " << f.calcEK() << ", " << f.calcSDK() << endl;
    cout << "# GWP normalization: " << f.calcNorm() << endl;
}


void testNonstationaryGWP() {
    size_t nx = 150, nk = 100;
    double lD = 2137, lC = 0;
    double k_max = 0.1;
    WignerSolver f(nx, lD, lC, nk, k_max);

    f.set_m(0.067);
    f.set_temp(300);
    f.set_epsilonR(13.1);

    f.setBoundCond(0);
    arma::vec gwp_params = {666.666, 67., 0.049, 0.00411}; 
    f.addWavePacket(gwp_params(0), gwp_params(1), gwp_params(2), gwp_params(3));
    double gwp_center = gwp_params(0);
    
    cout << "# Gaussian wave packet localization, blur" << endl;
    cout << "# x-space: " << f.calcEX() << ", " << f.calcSDX() << endl;
    cout << "# k-space: " << f.calcEK() << ", " << f.calcSDK() << endl;
    cout << "# GWP normalization: " << f.calcNorm() << endl << endl;
    
    f.set_dt(.1E-15/AU_s);
    double t = 0, t_total = 3E-14/AU_s;
    cout << "# t\tE[x]\tE[k]\tSD[x]\tSD[k]\tE[x]_an\tN\tJ" << endl;
    cout << "# fs\tau\tau\tau\tau\tau\tau\tau" << endl;
    while (t <= t_total) {
        t += f.get_dt();
        f.solveTimeDependentBTE();
        gwp_center += gwp_params(2)/f.get_m()*f.get_dt() ;
        cout << t*AU_s*1e15 << '\t'
             << f.calcEX() << '\t' << f.calcEK() << '\t'
             << f.calcSDX() << '\t' << f.calcSDK() << '\t'
             << gwp_center << '\t'
             << f.calcNorm() << '\t'
             << arma::as_scalar(arma::trapz(f.get_x(), f.calcCurrentDensity())) << endl;
    }

    f.saveDistFun();
}


// ------------------
// Main test function
// ------------------

 void testBTE() {
    // testZeroPotential();
    // testStationaryGWP();
    testNonstationaryGWP();
    // -------------------------------------------
    // Solve Poisson equation for linear potential
    // -------------------------------------------
    // Poisson1D p(f.get_nx(), f.get_dx());
    // p.set_boundary_conditions(f.get_uL(), f.get_uR());
    // p.set_epsilonR(f.get_epsilonR());
    // p.solve();
    // f.set_uC(p.get_uNew());

    // arma::mat force;
    // force.insert_cols(0, f.get_x());
    // force.insert_cols(1, f.get_uB());
    // force.insert_cols(2, f.get_uC());
    // force.insert_cols(3, f.get_u());
    // force.insert_cols(4, f.get_du());
    // force.print();

    // -----------------------
    // Add rectangular barrier
    // -----------------------
    // f.addRectBarr(0.1, 500, 100, 2);

    // ----------------------------------
    // Solve Boltzmann transport equation
    // ----------------------------------

    // arma::field<std::string> header(4);
    // arma::mat out_data;
    // out_data.insert_cols(0, arma::linspace(0, f.get_l(), f.get_nx())), header(0) = "x [au]";
    // out_data.insert_cols(1, f.get_u()*AU_eV), header(1) = "U [eV]";
    // out_data.insert_cols(2, f.get_currD()/arma::norm(f.get_currD())), header(2) = "J(x)/||J(x)||";
    // out_data.insert_cols(3, f.calcCD_X()/AU_cm3), header(3) = "n [cm^{-3}]";
	// // out_data.insert_cols(8, f.get_uB()*AU_eV), header(8) = "U^B [eV]";
	// // out_data.insert_cols(9, f.get_uC()*AU_eV), header(9) = "U^C [eV]";
    // out_data.save( arma::csv_name("output/test.csv", header) );

 }
