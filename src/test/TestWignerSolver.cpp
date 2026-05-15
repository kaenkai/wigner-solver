#include "WignerSolver.hpp"
#include "Poisson1D.hpp"
#include <armadillo>


// ------------
// Declarations 
// ------------

void testNormalization();
void testExpectedValues();


// -----------
// Definitions 
// -----------


/**
 * Zero potential test for verificating integrals
 * for uBias=0 and no dissipation density function should be independent of x
 * hence int dx f(k) = L*f(k) = L*BC(k)
 */
void testNormalization() {
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
    f.normalization();

    f.saveDistFun();

    arma::vec cdX = f.calcCD_X();
    arma::vec cdK = f.calcCD_K();
    arma::mat cdK_test;
    cdK_test.insert_cols(0, f.get_k());
    cdK_test.insert_cols(1, f.get_bc());
    cdK_test.insert_cols(2, cdK);
    cdK_test.insert_cols(3, f.get_bc()-cdK/f.get_l());  // for uBias =0 
    cdK_test.save("output/cdK_test.out", arma::raw_ascii);

    cout << "# Distribution function expected value in X: " << f.calcEX() << endl;
    cout << "# Distribution function expected value in K: " << f.calcEK() << endl;
    cout << "# Distribution function normalization: " << f.calcNorm() << endl;
    cout << "# Current density: " << f.calcCurrentDensity()*AU_A/AU_cm2 << " A/cm^2" << endl;
    cout << "# Carrier density in x space: " << cdX(0)/AU_cm3 << " cm^-3 at x=0 and " << cdX(nx-1)/AU_cm3 << " cm^-3 at x=L" << endl;
    cout << "# Min and max of distribution function: " << f.get_f().min() << ", " << f.get_f().max() << endl;
}


void testExpectedValues() {
    size_t nx = 150, nk = 150;
    double lD = 1000, lC = 0;
    double k_max = 0.1;
    WignerSolver f(nx, lD, lC, nk, k_max);

    f.addWavePacket(333.5, 20., 0.05, 0.005);
    
    f.saveDistFun();

    cout << "# Wavepacket expected value in X: " << f.calcEX() << endl;
    cout << "# Wavepacket expected value in K: " << f.calcEK() << endl;
    cout << "# Wavepacket normalization: " << f.calcNorm() << endl;
}


// ------------------
// Main test function
// ------------------

 void testBTE() {
    testNormalization();
    // testExpectedValues();
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
