#include "lib.hpp"
#include "WignerSolver.hpp"
#include "Poisson1D.hpp"
#include <armadillo>


// ------------
// Declarations 
// ------------


void testZeroPotential();
void testWavePacket();
void testLinearPotentialDrop();
void testCurrentVoltageCharacteristic();


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
    f.set_uL( calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp()) );
    f.set_uR( calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp()) );
    f.set_uBias(0./AU::eV);
    f.setBoundCond(1);

    f.solveBTE();
    f.saveDistFun();

    arma::vec cdX = f.calcCD_X();
    arma::mat cdX_test;
    cdX_test.insert_cols(0, f.get_x());
    cdX_test.insert_cols(1, cdX/AU::cm3);
    cdX_test.insert_cols(2, f.calcCurrentDensity()*AU::A/AU::cm2);
    cdX_test.save("output/cdX_test.out", arma::raw_ascii);

    arma::vec cdK = f.calcCD_K();
    arma::mat cdK_test;
    cdK_test.insert_cols(0, f.get_k());
    cdK_test.insert_cols(1, f.get_bc());
    cdK_test.insert_cols(2, cdK);
    cdK_test.insert_cols(3, f.get_bc()-cdK/f.get_l());
    cdK_test.save("output/cdK_test.out", arma::raw_ascii);

    std::cout << "# Distribution function expected value in X: " << f.calcEX() << std::endl;
    std::cout << "# Distribution function expected value in K: " << f.calcEK() << std::endl;
    std::cout << "# Distribution function normalization: " << f.calcNorm() << std::endl;
    std::cout << "# Carrier density in x space: " << cdX(0)/AU::cm3 << " cm^-3 at x=0 and " << cdX(nx-1)/AU::cm3 << " cm^-3 at x=L" << std::endl;
    std::cout << "# Min and max of distribution function: " << f.get_f().min() << ", " << f.get_f().max() << std::endl;
}


void testWavePacket() {
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
    
    std::cout << "# Gaussian wave packet localization, blur" << std::endl;
    std::cout << "# x-space: " << f.calcEX() << ", " << f.calcSDX() << std::endl;
    std::cout << "# k-space: " << f.calcEK() << ", " << f.calcSDK() << std::endl;
    std::cout << "# GWP normalization: " << f.calcNorm() << std::endl << std::endl;
    
    f.set_dt(.1E-15/AU::s);
    double t = 0, t_total = 3E-14/AU::s;
    std::cout << "# t\tE[x]\tE[k]\tSD[x]\tSD[k]\tE[x]_an\tN\tJ" << std::endl;
    std::cout << "# fs\tau\tau\tau\tau\tau\tau\tau" << std::endl;
    while (t <= t_total) {
        t += f.get_dt();
        f.solveTimeDependentBTE();
        gwp_center += gwp_params(2)/f.get_m()*f.get_dt() ;
        std::cout << t*AU::s*1e15 << '\t'
             << f.calcEX() << '\t' << f.calcEK() << '\t'
             << f.calcSDX() << '\t' << f.calcSDK() << '\t'
             << gwp_center << '\t'
             << f.calcNorm() << '\t'
             << arma::as_scalar(arma::trapz(f.get_x(), f.calcCurrentDensity())) << std::endl;
    }

    f.saveDistFun();
}

void testLinearPotentialDrop() {
    size_t nx = 100, nk = 100;
    double lD = 2137, lC = 0;
    double k_max = 0.15;
    WignerSolver f(nx, lD, lC, nk, k_max);

    f.set_m(0.067);
    f.set_temp(300);
    f.set_epsilonR(13.1);
    f.set_uL( calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp()) );
    f.set_uR( calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp()) );
    f.set_uBias( 1/AU::eV );
    f.setBoundCond(1);

    f.set_uC(
        arma::linspace(
            f.get_uL()+f.get_uBias(), 
            f.get_uR(), 
            f.get_nx()
        )
    );

    f.solveBTE();
    f.saveDistFun();

    arma::mat test;
    test.insert_cols(0, f.get_x());
    test.insert_cols(1, f.get_uC());
    test.insert_cols(2, f.calcCD_X()/AU::cm3);
    test.insert_cols(3, f.calcCurrentDensity()*AU::A/AU::cm2);
    test.insert_cols(4, arma::min(f.get_f(), 1));
    test.save("output/test.out", arma::raw_ascii);

    std::cout << "Current density: " << arma::as_scalar(arma::trapz(f.get_x(), f.calcCurrentDensity()))/f.get_l()*AU::A/AU::cm2 << " A/cm^2" << std::endl;
    std::cout << "Density function max: " << f.get_f().max() << std::endl;
    std::cout << "Density function min: " << f.get_f().min() << ", min/max: " << f.get_f().min()/f.get_f().max() << std::endl;
}


void testCurrentVoltageCharacteristic() {
    size_t nx = 100, nk = 100;
    double lD = 2137, lC = 0;
    double k_max = 0.15;
    WignerSolver f(nx, lD, lC, nk, k_max);

    f.set_m(0.067);
    f.set_temp(300);
    f.set_epsilonR(13.1);
    f.set_uL( calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp()) );
    f.set_uR( calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp()) );

    for (double v: arma::linspace(0, 1, 100)) {
        f.set_uBias( v/AU::eV );
        f.setBoundCond(1);
        f.solveBTE();
        std::cout << v << '\t'
            << arma::as_scalar(arma::trapz(
                f.get_x(),
                f.calcCurrentDensity()
            ))/f.get_l() * AU::A/AU::cm2 << std::endl;
    }
}


// ------------------
// Main test function
// ------------------

 void testBTE() {
    // testZeroPotential();
    // testWavePacket();
    // testLinearPotentialDrop();
    testCurrentVoltageCharacteristic();
 }
