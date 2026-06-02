#include "lib.hpp"
#include "WignerSolver.hpp"
#include "Poisson1D.hpp"
#include <armadillo>
#include <cmath>


// ------------
// Declarations 
// ------------


void testZeroPotential();
void testWavePacket();
void testLinearPotentialDrop();
void testCurrentVoltageCharacteristic();
void testSingleBarrier();


// -----------
// Definitions 
// -----------


/**
 * Zero potential test for verificating integrals
 * for uBias=0 and no dissipation density function should be independent of x
 * hence int dx f(k) = L*f(k) = L*BC(k)
 */
void testZeroPotential() {
    // ----------------------------
	// System/simulation parameters
	// ----------------------------
    size_t nx = 150, nk = 150;
    double l = 1500;
    double k_max = 0.1;
    WignerSolver f(nx, l, nk, k_max);

    f.set_m(0.067);
    f.set_temp(300);
    f.set_epsilonR(13.1);

    // ---------------------------------
	// Fermi level in left/right contact
    // ---------------------------------
    double muF = calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp());
    f.set_uL( muF + 0./AU::eV );
    f.set_uR( muF );

    // ---------------------------------------------------------------------------------------
    // Boundary conditions
    // 0 -> 0 (closed system, no carrier inflow)
	// 1 -> Supply function, 2:4 -> Supply function convolution
	// -1 -> Gaussian, -2:-4 -> Gauss function convolution
    // where convolution is with Lorentzian (2, -2), Gaussian (3, -3) or Voigt (4, -4) profile
	// ---------------------------------------------------------------------------------------
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
    double l = 2137;
    double k_max = 0.1;
    WignerSolver f(nx, l, nk, k_max);

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
    double l = 2137;
    double k_max = 0.15;
    WignerSolver f(nx, l, nk, k_max);

    f.set_m(0.067);
    f.set_temp(300);
    f.set_epsilonR(13.1);

    double muF = calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp());
    f.set_uL( muF + 0./AU::eV );
    f.set_uR( muF );
    f.setBoundCond(1);

    f.set_uC( arma::linspace(f.get_uL(), f.get_uR(), f.get_nx()) );

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
    double l = 2137;
    double k_max = 0.15;
    WignerSolver f(nx, l, nk, k_max);

    f.set_m(0.067);
    f.set_temp(300);
    f.set_epsilonR(13.1);

    double muF = calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp());
    f.set_uR( muF );

    for (double v: arma::linspace(0, 1, 100)) {
        f.set_uL( muF + v/AU::eV );
        f.setBoundCond(1);
        f.solveBTE();
        std::cout << v << '\t'
            << arma::as_scalar(arma::trapz(
                f.get_x(),
                f.calcCurrentDensity()
            ))/f.get_l() * AU::A/AU::cm2 << std::endl;
    }
}


void testSingleBarrier() {
    size_t nx = 89, nk = 89;
    double l = 2137;
    double k_max = 0.15;
    WignerSolver f(nx, l, nk, k_max);

    f.set_m(0.067);
    f.set_temp(300);
    f.set_epsilonR(13.1);

    double muF = calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp());
    f.set_uL( muF + 0./AU::eV );
    f.set_uR( muF );
    f.setBoundCond(1);

    f.set_uC(
        arma::linspace(
            f.get_uL(), 
            f.get_uR(), 
            f.get_nx()
        )
    );
    f.addRectBarr(0.05/AU::eV, 1000, 200, 2);

    f.solveBTE();
    f.saveDistFun();

    // arma::vec du = f.get_du();
    // // double sig = 0.3;
    // size_t m = 100;
    // // arma::vec kernel = arma::linspace(-10, 10, m).transform([sig](double x){return exp(-x*x/sig/sig/2.)/std::sqrt(2*M_PI)/sig;});
    // arma::vec kernel = arma::linspace(-10, 10, m).transform([](double x){return 1/3.+0*x;});
    // kernel.print();

    arma::mat test;
    test.insert_cols(0, f.get_x());
    test.insert_cols(1, f.get_uC()+f.get_uB());
    test.insert_cols(2, f.calcCD_X()/AU::cm3);
    test.insert_cols(3, f.calcCurrentDensity()*AU::A/AU::cm2);
    test.insert_cols(4, arma::min(f.get_f(), 1));
    test.insert_cols(5, f.get_du());
    test.insert_cols(6, calcFirstDer(f.get_u(), f.get_dx()));
    // test.insert_cols(7, arma::conv(du, kernel, "same"));
    test.save("output/test.out", arma::raw_ascii);

    std::cout << "Current density: " << arma::as_scalar(arma::trapz(f.get_x(), f.calcCurrentDensity()))/f.get_l()*AU::A/AU::cm2 << " A/cm^2" << std::endl;
    std::cout << "Density function max: " << f.get_f().max() << std::endl;
    std::cout << "Density function min: " << f.get_f().min() << ", min/max: " << f.get_f().min()/f.get_f().max() << std::endl;
}


// ----------------------
// Main BTE test function
// ----------------------

 void testBTE() {
    // testZeroPotential();
    // testWavePacket();
    // testLinearPotentialDrop();
    // testCurrentVoltageCharacteristic();
    testSingleBarrier();
 }


/**
 * Boltzmann/Wigner-Poisson test sandbox
 * to verify time step CFL condition is used
 * @see https://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy_condition
 */
void testBoltzmannPoisson() {
    size_t nx = 100, nk = 100;
    double l = 100/AU::nm;
    double k_max = 0.15;
    WignerSolver f(nx, l, nk, k_max);

    f.set_m(0.067);
    f.set_temp(300);
    f.set_epsilonR(13.1);
    
    double muF = calcFermiEn(2E18*AU::cm3, f.get_m(), f.get_temp());
    double uBias = 0.2/AU::eV;
    f.set_uL( muF + uBias );
    f.set_uR( muF );
    f.setBoundCond(1);

    f.setDopingProfile(2e18*AU::cm3, 20/AU::nm, 0.02);
    // f.addRectBarr(0.2/AU::eV, 50/AU::nm, 5/AU::nm, 2);

    Poisson1D p(f.get_nx(), f.get_dx());                            // Setting up Poisson solver
    p.set_boundary_conditions(uBias/2., -uBias/2.); // Dirichlet BC
    p.set_epsilonR(f.get_epsilonR()), p.set_temp(f.get_temp());     // Permittivity and temperature

    double mx = 1;
    double dt = .2E-15/AU::s;
    f.set_dt(dt);

    double cfl = f.get_dt()/f.get_dx()*k_max/f.get_m();
    std::cout << "# Courant–Friedrichs–Lewy (CFL) condition: " << cfl << std::endl;

    size_t veryImportantCounter = 0;
    double maxPotChange;
    arma::vec newRho(f.get_nx(), arma::fill::zeros);
    for (size_t i = 0; i < 10000; ++i) {
        p.set_rho(newRho);
        p.solve();
        f.set_uC( p.get_uNew() );
        // f.solveBTE();
        f.solveTimeDependentBTE();
        newRho = p.get_rho()*(1-mx) + (f.get_nD()-f.calcCD_X())*mx;
        maxPotChange = arma::abs(p.get_du()).max();
        std::cout 
            // << i << '\t'
            << i*dt*AU::s*1e12 << '\t'
            << arma::sum(f.calcCurrentDensity()) * f.get_dx() * AU::A/AU::cm2 << '\t'
            << f.calcNorm() << '\t'
            << maxPotChange << std::endl;
        // right now the criterion is that max potential change has to be lower than 1e-7 eV  
        veryImportantCounter += maxPotChange < 1e-8 ? 1 : 0;
        if (veryImportantCounter > 9) break;
    }

    std::cout << "# Density function max: " << f.get_f().max() << std::endl;
    std::cout << "# Density function min: " << f.get_f().min() << ", min/max: " << f.get_f().min()/f.get_f().max() << std::endl;

    arma::mat test;
    test.insert_cols(0, f.get_x()*AU::nm);
    test.insert_cols(1, (f.get_uC()+f.get_uB())*AU::eV);
    test.insert_cols(2, (f.get_nD()-f.calcCD_X())/AU::cm3);
    test.insert_cols(3, f.get_nD()/AU::cm3);
    test.insert_cols(4, f.calcCD_X()/AU::cm3);
    test.insert_cols(5, (f.calcCurrentDensity())*AU::A/AU::cm2);
    test.save("output/test.out", arma::raw_ascii);

    f.saveDistFun();
}