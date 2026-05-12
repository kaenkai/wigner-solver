#include "Poisson1D.hpp"

using namespace AtomicUnits;


/**
 * Solves 1D Poisson equation with Dirichlet boundary conditions for an uniform charge density
 *
 * Analytical solution:
 * $\phi(x) = \frac{\rho(x)}{2\epsilon}x(L-x) = \frac{1}{2}x(1-x)$,
 * for $L=1$, $\rho=1$, and $\epsilon=1$.
 */
void Poisson1D::testUniformCharge() {
    const size_t nx = 101; // grid size
    double len = 1;
    double h = len/(nx-1); // step size [AU]
    Poisson1D p(nx, h);

    // ------------------
    // Set charge density
    // ------------------
    double rho = 1;
    p.rho_.fill(-rho);
    p.epsilonR_ = 4*M_PI;
    p.set_boundary_conditions(0, 0);
    p.set_epsilonR(1);

    p.solve();
    arma::vec phi_num = -p.uNew_;

    // -------------------
    // Analytical solution
    // -------------------
    double epsilon = p.epsilonR_/4./M_PI;
    arma::vec phi_an = arma::linspace(0, len, nx)
        .transform(
            [rho, epsilon, len](double x){
                return -rho/2/epsilon*x*(len-x);
            }
        );
    arma::vec phi_err = phi_num-phi_an;

    // -------------
    // Print results
    // -------------
    cout << "# Poisson equation test for uniform charge density" << endl;
    cout << "# Left BC = " << p.get_dirichletL()*AU_eV
         << ", Right BC = " << p.get_dirichletR()*AU_eV
         << ", Grid points = " << p.get_nx() 
         << ", Spacing = " << p.get_h() << endl;
    cout << "# Charge = " << rho << endl;
    cout << "x phi_num phi_an phi_err res" << endl;
    for (size_t i = 1; i < nx-1; ++i) {
        cout << i*h << '\t'
             << phi_num(i) << '\t'
             << phi_an(i) << '\t'
             << phi_err(i) << '\t'
             << (phi_num(i-1)-2*phi_num(i)+phi_num(i+1))/h/h+p.rho_(i)/epsilon << endl;
    }
    cout << "# Potential minimum (numerical): " << arma::min(phi_num) << endl;
    cout << "# Potential minimum (analytical): " << arma::min(phi_an) << endl;
}


/**
 * Solves 1D Poisson equation with Dirichlet boundary conditions for an exponential charge density
 *
 * Analytical solution:
 * $\rho(x) = \epsilon\frac{k\pi}{L}^2\sin(k\pi x/L) = (k\pi)^2\sin(k\pi x)$,
 * $\pxi(x) = \sin(k\pi x/L) = \sin(k\pi x)$,
 * for $L=1$, $\rho=1$, and $\epsilon=1$.
 */
void Poisson1D::testExponentCharge(){
    const size_t nx = 101; // grid size
    double len = 1;
    double h = len/(nx-1); // step size [AU]
    Poisson1D p(nx, h);

    // ------------------
    // Set charge density
    // ------------------
    int alpha = 1;
    p.set_boundary_conditions(1, std::exp(alpha));
    p.set_epsilonR(4*M_PI);
    double epsilon = p.epsilonR_/4./M_PI;
    p.rho_ = arma::linspace(0, len, nx)
        .transform(
            [alpha, len, epsilon](double x) -> double {
                return -epsilon*std::pow(alpha/len, 2)*std::exp(alpha*x/len);
            }
        );
    p.solve();
    arma::vec phi_num = -p.uNew_;

    // -------------------
    // Analytical solution
    // -------------------
    arma::vec phi_an = arma::linspace(0, len, nx)
        .transform(
            [alpha, len](double x) -> double {
                return -std::exp(alpha*x/len);
            }
        );
    arma::vec phi_err = phi_num-phi_an;

    // -------------
    // Print results
    // -------------
    cout << "# Poisson equation test for uniform charge density" << endl;
    cout << "# Left BC = " << p.get_dirichletL()*AU_eV
         << ", Right BC = " << p.get_dirichletR()*AU_eV
         << ", Grid points = " << p.get_nx() 
         << ", Spacing = " << p.get_h() << endl;
    cout << "x phi_num phi_an phi_err res" << endl;
    for (size_t i = 1; i < nx-1; ++i) {
        cout << i*h << '\t'
             << phi_num(i) << '\t'
             << phi_an(i) << '\t'
             << phi_err(i) << '\t'
             << (phi_num(i-1)-2*phi_num(i)+phi_num(i+1))/h/h+p.rho_(i)/epsilon << endl;
    }
    cout << "# Potential minimum (numerical): " << arma::min(phi_num) << endl;
    cout << "# Potential minimum (analytical): " << arma::min(phi_an) << endl;
}


/**
 * Solves 1D Poisson equation with Dirichlet boundary conditions for a sinusoidal charge density
 *
 * Analytical solution:
 * $\rho(x) = \epsilon\frac{k\pi}{L}^2\sin(k\pi x/L) = (k\pi)^2\sin(k\pi x)$,
 * $\pxi(x) = \sin(k\pi x/L) = \sin(k\pi x)$,
 * for $L=1$, $\rho=1$, and $\epsilon=1$.
 * 
 * @return maximum error between numerical and analytical solution
 */
void Poisson1D::testSineCharge(){
    // ------------------
    // Set charge density
    // ------------------
    const size_t nx = 101; // grid size
    double len = 1;
    double h = len/(nx-1); // step size [AU]
    Poisson1D p(nx, h);
    p.set_boundary_conditions(0, 0);
    p.set_epsilonR(4*M_PI);
    double epsilon = p.epsilonR_/4./M_PI;
    int k = 6;
    p.rho_ = arma::linspace(0, len, p.nx_)
        .transform(
            [k, len, epsilon](double x) -> double {
                return -epsilon*std::pow(k*M_PI/len, 2)*std::sin(k*M_PI*x/len);
            }
        );
    p.rho_ += arma::vec(p.nx_, arma::fill::randn)
        * epsilon*std::pow(k*M_PI/len, 2) * 1e-2;  // add random noise

    p.solve();
    arma::vec phi_num = -p.uNew_;

    // -------------------
    // Analytical solution
    // -------------------
    arma::vec phi_an = arma::linspace(0, len, p.nx_)
        .transform(
            [k, len](double x) -> double {
                return -std::sin(k*M_PI*x/len);
            }
        );
    arma::vec phi_err = (phi_num-phi_an);
    phi_err.transform([](double x) -> double {return std::abs(x);});

    // -------------
    // Print results
    // -------------
    cout << "# Poisson equation test for uniform charge density" << endl;
    cout << "# Left BC = " << p.get_dirichletL()*AU_eV
         << ", Right BC = " << p.get_dirichletR()*AU_eV
         << ", Grid points = " << p.get_nx() 
         << ", Spacing = " << p.get_h() << endl;
    cout << "x\trho\tphi_num\tphi_an\tphi_err\tres" << endl;
    for (size_t i = 1; i < nx-1; ++i) {
        cout << i*h << '\t'
             << p.rho_(i) << '\t'
             << phi_num(i) << '\t'
             << phi_an(i) << '\t'
             << phi_err(i) << '\t'
             << (phi_num(i-1)-2*phi_num(i)+phi_num(i+1))/h/h+p.rho_(i)/epsilon << endl;
    }
    cout << "# Potential norm (numerical): " << arma::norm(phi_num) << endl;
    cout << "# Potential norm (analytical): " << arma::norm(phi_an) << endl;
}


/**
 * Solves 1D Poisson equation with Dirichlet boundary conditions for a sinusoidal charge density
 *
 * Analytical solution:
 * $\rho(x) = \epsilon\frac{k\pi}{L}^2\sin(k\pi x/L) = (k\pi)^2\sin(k\pi x)$,
 * $\pxi(x) = \sin(k\pi x/L) = \sin(k\pi x)$,
 * for $L=1$, $\rho=1$, and $\epsilon=1$.
 * 
 * @return maximum error between numerical and analytical solution
 */
double Poisson1D::testSine(){
    // ------------------
    // Set charge density
    // ------------------
    this->set_boundary_conditions(0, 0);
    this->set_epsilonR(4*M_PI);
    double epsilon = epsilonR_/4./M_PI;
    int k = 6;
    double len = (nx_-1)*h_;
    rho_ = arma::linspace(0, len, nx_)
        .transform(
            [k, len, epsilon](double x) -> double {
                return -epsilon*std::pow(k*M_PI/len, 2)*std::sin(k*M_PI*x/len);
            }
        );

    this -> solve();
    arma::vec phi_num = -uNew_;

    // -------------------
    // Analytical solution
    // -------------------
    arma::vec phi_an = arma::linspace(0, len, nx_)
        .transform(
            [k, len](double x) -> double {
                return -std::sin(k*M_PI*x/len);
            }
        );
    arma::vec phi_err = (phi_num-phi_an);
    phi_err.transform([](double x) -> double {return std::abs(x);});

    return arma::max(phi_err);
}


/**
 * Tests the convergence of the numerical solution on a sequence of grids
 * with decreasing step size h, the error is expected to decrease quadratically with h
 * (error ~ h^2) for a second-order method.
 */
void Poisson1D::testGrid() {
    double len = 1;
    Poisson1D p(51, len/(51-1));
    double err_prev = p.testSine(), err=0, h_prev = p.get_h(), h=0;
    arma::uvec grid_sizes = {101, 201, 401, 801};
    for (auto nx : grid_sizes) {
        h = len / (nx - 1);
        p = Poisson1D(nx, h);
        err = p.testSine();
        double ratio = err_prev / err;
        double order = std::log(ratio) / std::log(h_prev / h);
        cout << "nx=" << nx 
             << " h=" << p.h_ 
             << " err=" << err      // should decrease by a factor of 4 when h is halved
             << " ratio=" << ratio  // should be ~4, because error should decrease by a factor of 4 when h is halved
             << " order=" << order  // should be ~2, because the method is second-order accurate
             << " err/h^2=" << err/(h*h)
             << endl;
        err_prev = err;
        h_prev = h;
    }
}


/**
 * Tests the self-consistency of the Poisson solver by comparing the charge density
 * with the second derivative of the potential
*/
void Poisson1D::testSelfConsistency() {
    const size_t nx = 101; // grid size
    double len = 1;
    double h = len/(nx-1); // step size [AU]
    Poisson1D p(nx, h);
    p.set_epsilonR(4*M_PI);
    double epsilon = p.epsilonR_/4./M_PI;
    int k = 6;
    arma::vec rho = arma::linspace(0, len, nx)
        .transform(
            [k, len, epsilon](double x) -> double {
                return -epsilon*std::pow(k*M_PI/len, 2)*std::sin(k*M_PI*x/len);
            }
        );
    p.rho_ = rho + arma::vec(p.nx_, arma::fill::randn)
        * epsilon*std::pow(k*M_PI/len, 2) * 1e-2;  // add random noise
    p.solve();
    arma::vec phi_ref = -p.uNew_;
    arma::vec phi_num;
    for (size_t i = 100000; --i;) {
        p.rho_ = calcSecondDer(p.uNew_, h)*epsilon;
        p.solve();
        phi_num = -p.uNew_;
        cout << "Potential p-norm: " << arma::norm(phi_num, "inf")
             << ", solution consistency (phi_ref-phi): " << arma::norm(phi_ref - phi_num, "inf") << endl;
    }
    cout << "Error after 100_000 iterations: " 
         << arma::norm(phi_ref - phi_num, "inf") << endl;
}


/**
 * Solves 1D Poisson equation with Dirichlet boundary conditions for
 * an uniformly charged infinite plane located at 0.
 */
void Poisson1D::testChargedPlane() {
    // Create grid from -100 to 100 nm with 1nm spacing
    Poisson1D p(201, 1/AU_nm);
    
    // Set sheet charge density to -0.001 C/m^2 
    double sigma = -1E-3 * AU_m2/E0;
    
    // Place the sheet charge at x=0 (middle of grid)
    p.rho_(p.get_nx()/2) = sigma/(p.get_h());

    // Position grid (-100:100 nm)
    arma::vec x = arma::linspace(-100/AU_nm, 100/AU_nm, p.get_nx());
    
    // Set boundary conditions:
    // At x=0nm: V = 0 
    // At x=100nm: V = -5.65 eV
    p.set_boundary_conditions(-5.647/AU_eV, -5.647/AU_eV);
    
    p.solve();
    arma::vec phi_num = -p.uNew_;  // convert to potential by negating uNew_

    // Output results
    cout << "# Sheet charge density = " << sigma*E0/AU_m2 << " C/m^2" 
         << ", Left BC (x = -100 nm) = " << p.get_dirichletL()*AU_eV << " eV"
         << ", Right BC (x = +100 nm) = " << p.get_dirichletR()*AU_eV << " eV"
         << ", Grid points = " << p.get_nx() 
         << ", Spacing = " << p.get_h()*AU_nm << " nm" << endl;
    // Output x [nm], potential [eV], charge density [C/cm^3]
    cout << "x[nm]\trho[C/cm^3]\tV[eV]" << endl;
    for (size_t i = 0; i < p.get_nx(); ++i) {
        cout << x(i)*AU_nm << '\t' << p.rho_(i)*E0/AU_cm3 << '\t' << phi_num(i)*AU_eV << endl;
    }
}
