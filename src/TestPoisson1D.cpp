#include "lib.hpp"
#include "Poisson1D.hpp"

using namespace AtomicUnits;


/**
 * Static method for tesing Poisson solver.
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
    // -----------------
    double rho = 1;
    p.rho_.fill(rho);
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
                return rho/2/epsilon*x*(len-x);
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
    cout << "# Potential max value (numerical): " << arma::max(phi_num) << endl;
    cout << "# Potential max value (analytical): " << arma::max(phi_an) << endl;
}

/**
 * Static method for tesing Poisson solver.
 * Solves 1D Poisson equation with Dirichlet boundary conditions for an sinusoidal charge density
 *
 * Analytical solution:
 * $\rho(x) = \epsilon\frac{k\pi}{L}^2\sin(k\pi x/L) = (k\pi)^2\sin(k\pi x)$,
 * $\pxi(x) = \sin(k\pi x/L) = \sin(k\pi x)$,
 * for $L=1$, $\rho=1$, and $\epsilon=1$.
 */
void Poisson1D::testSinus(){
    const size_t nx = 101; // grid size
    double len = 1;
    double h = len/(nx-1); // step size [AU]
    Poisson1D p(nx, h);

    // ------------------
    // Set charge density
    // -----------------
    p.set_boundary_conditions(0, 0);
    p.set_epsilonR(4*M_PI);
    double epsilon = p.epsilonR_/4./M_PI;
    int k = 4;
    p.rho_ = arma::linspace(0, len, nx)
        .transform(
            [k, len, epsilon](double x)->double{
                return epsilon*std::pow(k*M_PI/len, 2)*std::sin(k*M_PI*x/len);
            }
        );

    p.solve();
    arma::vec phi_num = -p.uNew_;

    // -------------------
    // Analytical solution
    // -------------------
    arma::vec phi_an = arma::linspace(0, len, nx)
        .transform(
            [k, len](double x)->double{
                return std::sin(k*M_PI*x/len);
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
    cout << "# Potential max value (numerical): " << arma::max(phi_num) << endl;
    cout << "# Potential max value (analytical): " << arma::max(phi_an) << endl;
}


/**
 * Static method for tesing Poisson solver.
 * Solves 1D Poisson equation with Dirichlet boundary conditions for an sinusoidal charge density
 *
 * Analytical solution:
 * $\rho(x) = \epsilon\frac{k\pi}{L}^2\sin(k\pi x/L) = (k\pi)^2\sin(k\pi x)$,
 * $\pxi(x) = \sin(k\pi x/L) = \sin(k\pi x)$,
 * for $L=1$, $\rho=1$, and $\epsilon=1$.
 */
void Poisson1D::testExponent(){
    const size_t nx = 101; // grid size
    double len = 1;
    double h = len/(nx-1); // step size [AU]
    Poisson1D p(nx, h);

    // ------------------
    // Set charge density
    // -----------------
    int alpha = 1;
    p.set_boundary_conditions(-1, -std::exp(alpha));
    p.set_epsilonR(4*M_PI);
    double epsilon = p.epsilonR_/4./M_PI;
    p.rho_ = arma::linspace(0, len, nx)
        .transform(
            [alpha, len, epsilon](double x)->double{
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
            [alpha, len](double x)->double{
                return std::exp(alpha*x/len);
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
    cout << "# Potential max value (numerical): " << arma::max(phi_num) << endl;
    cout << "# Potential max value (analytical): " << arma::max(phi_an) << endl;
}


/**
 * Static method for tesing Poisson solver.
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

    // Output results
    cout << "# Sheet charge density = " << sigma*E0/AU_m2 << " C/m^2" 
         << ", Left BC (x = -100 nm) = " << p.get_dirichletL()*AU_eV << " eV"
         << ", Right BC (x = +100 nm) = " << p.get_dirichletR()*AU_eV << " eV"
         << ", Grid points = " << p.get_nx() 
         << ", Spacing = " << p.get_h()*AU_nm << " nm" << endl;
    // Output x [nm], potential [eV], charge density [C/cm^3]
    for (size_t i = 0; i < p.get_nx(); ++i) {
        cout << x(i)*AU_nm << ' ' << p.uNew_(i)*AU_eV << ' ' << p.rho_(i)*E0/AU_cm3 << endl;
    }
}
