#include "lib.hpp"
#include <armadillo>
#include <cmath>
#include "poisson1D.hpp"

using namespace AtomicUnits;


/**
 * Wrapper for Poisson equation solver, available solvers: solve_gummel, solve_tridiag
 * TODO: switching method through a parameter
 */
void Poisson1D::solve() { solve_tridiag(); }  


/**
 * Solving Poisson equation using Gummel algorithm
 * @see Biegel, B. A. & Plummer, J. D. Phys. Rev. B 54, 8070–8082 (1996).
 * @todo fix why solution different then for solve_tridiag()
 */
void Poisson1D::solve_gummel() {

    double epsilon = epsilonR_/4./M_PI;
    double c = h_*h_/epsilon;

    // double phi_L = uOld_(1)-2*uOld_(0)+dirichletL_;
    // double phi_R = uOld_(nx_-2)-2*uOld_(nx_-1)+dirichletR_;

    arma::vec uD(nx_, arma::fill::zeros);
    for (size_t i=1; i<nx_-1; ++i)
        uD(i) = uOld_(i+1)-2.*uOld_(i)+uOld_(i-1);
    uD(0) = dirichletL_, uD(nx_-1) = dirichletR_;

    pFun_.zeros();
    for (size_t i=0; i<nx_; ++i)
        pFun_(i) = -(uD(i) - c*rho_(i));  // '-' because in equation A*x =-pFun

    dPu_ = arma::sp_mat(nx_, nx_);
    for (size_t i=0; i<nx_; ++i) {
        for (size_t j=0; j<nx_; ++j) {
            if (i-1 == j || i+1 == j)
                dPu_(i, j) = 1;
            else if (j == i)
                dPu_(i, j) = -2;  // - c*nE_(i)/temp_;
        }
    }

    arma::vec x(nx_);
    arma::superlu_opts settings;
    settings.symmetric = true;
    // settings.refine = superlu_opts::REF_EXTRA;
    arma::spsolve(x, dPu_, pFun_, "superlu", settings);
    uNew_ = uOld_ + x;

}


/**
 * Solving Poisson equation using tridiagonal matrix
 * matrix equation is solved for potential energy
 */
void Poisson1D::solve_tridiag() {

    double epsilon = epsilonR_/4./M_PI;

    arma::sp_mat A(nx_, nx_);
    for (size_t i=nx_; i--;) {
        A(i,i) = -2;
        if (i > 0) A(i,i-1) = 1;
        if (i < nx_-1) A(i,i+1) = 1;
    }
    arma::vec d = rho_*h_*h_/epsilon;  // -rho_*h_*h_/epsilon when solving for potential
    arma::vec x(nx_, arma::fill::zeros);  // A*x = d

    // Dirichlet BC
    // Left boundary
    A.row(0).zeros();
    A(0,0) = 1;
    d(0) = dirichletL_;

    // Right boundary
    A.row(nx_-1).zeros();
    A(nx_-1,nx_-1) = 1;
    d(nx_-1) = dirichletR_;

    A.raw_print();

    // von Neumann BC
    /// TODO: implement von Neumann BC
    // if (pBC_vN_) {
    //     A(1, 1) = -1, A(nx_-1, nx_-1) = -1;
    //     d(0) -= h_*neumannL_, d(nx_-1) += h_*neumannR_;
    // }
    
    arma::superlu_opts settings;
    settings.symmetric = true;
    arma::spsolve(x, A, d, "superlu", settings);
    uNew_ = x;
    du_ = uNew_ - uOld_;
}


/**
 * Static method for tesing Poisson solver.
 * Solves 1D Poisson equation with Dirichlet boundary conditions for an uniform charge density
 *
 * Analytical solution:
 * $\phi(x) = \frac{\rho(x)}{2\epsilon}(Lx-x^2)$
 */
void Poisson1D::testUniformCharge() {
    size_t nx = 101;    // grid size
    double len = 1;
    double h = len/(nx-1); // step size [AU]
    Poisson1D p(nx, h);

    // Set uniform charge density to -0.0001 AU
    double rho = 1;
    p.rho_.fill(rho);
    p.epsilonR_ = 4*M_PI;
    p.set_boundary_conditions(0, 0);
    p.set_epsilonR(1);

    p.solve();
    arma::vec phi_num = -p.uNew_;

    // Analytical solution
    double epsilon = p.epsilonR_/4./M_PI;
    arma::vec phi_an = arma::linspace(0, len, nx)
        .transform([rho, epsilon, len](double x){return rho/2/epsilon*x*(len-x);});
    arma::vec phi_err = p.uNew_-phi_an;

    // Output results
    cout << "# Poisson equation test for uniform charge density" << endl;
    cout << "# Left BC = " << p.get_dirichletL()*AU_eV
         << ", Right BC = " << p.get_dirichletR()*AU_eV
         << ", Grid points = " << p.get_nx() 
         << ", Spacing = " << p.get_h() << endl;
    cout << "# Charge = " << rho << endl;
    // Output x [nm], charge density [AU, C/cm^3], potential_numeric [V], potential_analytical [V]
    cout << "x rhs phi_num phi_an phi_err phi_num/phi_an res" << endl;
    for (size_t i = 1; i < p.get_nx()-1; ++i) {
        cout << i*h << ' ' <<
                -p.rho_(i)*h*h/epsilon << ' ' <<
                phi_num(i) << ' ' <<
                phi_an(i) << ' ' <<
                phi_err(i) << ' ' <<
                phi_num(i)/phi_an(i) << ' ' <<
                (phi_num(i-1)-2*phi_num(i)+phi_num(i+1))/h/h+p.rho_(i)/epsilon << endl;
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
