#include "Poisson1D.hpp"

/**
 * Wrapper for Poisson equation solver
 * @todo implement more efficient solvers and iterative methods (i.e. Gummel iterati    on)
 * @todo implement von Neumann BC and mixed BC
 * @see Biegel, B. A. & Plummer, J. D. Phys. Rev. B 54, 8070–8082 (1996).
 */
void Poisson1D::solve() { solve_tridiag(); }


/**
 * Solving Poisson equation using tridiagonal matrix
 * warning: matrix equation is solved for potential energy (hence no minus sign in the rhs), not potential
 */
void Poisson1D::solve_tridiag() {

    // --------------------------------
    // Set up discrete Poisson equation
    // --------------------------------

    double epsilon = epsilonR_/4./M_PI;
    arma::sp_mat A(nx_, nx_);
    for (size_t i=nx_; i--;) {
        A(i,i) = -2;
        if (i > 0) A(i,i-1) = 1;
        if (i < nx_-1) A(i,i+1) = 1;
    }
    arma::vec d = rho_*h_*h_/epsilon;  // -rho_*h_*h_/epsilon when solving for potential
    arma::vec x(nx_, arma::fill::zeros);  // A*x = d

    // ------------
    // Dirichlet BC
    // ------------

    A.row(0).zeros();
    A(0,0) = 1;
    d(0) = dirichletL_;

    A.row(nx_-1).zeros();
    A(nx_-1,nx_-1) = 1;
    d(nx_-1) = dirichletR_;

    // ----------------------
    // Solve Poisson equation
    // ----------------------s
    
    arma::superlu_opts settings;
    settings.symmetric = true;
    arma::spsolve(x, A, d, "superlu", settings);
    uNew_ = x;
    du_ = uNew_ - uOld_;
}
