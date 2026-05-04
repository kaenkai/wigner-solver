#include "lib.hpp"
#include "Poisson1D.hpp"

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
