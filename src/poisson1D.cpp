#include "lib.hpp"
#include "poisson1D.hpp"

void Poisson1D::solve() { solve_tridiag(); }  // solve_gummel solve_tridiag

void Poisson1D::solve_gummel() {
    // Solving Poisson equation using Gummel algorithm

    double epsilon = epsilonR_/4./M_PI;
    double c = h_*h_/epsilon;

    // P_i

    double phi_L = uOld_(1)-2*uOld_(0)+dirichletL_;
    double phi_R = uOld_(nx_-2)-2*uOld_(nx_-1)+dirichletR_;

    arma::vec uD(nx_, arma::fill::zeros);
    for (size_t i=1; i<nx_-1; ++i)
        uD(i) = uOld_(i+1)-2.*uOld_(i)+uOld_(i-1);
    uD(0) = phi_L, uD(nx_-1) = phi_R;

    pFun_.zeros();
    for (size_t i=0; i<nx_; ++i)
        pFun_(i) = -(uD(i) - c*rho_(i));  // '-' because in equation A*x =-pFun

    // dP_i/du_j

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

    du_ = x*beta_;
    uNew_ = uOld_ + du_;

}


void Poisson1D::solve_tridiag() {
    // Solving Poisson equation using tridiagonal matrix

    double epsilon = epsilonR_/4./M_PI;

    arma::sp_mat A(nx_, nx_);
    for (size_t i=nx_; i--;) {
        A(i,i) = -2;
        if (i > 0) A(i,i-1) = 1;
        if (i < nx_-1) A(i,i+1) = 1;
    }

    // arma::vec b(nx_, arma::fill::ones); b.fill(-2.);  // Diagonal
    // arma::vec c(nx_, arma::fill::ones);  // Upper diagonal
    // arma::vec a(nx_, arma::fill::ones);  // Lower diagonal
    arma::vec d = rho_;
    arma::vec x(nx_, arma::fill::zeros);  // A*x = d

    for (size_t i=nx_; i--;)
        d(i) *= h_*h_/epsilon;

    // Dirichlet BC
    if (pBC_D_)
        d(0) -= dirichletL_, d(nx_-1) -= dirichletR_;

    // von Neumann BC
    // @TODO: Fix it!
    if (pBC_vN_) {
        A(1, 1) = -1, A(nx_-1, nx_-1) = -1;
        d(0) -= h_*neumannL_, d(nx_-1) += h_*neumannR_;
    }

    arma::superlu_opts opts;
    opts.symmetric = true;
    opts.equilibrate = false;
    opts.permutation = arma::superlu_opts::COLAMD;
    opts.refine = arma::superlu_opts::REF_EXTRA;
    // opts.allow_ugly  = false;
    opts.pivot_thresh = 0;

    arma::spsolve(x, A, d, "superlu", opts);

    // arma::vec r, s;
    // r = A*x-d, s = arma::abs(A)*arma::abs(x)+arma::abs(d);
    // double berr = max(abs(r)/s);
    // cout<<"BERR = "<<berr<<endl;

    // c(0) /= b(0);
    // d(0) /= b(0);
    //
    // double m;
    // for (size_t i=1; i<nx_; ++i) {
    //     m = 1./(b(i) - a(i)*c(i-1));
    //     c(i) *= m;
    //     d(i) = (d(i) - a(i)*d(i-1)) * m;
    // }
    //
    // x(nx_-1) = d(nx_-1);
    // for (size_t i=nx_-2; i--;)
    //     x(i) = d(i) - c(i)*x(i+1);
    // x(0) = d(0) - c(0)*x(1);

    // uNew_ = (1-beta_)*uOld_ + beta_*x;  // mixing old and new potential
    uNew_ = x;
    du_ = uNew_ - uOld_;

}


void Poisson1D::testPoisson() {

    double sigma = -1e-3; // 10^-3 C/m^2

    rho_(nx_/2) = sigma/(h_*AU_m)*1e-6;
    rho_(nx_/2) *= AU_cm3/E0;  // SI -> AU

    dirichletL_ = -5.64717/AU_eV;
    dirichletR_ = -5.64717/AU_eV;
    epsilonR_ = 1;

    solve();

    cout<<"# sigma = "<<sigma<<" dirichletL "<<dirichletL_<<" dirichletR "<<dirichletR_<<" n "<<nx_<<" h "<<h_<<endl;
    for (size_t i = 0; i < nx_; ++i)
        cout<<i*h_*AU_nm<<' '<<uNew_(i)*AU_eV<<' '<<rho_(i)*E0/AU_cm3<<endl;

    // Run example:
    // Poisson1D p(200, 1./AU_nm);
    // p.testPoisson();

}
