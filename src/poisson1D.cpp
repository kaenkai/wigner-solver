#include "lib.hpp"
#include "poisson1D.hpp"

/// Wrapper for Poisson equation solver
void Poisson1D::solve() { solve_tridiag(); }  
// solve_gummel solve_tridiag
/// TODO: switching method through a parameter

/**
 * Solving Poisson equation using Gummel algorithm
 * @see Biegel, B. A. & Plummer, J. D. Phys. Rev. B 54, 8070–8082 (1996).
 */
void Poisson1D::solve_gummel() {

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
    uNew_ = uOld_ + x;

}

/// Solving Poisson equation using tridiagonal matrix
void Poisson1D::solve_tridiag() {

    double epsilon = epsilonR_/4./M_PI;

    arma::sp_mat A(nx_, nx_);
    for (size_t i=nx_; i--;) {
        A(i,i) = -2;
        if (i > 0) A(i,i-1) = 1;
        if (i < nx_-1) A(i,i+1) = 1;
    }
    arma::vec d = rho_;
    arma::vec x(nx_, arma::fill::zeros);  // A*x = d

    for (size_t i=nx_; i--;)
        d(i) *= h_*h_/epsilon;

    // Dirichlet BC
    d(0) -= dirichletL_, d(nx_-1) -= dirichletR_;

    // von Neumann BC
    /// TODO: implement von Neumann BC
    // if (pBC_vN_) {
    //     A(1, 1) = -1, A(nx_-1, nx_-1) = -1;
    //     d(0) -= h_*neumannL_, d(nx_-1) += h_*neumannR_;
    // }

    arma::superlu_opts opts;
    opts.symmetric = true;
    opts.equilibrate = false;
    opts.permutation = arma::superlu_opts::COLAMD;
    opts.refine = arma::superlu_opts::REF_EXTRA;
    // opts.allow_ugly  = false;
    opts.pivot_thresh = 0;

    arma::spsolve(x, A, d, "superlu", opts);
    uNew_ = x;
    du_ = uNew_ - uOld_;

}

/**
 * Static method for tesing Poisson solver.
 * Solves 1D Poisson equation with Dirichlet boundary conditions for
 * an uniformly charged infinite plane located at 0.
 */
void Poisson1D::testPoisson() {
    // Create grid from 0 to 100nm with 1nm spacing
    Poisson1D p(201, 1/AU_nm);
    
    // Set sheet charge density to -0.001 C/m^2 
    double sigma = -1E-3 * AU_m2/E0;
    
    // Place the sheet charge at x=0 (middle of grid)
    p.rho_(p.get_nx()/2) = sigma/(p.get_h());
    arma::vec x = arma::linspace(-100/AU_nm, 100/AU_nm, p.get_nx());
    
    // Set boundary conditions:
    // At x=0nm: V = 0 
    // At x=100nm: V = -5.65 eV
    p.set_boundary_conditions(-5.647/AU_eV, -5.647/AU_eV);
    
    p.solve();

    // Output results
    cout << "# Sheet charge density = " << sigma*E0/AU_m2 << " C/m^2" 
         << ", Left BC (x=-100nm) = " << p.get_dirichletL()*AU_eV << " eV"
         << ", Right BC (x=+100nm) = " << p.get_dirichletR()*AU_eV << " eV"
         << ", Grid points = " << p.get_nx() 
         << ", Spacing = " << p.get_h()*AU_nm << " nm" << endl;
    // Output x [nm], potential [eV], charge density [C/cm^3]
    for (size_t i = 0; i < p.get_nx(); ++i) {
        cout << x(i)*AU_nm << ' ' << p.uNew_(i)*AU_eV << ' ' << p.rho_(i)*E0/AU_cm3 << endl;
    }

    // Linear fit results
    arma::vec fit = arma::polyfit(x, p.uNew_, 1);
    cout<<"Linear fit results: "<<fit(0)*AU_eV/AU_nm<<" "<<fit(1)*AU_eV<<endl;

}
