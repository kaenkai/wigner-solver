#include "lib.hpp"
#include "WignerSolver.hpp"


/** 
 * Solve stationary Wigner equation
 * @todo OpenMP parallel calculations
 */
void WignerSolver::solveBTE() {
    u_ = uB_ + uC_;
    du_ = calcFirstDer(u_, dx_);
    d3u_ = calcThirdDer(u_, dx_);
    a_.zeros(), b_.zeros();

    // #pragma omp parallel for collapse(2) shared(a_, b_)
    for (size_t i=0; i<nx_; ++i) {
        for (size_t j=0; j<nk_; ++j) {
            diffusionTerm(i, j, -1);
            driftTerm(i, j, -1);
            scatteringTerm(i, j, -1);
            if (useQC_) quantumCorrTerm(i, j, -1);
        }  // end j loop
    }  // end i loop

    arma::vec x(nxk_, arma::fill::zeros);

    arma::spsolve(x, a_, b_, "superlu");

    f_.zeros();
    for (size_t i=0; i<nx_; ++i)
        for (size_t j=0; j<nk_; ++j)
            f_(i,j) = x(i*nk_+j);
}


/** 
 * Solve stationary Wigner equation
 * @todo review armadillo options for sparse matrix solvers
 * @todo OpenMP parallel calculations
 */
void WignerSolver::solveWTE() {
    u_ = uB_ + uC_;
    du_ = calcFirstDer(u_, dx_);
    d3u_ = calcThirdDer(u_, dx_);
    a_.zeros(), b_.zeros();

    // #pragma omp parallel for collapse(2) shared(a_, b_)
    for (size_t i=0; i<nx_; ++i) {
        for (size_t j=0; j<nk_; ++j) {
            diffusionTerm(i, j, -1);
            nonLocalPotentialTerm(i, j);
            scatteringTerm(i, j, -1);
        }  // end j loop
    }  // end i loop

    // Setting up solver options
    // arma::superlu_opts opts;
    // opts.symmetric = true;
    // opts.equilibrate = false;
    // opts.permutation = arma::superlu_opts::COLAMD;
    // opts.refine = arma::superlu_opts::REF_EXTRA;
    // opts.allow_ugly  = false;
    // opts.pivot_thresh = 0;

    arma::vec x(nxk_, arma::fill::zeros);
    arma::spsolve(x, a_, b_, "superlu");

    f_.zeros();
    for (size_t i=0; i<nx_; ++i)
        for (size_t j=0; j<nk_; ++j)
            f_(i,j) = x(i*nk_+j);
}


/** 
 * Solve time dependent BTE
 * @todo review armadillo options for sparse matrix solvers
 * @todo move solution to solveBTE() and solveWTE()
 * @deprecated this function is not used anymore
 */
void WignerSolver::solveTimeDependentBTE() {
    // setEquilibriumFunction();
    u_ = uB_ + uC_;
    du_ = calcFirstDer(u_, dx_);
    d3u_ = calcThirdDer(u_, dx_);
    a_.zeros(), b_.zeros();

    // #pragma omp parallel for collapse(2)
    for (size_t i=nx_; i--;) {
        for (size_t j=nk_; j--;) {
            a_(i*nk_+j, i*nk_+j) += 1.;
            b_(i*nk_+j) = 2.*f_(i, j);
            diffusionTerm(i, j, dt_);
            driftTerm(i, j, dt_);
            scatteringTerm(i, j, dt_);
            if (useQC_) quantumCorrTerm(i, j, dt_);
        }  // end j loop
    }  // end i loop

    // arma::superlu_opts opts;
    // opts.symmetric = true;
    // opts.equilibrate = false;
    // opts.permutation = arma::superlu_opts::COLAMD;
    // opts.refine = arma::superlu_opts::REF_EXTRA;
    // opts.allow_ugly  = false;
    // opts.pivot_thresh = 0;

    arma::vec x(nxk_, arma::fill::zeros);
    arma::spsolve(x, a_, b_, "superlu");

    for (size_t i=0; i<nx_; ++i)
        for (size_t j=0; j<nk_; ++j)
            f_(i,j) = x(i*nk_+j) - f_(i,j);
}


/**
 * Diffusion term, HDS22 rule is used
 * HDS22 = (alpha*CDS2+beta*UDS2)/(alpha+beta)
 * @param i, j grid point indices
 * @param dt time step, if dt < 0, diffusion term is stationary, time dependent otherwise
 */
void WignerSolver::diffusionTerm(size_t i, size_t j, double dt) {
    size_t r = i*nk_ + j;
    int alpha = 2, beta = 1;
    double C = k_(j)/m_/dx_/(alpha+beta)/2. * (dt > 0 ? dt/2. : 1);
    if (k_(j) < 0.) {
        if (i == 0) {
            a_(r, r)             += -3.*k_(j)/m_/dx_/2.;
            a_(r, (i+1)*nk_ + j) += 4.*k_(j)/m_/dx_/2.;
            a_(r, (i+2)*nk_ + j) += -k_(j)/m_/dx_/2.;
        }
        else if (i == nx_-1) {
            a_(r, (i-1)*nk_ + j) += -alpha*C;
            a_(r, r)             += -3.*beta*C;
            b_(r)                += -(alpha+3.*beta)*C*bc_(j);
        }
        else if (i == nx_-2) {
            a_(r, (i-1)*nk_ + j) += -alpha*C;
            a_(r, r)             += -3.*beta*C;
            a_(r, (i+1)*nk_ + j) += (alpha+4.*beta)*C;
            b_(r)                += beta*C*bc_(j);
        }
        else {
            a_(r, (i-1)*nk_ + j) += -alpha*C;
            a_(r, r)             += -3.*beta*C;
            a_(r, (i+1)*nk_ + j) += (alpha+4.*beta)*C;
            a_(r, (i+2)*nk_ + j) += -beta*C;
        }
    }
    else if (k_(j) > 0.) {
        if (i == nx_-1) {
            a_(r, r)             += 3.*k_(j)/m_/dx_/2.;
            a_(r, (i-1)*nk_ + j) += -4.*k_(j)/m_/dx_/2.;
            a_(r, (i-2)*nk_ + j) += k_(j)/m_/dx_/2.;
        }
        else if (i == 0) {
            a_(r, (i+1)*nk_ + j) += alpha*C;
            a_(r, r)             += 3.*beta*C;
            b_(r)                += (alpha+3.*beta)*C*bc_(j);
        }
        else if (i == 1) {
            a_(r, (i+1)*nk_ + j) += alpha*C;
            a_(r, r)             += 3.*beta*C;
            a_(r, (i-1)*nk_ + j) += -(alpha+4.*beta)*C;
            b_(r)                += -beta*C*bc_(j);
        }
        else {
            a_(r, (i+1)*nk_ + j) += alpha*C;
            a_(r, r)             += 3.*beta*C;
            a_(r, (i-1)*nk_ + j) += -(alpha+4.*beta)*C;
            a_(r, (i-2)*nk_ + j) += beta*C;
        }
    }
}


/** 
 * Drift term, UDS1 is used
 * @param i, j grid point indices
 * @param dt time step, if dt <= 0 the term is stationary, time dependent otherwise
 * @todo consider other boundary condition 
 */
void WignerSolver::driftTerm(size_t i, size_t j, double dt) {
    size_t r = i*nk_ + j;
    double F = -du_(i);
    double C = F/dk_ * (dt > 0 ? dt/2. : 1);
    if (F > 0) {
        a_(r, r) += C;
        if (j > 0)
            a_(r, r-1) += -C;
    }
    else if (F < 0) {
        a_(r, r) += -C;
        if (j < nk_-1)
            a_(r, r+1) += C;
    }
}


/**
 * Non-local potential in WTE
 * @param i, j grid point indices
 * @todo verify and implement time dependence
 */
void WignerSolver::nonLocalPotentialTerm(size_t i, size_t j) {
    size_t r = i*nk_ + j;
    size_t v;
    double sum, u1, u2;
    double C;
    // #pragma omp parallel for
    for (size_t l=0; l<nk_; l++){
        v = i*nk_+l;
        if (v <= r) {
            sum = 0;
            for (size_t q=0; q<nk2_; q++){
                if (i+q > nx_-1)
                    u1 = u_(nx_-1);
                else
                    u1 = u_(i+q);
                if ( int(i-q) < 0)
                    u2 = u_(0);
                else
                    u2 = u_(i-q);
                sum += sin_(j,l*nk2_+q) * (u1 - u2);  // sin(2*M_PI/nk_*q*(j-l))
            }  // end q loop
            C = 2./float(nk_) * sum;
            if (v < r){  // To avoid taking diagonal term twice
                a_(r, v) += C;
                a_(v, r) += -C;
            }
            else
                a_(r, r) += C;
        }  // end if
    }  // end l loop
}


/// Scattering term
void WignerSolver::scatteringTerm(size_t i, size_t j, double dt) {
    size_t r = i*nk_ + j;
    double cR = scR_, cM = scM_, cL = lambda_/dk_/dk_;
    if (dt > 0) cR *= dt/2., cM *= dt/2., cL *= dt/2.;
    // #################### scR term ####################
    a_(r, r) += cR;
    b_(r) += feq_(i, j)*cR;
    // #################### scM term ####################
    a_(r, r) += cM;
    a_(r, i*nk_+(nk_-j-1)) += -cM;
    // #################### Lambda term ####################
    a_(r, r) += 2*cL;
    if (j==0)
        a_(r, r+1) += -cL;
    else if (j==nk_-1)
        a_(r, r-1) += -cL;
    else{
        a_(r, r-1) += -cL;
        a_(r, r+1) += -cL;
    }
    // #################### gamma term ####################
    double C = -scF_;
    if (dt > 0) C *= dt/2.;
    double F = -du_(i);  // classical force equal to -du/dx
    // UDS1
    // C *= k_(j)/dk_;
    // a_(r,r) += C;
    // if (F > 0){
    //     if (j == 0)
    //         a_(r, r) += C;
    //     else{
    //         a_(r, r) += C;
    //         a_(r, r-1) += -C;
    //     }
    // }
    // else if (F <= 0){
    //     if (j == nk_-1)
    //         a_(r, r) += -C;
    //     else{
    //         a_(r, r) += -C;
    //         a_(r, r+1) += C;
    //     }
    // }
    // HDS22
    C *= k_(j)/dk_;
    a_(r,r) += C;
    double alpha = 2., beta = 1.;
    double D = C/(alpha+beta);
    if (F <= 0) {
        if (j==0) {
            b_(r) += fermiDirac(kmax_)*alpha*D;
            a_(r, r) += -3.*beta*D;
            a_(r, r+1) += (alpha+4.*beta)*D;
            a_(r, r+2) += -beta*D;
            // a_(r, r) += -3.*C;
            // a_(r, r+1) += 4.*C;
            // a_(r, r+2) += -C;
        }
        else if (j==nk_-1) {
            a_(r, r-1) += -alpha*D;
            a_(r, r) += -3.*beta*D;
            b_(r) += -fermiDirac(kmax_)*(alpha+3.*beta)*D;
            // a_(r, r) += -3.*C;
            // b_(r) += -3.*fermiDirac(kmax_)*C;
        }
        else if (j==nk_-2) {
            a_(r, r-1) += -alpha*D;
            a_(r, r) += -3.*beta*D;
            a_(r, r+1) += (alpha+4.*beta)*D;
            b_(r) += fermiDirac(kmax_)*beta*D;
            // a_(r, r) += -3.*C;
            // a_(r, r+1) += 4.*C;
            // b_(r) += fermiDirac(kmax_)*C;
        }
        else {
            a_(r, r-1) += -alpha*D;
            a_(r, r) += -3.*beta*D;
            a_(r, r+1) += (alpha+4.*beta)*D;
            a_(r, r+2) += -beta*D;
        }
    }
    else if (F > 0) {
        if (j==nk_-1) {
            b_(r) += -fermiDirac(-kmax_)*alpha*D;
            a_(r, r) += 3.*beta*D;
            a_(r, r-1) += -(alpha+4.*beta)*D;
            a_(r, r-2) += beta*D;
            // a_(r, r) += 3.*C;
            // a_(r, r-1) += -4.*C;
            // a_(r, r-2) += C;
        }
        else if (j==0) {
            a_(r, r+1) += alpha*D;
            a_(r, r) += 3.*beta*D;
            b_(r) += fermiDirac(-kmax_)*(alpha+3.*beta)*D;
            // a_(r, r) += 3.*C;
            // b_(r) += 3.*fermiDirac(-kmax_)*C;
        }
        else if (j==1) {
            a_(r, r+1) += alpha*D;
            a_(r, r) += 3.*beta*D;
            a_(r, r-1) += -(alpha+4.*beta)*D;
            b_(r) += -fermiDirac(-kmax_)*beta*D;
            // a_(r, r) += 3.*C;
            // a_(r, r-1) += -4.*C;
            // b_(r) += -fermiDirac(-kmax_)*C;
        }
        else {
            a_(r, r+1) += alpha*D;
            a_(r, r) += 3.*beta*D;
            a_(r, r-1) += -(alpha+4.*beta)*D;
            a_(r, r-2) += beta*D;
        }
    }
}


/// Quantum correction term
void WignerSolver::quantumCorrTerm(size_t i, size_t j, double dt) {
    // Fills Boltzmann equation matrix with drift terms
    size_t r = i*nk_ + j;
    double C = d3u_(i)/dk_/dk_/dk_/24.;
    if (dt > 0) C *= dt/2.;
    // CDS1
    if (j==nk_-1)
    {
        a_(r, r-1) += C;
        a_(r, r-2) += -C/2.;
    }
    else if (j==nk_-2)
    {
        a_(r, r+1) += -C;
        a_(r, r-1) += C;
        a_(r, r-2) += -C/2.;
    }
    else if (j==0)
    {
        a_(r, r+2) += C/2.;
        a_(r, r+1) += -C;
    }
    else if (j==1)
    {
        a_(r, r+2) += C/2.;
        a_(r, r+1) += -C;
        a_(r, r-1) += C;
    }
    else
    {
        a_(r, r+2) += C/2.;
        a_(r, r+1) += -C;
        a_(r, r-1) += C;
        a_(r, r-2) += -C/2.;
    }
}


/// Solves 1D Schrodinger equation
void WignerSolver::solveSchrEq() {
    u_ = uB_ + uC_;
    size_t n = nx_;
    double m = m_, l = l_;
    double a = -1./2./m/dx_/dx_;
    arma::mat A(n, n);
    for (size_t i=n; i--;) {
        A(i,i) = -2*a+u_(i);
        if (i>0) A(i,i-1) = a;
        if (i<n-1) A(i,i+1) = a;
    }
    // VN BC
    // A(0, 0) = -a+u_(0), A(n-1, n-1) = -a+u_(n-1);
    A.brief_print("A:");
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, A);
    // for (size_t i=n; i--;)
    //     eigval.col(i) /= arma::max(eigval.col(i));
    eigvec.each_col( [](arma::vec& c){ c /= arma::sum(c); } );
    // arma::vec ne = arma::abs(eigvec.col(0))*arma::abs(eigvec.col(0));
    //
    // Results
    (eigval*AU::eV).brief_print("Eigenvalues [eV:");
    // std::cout<<eigvec<<std::endl;
    // eigvec(0).print("Eigenvector(0):");
    eigvec.save( arma::csv_name("eigvec.csv") );
    std::cout<<"Normalization: "<<arma::sum(eigvec.col(0))<<std::endl;
    arma::mat out_data;
    arma::field<std::string> header(3);
    out_data.insert_cols(0, x_), header(0) = "x [au]"; // col. 1
    out_data.insert_cols(1, u_*AU::eV), header(1) = "U [eV]";  // col. 2
    out_data.insert_cols(2, eigvec.col(0)), header(2) = "Psi";  // col. 3
    out_data.save( arma::csv_name("SchrEqResults.csv", header) );
    // std::cout<<eigvec(0)*eigval(0)<<A*eigvec(0)<<std::endl;
    //
    // Theory
    // double l = l_, m = m_;
    auto an_eigval = [l, m](int i) { return i*i*M_PI*M_PI/2./m/l/l; };
    std::cout<<"Eigen values, theory: "<<an_eigval(1)*AU::eV<<' '<<an_eigval(2)*AU::eV<<' '<<an_eigval(3)*AU::eV<<" eV"<<std::endl;
}
