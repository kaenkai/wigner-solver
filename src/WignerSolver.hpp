#ifndef WIGNERSOLVER_HPP
#define WIGNERSOLVER_HPP

#include "lib.hpp"
#include <iostream>
#include <cmath>


/**
 * Wigener function class
 */
class WignerSolver{
    size_t nx_;             // number of x-space grid points
    double l_;              // system length
    double dx_;             // x-space step size
    size_t nk_;             // number of k-space grid points
    double kmax_;           // k-space range
    double dk_;             // k-space step size
    size_t nk2_;            // nk_ / 2
    size_t nxk_;            // nx_ * nk_
    double dt_ = 0;         // time step
    double m_ = 1;          // effective mass in the device
    double temp_ = 300;     // contacts temperature [K]
    double uR_ = 1;         // Fermi energy in right contact
    double uL_ = 1;         // Fermi energy in left contact
    double epsilonR_ = 1;   // relative permittivity
    double scR_ = 0;        // dissipation
    double scM_ = 0;        // momentum randomization
    double scG_ = 0;        // contacts scattering rate
    double scF_ = 0;        // friction
    double lambda_ = 0;     // localization rate
    int bcType_ = 1;        // boundary condition type, 1 -> Supply function
    bool useQC_;            // whether to use quantum correction term (third derivative of potential)

    arma::mat f_;           // Wigner function
    arma::mat feq_;         // Equilibrium Wigner function
    arma::vec u_;           // Potential energy
    arma::vec uC_;          // Hartree potential
    arma::vec uB_;          // Conduction band offset
    arma::vec du_;          // Potential derivative
    arma::vec d3u_;         // Potential third derivative
    arma::vec bc_;          // Boundary condition
    arma::vec x_;           // Position values
    arma::vec k_;           // Wave vector values
    arma::mat sin_;         // Sine function values
    arma::vec nD_;          // Doping profile

    arma::sp_mat a_;        // Coefficient matrix for linear system
    arma::vec b_;           // Right-hand side vector for linear system

public:

    // Default constructor
    WignerSolver() :
        nx_ (100),
        l_ (100./AU::nm),
        dx_ (l_/float(nx_-1)),
        nk_ (100),
        kmax_ (M_PI/2./dx_),
        dk_ (2.*kmax_/float(nk_)),
        nk2_ (size_t(nk_/2.)),
        nxk_ (nx_*nk_),
        f_(arma::mat(nx_, nk_)),
        feq_(arma::mat(nx_, nk_)),
        u_(arma::vec(nx_, arma::fill::zeros)),
        uC_(arma::vec(nx_, arma::fill::zeros)),
        uB_(arma::vec(nx_, arma::fill::zeros)),
        du_(arma::vec(nx_, arma::fill::zeros)),
        d3u_(arma::vec(nx_, arma::fill::zeros)),
        bc_(arma::vec(nk_, arma::fill::zeros)),
        x_(arma::vec(nx_, arma::fill::zeros)),
        k_(arma::vec(nk_, arma::fill::zeros)),
        sin_(arma::mat(nk_,nk_*nk2_)),
        nD_(arma::vec(nx_, arma::fill::zeros)),
        a_(arma::sp_mat(nxk_, nxk_)),
        b_(arma::vec(nxk_, arma::fill::zeros))
        {
        std::cout<<"## Start: WignerSolver default constructor"<<std::endl;
        // ########## Configuration space array values ##########
        std::cout<<"# Setting up configuration space array values"<<std::endl;
        for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;
        // ########## Wave vector space array values ##########
        std::cout<<"# Setting up wave vector space array values"<<std::endl;
        for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);
        // ########## NLP sinus values ##########
        // #pragma omp parallel for collapse(3)
        std::cout<<"# Setting up NLP sine values"<<std::endl;
        for (size_t j=0; j<nk_; ++j)
                for (size_t g=0; g<nk_; g++)
                        for (size_t h=0; h<nk2_; h++)
                                sin_(j,g*nk2_+h) = sin(2*M_PI/nk_*h*(j-g));
        std::cout<<"## End: WignerSolver default constructor"<<std::endl;
    }  // End of constructor

    WignerSolver(size_t nx, double l, size_t nk, double kmax) :
        nx_ (nx),
        l_ (l),
        dx_ (l_/double(nx_-1)),
        nk_ (nk),
        kmax_ (kmax > 0 ? kmax : M_PI/2./dx_),
        dk_ (2.*kmax_/double(nk_-1)),
        nk2_ (size_t(nk_/2.)),
        nxk_ (nx_*nk_),
        f_(arma::mat(nx_, nk_)),
        feq_(arma::mat(nx_, nk_)),
        u_(arma::vec(nx_, arma::fill::zeros)),
        uC_(arma::vec(nx_, arma::fill::zeros)),
        uB_(arma::vec(nx_, arma::fill::zeros)),
        du_(arma::vec(nx_, arma::fill::zeros)),
        d3u_(arma::vec(nx_, arma::fill::zeros)),
        bc_(arma::vec(nk_, arma::fill::zeros)),
        x_(arma::vec(nx_, arma::fill::zeros)),
        k_(arma::vec(nk_, arma::fill::zeros)),
        sin_(arma::mat(nk_,nk_*nk2_)),
        nD_(arma::vec(nx_, arma::fill::zeros)),
        a_(arma::sp_mat(nxk_, nxk_)),
        b_(arma::vec(nxk_, arma::fill::zeros))
        {
        std::cout << "# -------------------------------" << std::endl;
        std::cout << "# Start: WignerSolver constructor" << std::endl;
        std::cout << "# Setting up configuration space array values" << std::endl;
        for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;
        std::cout << "# Setting up wave vector space array values" << std::endl;
        for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);
        // ----------------
        // NLP sinus values
        // ----------------
        // #pragma omp parallel for collapse(3)
        std::cout << "# Setting up NLP sine values" << std::endl;
        for (size_t j=0; j<nk_; ++j)
                for (size_t g=0; g<nk_; g++)
                        for (size_t h=0; h<nk2_; h++)
                                sin_(j,g*nk2_+h) = sin(2*M_PI/nk_*h*(j-g));
        std::cout << "# End: WignerSolver constructor" << std::endl;
        std::cout << "# -------------------------------\n" << std::endl;
    }  // End of constructor

    ~WignerSolver(){}

    // Getters

    size_t get_nx() { return this -> nx_; }
    size_t get_nk() { return this -> nk_; }
    double get_dk() { return this -> dk_; }
    double get_dx() { return this -> dx_; }
    double get_l() { return this -> l_; }
    double get_m() { return this -> m_; }
    double get_temp() { return this -> temp_; }
    double get_epsilonR() { return this -> epsilonR_; }
    double get_uL() { return this -> uL_; }
    double get_uR() { return this -> uR_; }
    double get_dt(){ return this -> dt_; }
    double get_scR() { return this -> scR_; }
    double get_scM() { return this -> scM_; }
    double get_scF() { return this -> scF_; }
    double get_scG() { return this -> scG_; }
    double get_lambda() { return this -> lambda_; }
    arma::vec get_x() { return this -> x_; }
    arma::vec get_k() { return this -> k_; }
    arma::vec get_u() { return this -> u_; }
    arma::vec get_uB() { return this -> uB_; }
    arma::vec get_uC() { return this -> uC_; }
    arma::vec get_du() { return this -> du_; }
    arma::vec get_d3u() { return this -> d3u_; }
    arma::vec get_nD() { return this -> nD_; }
    arma::vec get_bc() { return this -> bc_; }
    arma::mat get_f() { return this -> f_; }

    // Setters

    void set_m(double m) { this -> m_ = m; }
    void set_temp(double temp) { this -> temp_ = temp; }
    void set_epsilonR(double epsilonR) { this -> epsilonR_ = epsilonR; }
    void set_dt(double dt) { this -> dt_ = dt; }
    void set_scR(double scR) { this -> scR_ = scR; }
    void set_scM(double scM) { this -> scM_ = scM; }
    void set_scF(double scF) { this -> scF_ = scF; }
    void set_scG(double scG) { this -> scG_ = scG; }
    void set_lambda(double lambda) { this -> lambda_ = lambda; }
    void set_useQC(bool useQC) { this -> useQC_ = useQC; }
    void set_uL(double uL) { this -> uL_ = uL; }
    void set_uR(double uR) { this -> uR_ = uR; }
    void set_uC(arma::vec uC) { this -> uC_ = uC; }
    void set_f(arma::mat f) { this -> f_ = f; }

    /** 
     * Sets step doping profile
     * @param nD doping concentration in contacts (in au)
     * @param lC contact length, symmetric at each side
     * @param s smoothing parameter (0 - rectangular profile, -> 1 - smoother profile)
     * @todo move this function to WignerTools.cpp, rename, create more complex doping profiles
     */
    void setDopingProfile(double nD, double lC, double s = 0.01){
        for (size_t i=0; i<nx_; ++i)
            nD_(i) = nD*(1+1/(1+exp((x_(i)-lC)/s/l_))-1/(1+exp((x_(i)-l_+lC)/s/l_)));
    }

    // -------------------------------------
    // Solvers for Wigner/Boltzmann equation
    // -------------------------------------

    void solveBTE();
    void solveWTE();
    void solveTimeDependentBTE();
    void solveSchrEq();
    void diffusionTerm(size_t, size_t, double);
    void driftTerm(size_t, size_t, double);
    void nonLocalPotentialTerm(size_t, size_t);
    void scatteringTerm(size_t, size_t, double);
    void quantumCorrTerm(size_t, size_t, double);

    // -------------------
    // Boundary conditions
    // -------------------
    
    void setBoundCond(int);                     // Boundary conditions
    double fermiDirac(double);                  // Fermi-Dirac distribution
    double supplyFunction(double);              // Supply function as function of wave vector
    double sf(double, double);                  // Supply function as function of energy (used for convolution with Lorentz/Gauss/Voigt profiles)
    double maxwellBoltzmann(double, double);
    double gaussian(double);
    double eqFun(double, double);
    double lorentz(double);                     // Lorentzian profile
    double gauss(double);                       // Gaussian profile
    double voigt(double);                       // (pseudo-)Voigt profile

    /** 
    * Sets up equilibrium function
    * @param feq distribution function to be set as an equilibrium function
    */
    void setEquilibriumFunction(arma::mat feq){ this -> feq_ = feq; }

    void normalization() {
        double f_norm = this -> calcNorm();
        std::cout << "# distribution function normalization, norm: " << f_norm << std::endl;
        this -> f_ /= f_norm;
    }

    // ----------------------------
    // Functions in WignerTools.cpp
    // ----------------------------
    double calcNorm();
    double calcEX();
    double calcEK();
    double calcEX2();
    double calcEK2();
    double calcSDX();
    double calcSDK();
    arma::vec calcCD_X();           // Carrier density in x-space
    arma::vec calcCD_K();           // Carrier density in k-space
    arma::vec calcCurrentDensity(); // Current density according to W. R. Frensley, Phys. Rev. B 36, 1570 (1987)
    void addGaussBarr(double, double, double);
    void addRectBarr(double, double, double, double);
    void addWavePacket(double, double, double, double);
    double nC(double, double);
    double fermiInt(double, double);

    // ------------
    // IO functions
    // ------------
    void printParam();
    void saveDistFun();
};

#endif
