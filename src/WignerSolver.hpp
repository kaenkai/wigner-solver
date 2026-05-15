#ifndef WIGNERSOLVER_HPP
#define WIGNERSOLVER_HPP

#include "lib.hpp"

using namespace AtomicUnits;


/**
 * Wigener function class
 */
class WignerSolver{
    size_t nx_;             // number of x-space grid points
    double lD_;             // device size
    double lC_;             // contacts size
    double l_;              // total length
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
    double uBias_ = 0;      // bias voltage [eV]
    double rR_ = 0;         // scattering rate (1/tau)
    double rM_ = 0;         // scattering rate (1/tau)
    double rG_ = 0;         // scattering rate (1/tau)
    double Gamma_ = 0;      // scattering rate (1/tau)
    double rF_ = 0;         // scattering rate (1/tau)
    double lambda_ = 0;     // localization rate
    int bcType_ = 1;        // boundary condition type, 1 -> Supply function
    bool useQC_;            // whether to use quantum correction term (third derivative of potential)

    arma::mat f_;       // Wigner function
    arma::mat feq_;     // Equilibrium Wigner function
    arma::mat f0_;      // Wigner function before time evolution
    arma::mat fL_;      // Wigner function for el. from LEFT contact
    arma::mat fR_;      // Wigner function for el. from RIGHT contact
    arma::vec u_;       // Potential energy
    arma::vec uC_;      // Hartree potential
    arma::vec uB_;      // Conduction band offset
    arma::vec du_;      // Potential derivative
    arma::vec d3u_;     // Potential third derivative
    arma::vec bc_;      // Boundary condition
    arma::vec x_;       // Position values
    arma::vec k_;       // Wave vector values
    arma::mat sin_;     // Sine function values
    arma::vec nD_;      // Doping profile
    arma::vec currD_;   // Current density

    arma::sp_mat a_;
    arma::vec b_;

    arma::vec iv_i_, iv_v_, iv_iRange_, iv_n_;

public:

    // Default constructor
    WignerSolver() :
        nx_ (100),
        lD_ (60./AU_nm),
        lC_ (20./AU_nm),
        l_ (lD_ + 2*lC_),
        dx_ (l_/float(nx_-1)),
        nk_ (100),
        kmax_ (M_PI/2./dx_),
        dk_ (2.*kmax_/float(nk_)),
        nk2_ (size_t(nk_/2.)),
        nxk_ (nx_*nk_),
        f_(arma::mat(nx_, nk_)),
        feq_(arma::mat(nx_, nk_)),
        f0_(arma::mat(nx_, nk_)),
        fL_(arma::mat(nx_, nk_)),
        fR_(arma::mat(nx_, nk_)),
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
        currD_(arma::vec(nx_, arma::fill::zeros)),
        a_(arma::sp_mat(nxk_, nxk_)),
        b_(arma::vec(nxk_, arma::fill::zeros))
        {
        cout<<"## Start: WignerSolver default constructor"<<endl;
        // ########## Configuration space array values ##########
        cout<<"# Setting up configuration space array values"<<endl;
        for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;
        // ########## Wave vector space array values ##########
        cout<<"# Setting up wave vector space array values"<<endl;
        for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);
        // ########## NLP sinus values ##########
        // #pragma omp parallel for collapse(3)
        cout<<"# Setting up NLP sine values"<<endl;
        for (size_t j=0; j<nk_; ++j)
                for (size_t g=0; g<nk_; g++)
                        for (size_t h=0; h<nk2_; h++)
                                sin_(j,g*nk2_+h) = sin(2*M_PI/nk_*h*(j-g));
        cout<<"## End: WignerSolver default constructor"<<endl;
    }  // End of constructor

    WignerSolver(size_t i_nx, double i_lD, double i_lC, size_t i_nk, double i_kmax) :
        nx_ (i_nx),
        lD_ (i_lD),
        lC_ (i_lC),
        l_ (lD_ + 2*lC_),
        dx_ (l_/float(nx_-1)),
        nk_ (i_nk),
        kmax_ (i_kmax > 0 ? i_kmax : M_PI/2./dx_),
        dk_ (2.*kmax_/float(nk_)),
        nk2_ (size_t(nk_/2.)),
        nxk_ (nx_*nk_),
        f_(arma::mat(nx_, nk_)),
        feq_(arma::mat(nx_, nk_)),
        f0_(arma::mat(nx_, nk_)),
        fL_(arma::mat(nx_, nk_)),
        fR_(arma::mat(nx_, nk_)),
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
        currD_(arma::vec(nx_, arma::fill::zeros)),
        a_(arma::sp_mat(nxk_, nxk_)),
        b_(arma::vec(nxk_, arma::fill::zeros))
        {
        cout << "# -------------------------------" << endl;
        cout << "# Start: WignerSolver constructor" << endl;
        cout << "# Setting up configuration space array values" << endl;
        for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;
        cout << "# Setting up wave vector space array values" << endl;
        for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);
        // ----------------
        // NLP sinus values
        // ----------------
        // #pragma omp parallel for collapse(3)
        cout << "# Setting up NLP sine values" << endl;
        for (size_t j=0; j<nk_; ++j)
                for (size_t g=0; g<nk_; g++)
                        for (size_t h=0; h<nk2_; h++)
                                sin_(j,g*nk2_+h) = sin(2*M_PI/nk_*h*(j-g));
        cout << "# End: WignerSolver constructor" << endl;
        cout << "# -------------------------------\n" << endl;
    }  // End of constructor

    ~WignerSolver(){}

    size_t get_nx() { return this -> nx_; }
    size_t get_nk() { return this -> nk_; }
    double get_dk() { return this -> dk_; }
    double get_dx() { return this -> dx_; }
    double get_l() { return this -> l_; }
    double get_lD() { return this -> lD_; }
    double get_lC() { return this -> lC_; }
    double get_m() { return this -> m_; }
    double get_temp() { return this -> temp_; }
    double get_epsilonR() { return this -> epsilonR_; }
    double get_uL() { return this -> uL_; }
    double get_uR() { return this -> uR_; }
    double get_dt(){ return this -> dt_; }
    double get_rR() { return this -> rR_; }
    double get_rM() { return this -> rM_; }
    double get_rF() { return this -> rF_; }
    double get_rG() { return this -> rG_; }
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
    arma::vec get_currD() { return this -> currD_; }

    void set_m(double m) { this -> m_ = m; }
    void set_temp(double temp) { this -> temp_ = temp; }
    void set_epsilonR(double epsilonR) { this -> epsilonR_ = epsilonR; }

    void set_dt(double dt) { this -> dt_ = dt; }

    void set_rR(double rR) { this -> rR_ = rR; }
    void set_rM(double rM) { this -> rM_ = rM; }
    void set_rF(double rF) { this -> rF_ = rF; }
    void set_rG(double rG) { this -> rG_ = rG; }
    void set_lambda(double lambda) { this -> lambda_ = lambda; }
    void set_useQC(bool useQC) { this -> useQC_ = useQC; }

    void set_uL(double uL) { this -> uL_ = uL; }
    void set_uR(double uR) { this -> uR_ = uR; }
    void set_uBias(double uBias) { this -> uBias_ = uBias; }
    void set_uC(arma::vec uC) { 
        this -> uC_ = uC; 
        this -> u_ = uB_ + uC_;
    }

    void set_f(arma::mat f) { this -> f_ = f; }

    /** 
     * Sets up doping profile
     * @param nD - doping concentration in contacts [cm^-3]
     * @param s - smoothing parameter (0 - rectangular profile)
    */
    void set_doping_profile(double nD, double s = 0.01){
        for (size_t i=0; i<nx_; ++i)
            nD_(i) = nD*(1+1/(1+exp((x_(i)-lC_)/s/l_))-1/(1+exp((x_(i)-l_+lC_)/s/l_)));
    }

    double calcCurrentDensity();  // Current density
    double calcNorm();
    double calcEX();
    double calcEK();
    double calcEK2();
    double calcSDK();
    double calcSDX();
    arma::vec calcCD_X();
    arma::vec calcCD_K();

    void printParam();
    void saveDistFun();

    void solveBTE();
    void solveWTE();
    void solveTimeEv();
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
        cout << "# distribution function normalization, norm: " << f_norm << endl;
        this -> f_ /= f_norm;
    }

    // ----------------------------
    // Functions in WignerTools.cpp
    // ----------------------------
    void addGaussBarr(double, double, double);
    void addRectBarr(double, double, double, double);
    void addWavePacket(double, double, double, double);
    double wavePacket_TEV(double, double, double, double, double, double);
    double nC(double, double);
    double fermiInt(double, double);
};

#endif
