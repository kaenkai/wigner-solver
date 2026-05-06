#ifndef WIGNERSOLVER_HPP
#define WIGNERSOLVER_HPP

#include "lib.hpp"
#include <armadillo>

using namespace AtomicUnits;


/**
 * Wigener function class
 */
class WignerSolver{

    size_t nx_;         // number of x-space grid points
    double lD_, lC_;    // devie len., contacts len.
    double l_;          // total length
    double dx_;         // x-space step size (lattice constant)
    size_t nk_;         // number of k-space grid points
    double kmax_;       // k-space range
    double dk_;         // k-space step size (Brillouin zone / Nk)
    size_t nk2_, nxk_;  // x/k-space nr of steps

    double dt_ = 0;  // time step size (1 fs)

    double m_ = 1;                                          // effective mass in the device
    double temp_ = 300;                                     // contacts temperature [K]
    double uR_ = 1, uL_ = 1;                                // Fermi energy in right/left contact
    double cD_ = 1;                                         // dopant concentration in contacts [AU]
    double epsilonR_ = 1;                                   // relative permitivitty (for GaAs)

    double uBias_ = 0;                                      // bias voltage [eV]
    double rR_ = 0, rM_ = 0, rG_ = 0, Gamma_ = 0, rF_ = 0;  // Scattering rate (1/tau)
    double lambda_ = 0;                                     // Localization rate
    int bcType_ = 1;                                        // boundary condition type, 1 -> Supply function
    bool useQC_;                                            // whether to use quantum correction term (third derivative of potential)

    arma::mat f_;       // Wigner function
    arma::mat fEq_;     // Equilibrium Wigner function
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
    arma::vec cdX_;     // Carrier density in x / k
    arma::vec cdK_;
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
        fEq_(arma::mat(nx_, nk_)),
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
        cdX_(arma::vec(nx_, arma::fill::zeros)),
        cdK_(arma::vec(nk_, arma::fill::zeros)),
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
        fEq_(arma::mat(nx_, nk_)),
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
        cdX_(arma::vec(nx_, arma::fill::zeros)),
        cdK_(arma::vec(nk_, arma::fill::zeros)),
        nD_(arma::vec(nx_, arma::fill::zeros)),
        currD_(arma::vec(nx_, arma::fill::zeros)),
        a_(arma::sp_mat(nxk_, nxk_)),
        b_(arma::vec(nxk_, arma::fill::zeros))
        {
        cout<<"## Start: WignerSolver constructor"<<endl;
        cout<<"# Setting up configuration space array values"<<endl;
        for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;
        cout<<"# Setting up wave vector space array values"<<endl;
        for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);
        // ----------------
        // NLP sinus values
        // ----------------
        // #pragma omp parallel for collapse(3)
        cout<<"# Setting up NLP sine values"<<endl;
        for (size_t j=0; j<nk_; ++j)
                for (size_t g=0; g<nk_; g++)
                        for (size_t h=0; h<nk2_; h++)
                                sin_(j,g*nk2_+h) = sin(2*M_PI/nk_*h*(j-g));
        cout<<"## End: WignerSolver constructor"<<endl;
    }  // End of constructor

    ~WignerSolver(){}

    size_t get_nx() { return nx_; }
    size_t get_nk() { return nk_; }
    double get_dk() { return dk_; }
    double get_dx() { return dx_; }
    double get_l() { return l_; }
    double get_lD() { return lD_; }
    double get_lC() { return lC_; }
    double get_m() { return m_; }
    double get_temp() { return temp_; }
    double get_epsilonR() { return epsilonR_; }
    double get_cD() { return cD_; }
    double get_uL() { return uL_; }
    double get_uR() { return uR_; }
    double get_dt(){ return dt_; }
    double get_rR() { return rR_; }
    double get_rM() { return rM_; }
    double get_rF() { return rF_; }
    double get_rG() { return rG_; }
    double get_lambda() { return lambda_; }
    arma::vec get_x_arr() { return x_; }
    arma::vec get_k_arr() { return k_; }
    arma::vec get_u() { return u_; }
    arma::vec get_uB() { return uB_; }
    arma::vec get_uC() { return uC_; }
    arma::vec get_du() { return du_; }
    arma::vec get_d3u() { return d3u_; }
    arma::vec get_nD() { return nD_; }
    arma::vec get_bc() { return bc_; }
    arma::mat get_wf() { return f_; }
    arma::vec get_currD() { return currD_; }

    void set_m(double m) { m_ = m; }
    void set_temp(double temp) { temp_ = temp; }
    void set_cD(double cD) { cD_ = cD; }

    void set_dt(double dt) {
        dt_ = dt;
    }

    void set_rR(double rR) { rR_ = rR; }
    void set_rM(double rM) { rM_ = rM; }
    void set_rF(double rF) { rF_ = rF; }
    void set_rG(double rG) { rG_ = rG; }
    void set_lambda(double lambda) { lambda_ = lambda; }
    void set_epsilonR(double epsilonR) { epsilonR_ = epsilonR; }
    void set_useQC(bool useQC) { useQC_ = useQC; }
    void set_bcType(int bcType) { bcType_ = bcType; }

    void set_uL(double uL) { uL_ = uL; }
    void set_uR(double uR) { uR_ = uR; }
    void set_uBias(double uBias) {
        uBias_ = uBias;
        uL_ += uBias_;
    }
    void set_uC(arma::vec uC) { 
        uC_ = uC; 
        u_ = uB_ + uC_;
    }

    void set_wf(arma::mat f) { f_ = f; }

    /** 
     * Sets up doping profile
     * @param s - smoothing parameter (0 - rectangular profile)
    */
    void set_doping_profile(double s = 0.01){
        for (size_t i=0; i<nx_; ++i)
            nD_(i) = cD_*(1+1/(1+exp((x_(i)-lC_)/s/l_))-1/(1+exp((x_(i)-l_+lC_)/s/l_)));
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

    void setBoundCond();                        // Boundary conditions
    void setEquilibriumFunction(std::string);   // Calculates equilibrium function

    void diffusionTerm(size_t, size_t, double);
    void driftTerm(size_t, size_t, double);
    void nonLocalPotentialTerm(size_t, size_t);
    void scatteringTerm(size_t, size_t, double);
    void quantumCorrTerm(size_t, size_t, double);

    double fermiDirac(double);      // Fermi-Dirac distribution
    double supplyFunction(double);  // Supply function as function of wave vector
    double sf(double, double);      // Supply function as function of energy (used for convolution with Lorentz/Gauss/Voigt profiles)

    double maxwellBoltzmann(double);
    double gaussian(double);

    double eqFun(double, double);
    double lorentz(double);         // Lorentzian profile
    double gauss(double);           // Gaussian profile
    double voigt(double);           // (pseudo-)Voigt profile

    // ----------------------------
    // Functions in WignerTools.cpp
    // ----------------------------
    void addGaussBarr(double, double, double);
    void addRectBarr(double, double, double, double);
    void addWavePacket(double, double, double, double);
    double wavePacket_TEV(double, double, double, double, double, double);
    double nC(double, double);
    double fermiInt(double, double);
    double calcFermiEn(double, double, double);
};

#endif
