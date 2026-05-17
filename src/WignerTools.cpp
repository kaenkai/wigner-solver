#include "lib.hpp"
#include "WignerSolver.hpp"
#include <armadillo>
#include <cmath>

using namespace AtomicUnits;


/**
 * Calculates current density in accordance with Biegel PhD (p.132) and Frensley (1987)
 * current density, being a vector, is most appropriately and accurately defined at the centerpoint between position grid nodes. 
 * the discrete expression for current density depends on the form of the diffusion operator K .
 * @return current density
 * @see W. R. Frensley. Physical Review B, 36(3):1570–1580, 1987
 */ 
arma::vec WignerSolver::calcCurrentDensity() {
    arma::vec currD(nx_, arma::fill::zeros);
    int alpha = 2, beta = 1;
    for (size_t i=1; i<nx_-2; ++i) {
        for (size_t j=0; j<nk2_; ++j)  // k < 0
            currD(i) += 
                k_(j)*(
                    alpha*f_(i,j) +
                    (alpha+3*beta)*f_(i+1,j) -
                    beta*f_(i+2,j)
                );
        for (size_t j=nk2_; j<nk_; ++j)  // k > 0
            currD(i) += 
                k_(j)*(
                    alpha*f_(i+1,j) +
                    (alpha+3*beta)*f_(i,j) -
                    beta*f_(i-1,j)
                );
    }
    currD(0) = currD(1);
    currD(nx_-2) = currD(nx_-3);
    currD(nx_-1) = currD(nx_-2);
    currD *= dk_/m_/4./M_PI/(alpha+beta);
    return currD;
}


/**
 * Calculates carrier density in x-space
 * @return carrier density in x-space
 * @see https://arma.sourceforge.net/docs.html#trapz
 */
arma::vec WignerSolver::calcCD_X(){
    return arma::trapz(k_, f_, 1).as_col()/2./M_PI;
}


/**
 * Calculates carrier density in k-space
 * @return carrier density in k-space
 * @see https://arma.sourceforge.net/docs.html#trapz
 */
arma::vec WignerSolver::calcCD_K(){
    return arma::trapz(x_, f_, 0).as_col();
}


/**
 * Calculates density function norm (integral over whole space)
 * @todo resolve problems with normalization, should be equal to 1
 * @see https://arma.sourceforge.net/docs.html#trapz
 */
double WignerSolver::calcNorm(){
    arma::vec inner = arma::trapz(k_, f_, 1).as_col();
    arma::rowvec outer = arma::trapz(x_, inner);
    return arma::as_scalar(outer);
}


/**
 * Expected value of x
 * @return expected value in x-space
 * @see https://en.wikipedia.org/wiki/Expected_value
 * @see https://arma.sourceforge.net/docs.html#trapz
 */
double WignerSolver::calcEX(){
    double norm = calcNorm();
    if (norm < 1e-12) return 0;
    arma::vec inner = arma::trapz(k_, f_, 1).as_col();
    arma::rowvec outer = arma::trapz(x_, inner % x_, 0);
    return arma::as_scalar(outer)/norm;
}


/**
 * Expected value of k
 * @return expected value in k-space
 * @see https://en.wikipedia.org/wiki/Expected_value
 * @see https://arma.sourceforge.net/docs.html#trapz
 */
double WignerSolver::calcEK(){
    double norm = calcNorm();
    if (norm < 1e-12) return 0;
    arma::rowvec inner = arma::trapz(x_, f_, 0).as_row();
    arma::vec outer = arma::trapz(k_, inner % k_.as_row(), 1);
    return arma::as_scalar(outer)/norm; 
}


/**
 * Expected value of x squared
 * @return expected value in x-space
 * @see https://en.wikipedia.org/wiki/Expected_value
 * @see https://arma.sourceforge.net/docs.html#trapz
 */
double WignerSolver::calcEX2(){
    double norm = calcNorm();
    if (norm < 1e-12) return 0;
    arma::vec inner = arma::trapz(k_, f_, 1).as_col();
    arma::rowvec outer = arma::trapz(x_, inner % arma::pow(x_, 2), 0);
    return arma::as_scalar(outer)/norm;
}


/**
 * Expected value of k squared
 * @todo verify, use armadillo functions to simplify and speed up
 * @see https://en.wikipedia.org/wiki/Expected_value
 * @see https://arma.sourceforge.net/docs.html#trapz
 */
double WignerSolver::calcEK2(){
    double norm = calcNorm();
    if (norm < 1e-12) return 0;
    arma::rowvec inner = arma::trapz(x_, f_, 0).as_row();
    arma::vec outer = arma::trapz(k_, inner % arma::pow(k_, 2).as_row(), 1);
    return arma::as_scalar(outer)/norm; 
}


/**
 * Standard deviation of x
 * @see https://en.wikipedia.org/wiki/Standard_deviation
 */
double WignerSolver::calcSDX(){
    return sqrt(calcEX2()-calcEX()*calcEX());
}


/**
 * Standard deviation of k
 * @see https://en.wikipedia.org/wiki/Standard_deviation
 */
double WignerSolver::calcSDK(){
    return sqrt(calcEK2()-calcEK()*calcEK());
}


/** 
 * Adds gaussian barrier
 * @param u0 height
 * @param x0 position
 * @param sX width
 */
void WignerSolver::addGaussBarr(double u0 = 0.3/AU_eV, double x0 = 1000, double sX = 100) {
    for (size_t i=0; i<nx_; ++i)
        uB_(i) += exp(-(x_(i)-x0)*(x_(i)-x0)/sX/sX)*u0;
}


/** 
 * Adds rectangular barrier
 * @param u0 height
 * @param x0 position
 * @param wB width
 * @param p rectangularity parameter
 */
void WignerSolver::addRectBarr(double u0 = 0.3/AU_eV, double x0 = 1000, double wB = 100, double p = 10) {
    double sig = wB/2.;
    for (size_t i=0; i<nx_; ++i)
        uB_(i) += exp(-pow(x_(i)-x0, 2*p)/2./pow(sig, 2*p))*u0;
}


// --------------------
// Gaussian wave packet
// --------------------


/** 
 * Add gaussian wave packet to distribution function
 * @param x0 GWP position in x-space
 * @param delta_x GWP blur in x-space
 * @param k0 GWP position in k-space
 * @param delta_k GWP blur in k-space
 * @see BJS habilitation thesis, p. 62, eq. (2.44)
 */
void WignerSolver::addWavePacket(double x0, double delta_x, double k0, double delta_k) {
    cout << "# Setting up GWP with parameters:\n"
    << "# x0 = " << x0*AU_nm << " nm = " << x0 << " au, delta_x = " << delta_x*AU_nm << " nm = " << delta_x <<" au\n"
    << "# k0 = " << k0 << " au, delta_k = " << delta_k << " au\n" <<endl;
    for (size_t i=0; i<nx_; ++i) {
        for (size_t j=0; j<nk_; ++j)
            f_(i,j) += exp(
                - std::pow((k_(j)-k0)/delta_k, 2)/2
                - std::pow((x_(i)-x0)/delta_x, 2)/2
            ) * 2;
    }
}


// -------------------
// Boundary conditions
// -------------------


/**
 * Boundary conditions
 * @param bcType boundary condition type
 * @todo SIMPLIFY!!!
 * @todo implement armadillo convolution
 */
void WignerSolver::setBoundCond(int bcType){
    bcType_ = bcType;;
    if (bcType_ == 0)
        bc_.zeros();
    else if (bcType_ == 1)
        for (size_t j=0; j<nk_; ++j)
            bc_(j) = supplyFunction(k_(j));
    else if (bcType_ == -1)
        for (size_t j=0; j<nk_; ++j)
            bc_(j) = gaussian(k_(j));
    else {
        if (scG_ > 0) {
            if (bcType_ == 2 || bcType_ == -2)
                for (size_t j=0; j<nk_; ++j)
                    bc_(j) = lorentz(k_(j));
            else if (bcType_ == 3 || bcType_ == -3)
                for (size_t j=0; j<nk_; ++j)
                    bc_(j) = gauss(k_(j));
            else if (bcType_ == 4 || bcType_ == -4)
                for (size_t j=0; j<nk_; ++j)
                    bc_(j) = voigt(k_(j));
            else {
                cout << "# ERROR WHILE SETTING BOUNDARY CONDITIONS" << endl;
                cout << "# bcType_ = " << bcType_
                    << " IS WRONG BOUNDARY CONDITION TYPE INT" << endl;
                exit(0);
            }
        }
        else {
            cout << "# ERROR WHILE SETTING BOUNDARY CONDITIONS" << endl;
            cout << "# scG_ = " << scG_
                << " SCATTERING RATE SHOULD BE GREATER THAN 0" << endl;
            exit(0);
        }
    }
    arma::mat bc_out;
    bc_out.insert_cols(0, k_);
    bc_out.insert_cols(1, bc_);
    bc_out.save("output/bc.out", arma::raw_ascii);
}


/**
 * Supply function as a function of wave vector
 * @param k wave vector
 * @return supply function value
 */
double WignerSolver::supplyFunction(double k){
    double mu = k > 0 ? uL_ + uBias_ : uR_;
    double m = m_;
    double c = m/M_PI*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);     // [au]
    if (ex < 700)
        return c * log(exp(ex)+1);
    else
        return ex > 0 ? c * ex : 0;
}


/**
 * Gaussian as a function of wave vector
 * @param k wave vector
 * @return Gaussian value
 */
double WignerSolver::gaussian(double k) {
    return pow((2.*M_PI),-1/2.) * exp(-k*k/2.);
}


/**
 * Supply function as a function of energy
 * used for convolution with Lorentzian and Gauss profiles
 * @param mu chemical potential
 * @param energy energy at which supply function is calculated
 * @return supply function value
 */
inline double WignerSolver::sf(double mu, double energy) {
    double m = m_;
    double c = m/M_PI*KB/AU_eV*temp_, ex = -(energy-mu)/(KB/AU_eV*temp_);     // [au]
    if (ex < 700)
        return c * log(exp(ex)+1);
    else
        return ex > 0 ? c * ex : 0;
}


/**
 * Contact equilibrium function (used only in convolution)
 * for quantum boundary conditions -> supply function
 * for classical boundary conditions -> 1D Maxwell-Boltzmann distribution
 * @param mu chemical potential
 * @param energy energy at which equilibrium function is calculated
 * @return equilibrium function value
 * @todo review classical part of the function, decide on the form of the function
 * @see https://en.wikipedia.org/wiki/Maxwell-Boltzmann_distribution
 */
double WignerSolver::eqFun(double mu, double energy) {
    double kBT = KB/AU_eV*temp_;
    return bcType_ > 0 ? sf(mu, energy) : pow(2.*M_PI*kBT/m_,-1/2.) * exp(-energy/kBT);
}


/**
 * Lorentzian and supply function convolution
 * @param k wave vector at which convolution is calculated
 * @return Lorentzian and supply function convolution
 */
double WignerSolver::lorentz(double k) {
    double mu = k > 0 ? uL_ + uBias_ : uR_;
    double m = m_;
    double u = k*k/m/2., g = scG_/2.;
    double beta = 1/(KB/AU_eV*temp_);
    // ---------------
    // Lorentz profile
    // ---------------
    auto f = [g](double x) { return g/(x*x+g*g)/M_PI; };
    // -----------
    // Convolution
    // -----------
    size_t N = 1e4, i;
    double h = (40/beta+mu)/float(N);  // 40 from exp(x) -> 0 in SF
    double fb = 0;
    double x0, x1, xm1, xm2, xm3;
    for (i=1; i<N-1; i++){
        x0 = i*h, x1 = (i+1)*h;
        xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
        fb += 7*f(x0-u) * eqFun(mu, x0) +
            32*f(xm1-u) * eqFun(mu, xm1) +
            12*f(xm2-u) * eqFun(mu, xm2) +
            32*f(xm3-u) * eqFun(mu, xm3) +
            7*f(x1-u) * eqFun(mu, x1);
    }
    fb *= 2*h/4./45.;
    return fb;
}


/**
 * Gauss and supply function convolution
 * @param k wave vector at which convolution is calculated
 * @return Gauss and supply function convolution
 */
double WignerSolver::gauss(double k) {
    double mu = k > 0 ? uL_ + uBias_ : uR_;
    double m = m_;
    double u = k*k/m/2., g = scG_/2.;
    double beta = 1/(KB/AU_eV*temp_);
    // -------------
    // Gauss profile
    // -------------
    auto f = [g](double x) -> double { return 1/g/sqrt(2*M_PI)*exp(-x*x/2./g/g); };
    // -----------
    // Convolution
    // -----------
    int N = 1e4, i;
    double fb = 0;
    double h = (40/beta+mu)/float(N);  // 40 from exp(x) -> 0 in SF
    double x0, x1, xm1, xm2, xm3;
    for (i=1; i<N-1; i++){
        x0 = i*h, x1 = (i+1)*h;
        xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
        fb += 7*f(x0-u) * eqFun(mu, x0) +
            32*f(xm1-u) * eqFun(mu, xm1) +
            12*f(xm2-u) * eqFun(mu, xm2) +
            32*f(xm3-u) * eqFun(mu, xm3) +
            7*f(x1-u) * eqFun(mu, x1);
    }
    return fb*2*h/4./45.;
}


/**
 * (pseudo-)Voigt profile and supply function convolution
 * @param k wave vector at which convolution is calculated
 * @return (pseudo-)Voigt profile and supply function convolution
 * @see https://en.wikipedia.org/wiki/Voigt_profile 
 */
double WignerSolver::voigt(double k) {
    double mu = k > 0 ? uL_ + uBias_ : uR_;
    double m = m_;
    double u = k*k/m/2., g = scG_/2.;
    double beta = 1/(KB/AU_eV*temp_);
    // Calculating eta - mixing parameter in pseudo-Voigt profile
    double gammaL = 2*scG_/2., gammaG = 2*sqrt(2*log(2))*scG_/2.;
    double gamma = pow(pow(gammaG,5) +
        2.69269*pow(gammaG,4)*gammaL +
        2.42843*pow(gammaG,3)*pow(gammaL,2) +
        4.47163*pow(gammaG,2)*pow(gammaL,3) +
        0.07842*gammaG*pow(gammaL,4) + pow(gammaL,5), 0.2);
    double eta = 1.36603*(gammaL/gamma) - 0.47719*pow(gammaL/gamma,2) + 0.11116*pow(gammaL/gamma, 3);
    // ----------------------
    // (pseudo-)Voigt profile
    // ----------------------
    auto f = [eta, g](double x) -> double { 
        return eta*( g/(x*x+g*g)/M_PI ) +  // lorentz
               (1-eta)*( 1/g/sqrt(2*M_PI)*exp(-x*x/2./g/g) );  // gauss
        };
    // -----------
    // Convolution
    // -----------
    int N = 1e4, i;
    double fb = 0;
    double h = (40/beta+mu)/float(N);  // 40 from exp(x) -> 0 in SF
    double x0, x1, xm1, xm2, xm3;
    for (i=1; i<N-1; i++){
        x0 = i*h, x1 = (i+1)*h;
        xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
        fb += 7*f(x0-u) * eqFun(mu, x0) +
            32*f(xm1-u) * eqFun(mu, xm1) +
            12*f(xm2-u) * eqFun(mu, xm2) +
            32*f(xm3-u) * eqFun(mu, xm3) +
            7*f(x1-u) * eqFun(mu, x1);
    }
    return fb*2*h/4./45.;
}


// --------------
// FERMI INTEGRAL
// --------------


/** 
 * Density of states
 * @param m electron effective mass
 * @param T temperature
 */
inline double nC(double m, double T) { return 2.*pow(m*KB*T/AU_eV/2./M_PI, 3./2.); }  // [au]


/**
 * Calculates Fermi integral
 * @param n doping concentration
 * @param eta relative Fermi level: (Ec-Ef)/kb/T
 */
inline double fermiInt(double n, double eta) {
    auto fi = [n, eta](double x) { return pow(x, n)*exp(eta-x)/(exp(eta-x)+1); } ;
    size_t N = 100;
    double h = 1./float(N), j;
    double i1 = 0, i2 = 0;
    for (size_t i=1; i<N; ++i){
        j = i*h;
        i1 += (fi(j)+4*fi((j+1)*h)+fi((j+2)*h))*h/3.;
        i2 += pow( (fi(1/j/h)+4*fi(1/(j+1)/h)+fi(1/(j+2)/h))/(j*h), 2. )*h/3.;
        // i1 += (fi(j)+fi(j+h))*h/2.;
        // i2 += pow( (fi(1./j/h)+fi(1./(j+1)/h))/j/h, 2. )*h/2.;
    }
    double integral = i1+i2;
    // std::tgamma(n+1) <- calculates gamma function for n+1
    return 1./std::tgamma(n+1.)*integral;
}


/**
 * Calculates Fermi energy from electron density using Fermi integral
 * @param n0 electron density
 */
double calcFermiEn(double n0, double m, double T) {
    double eta = 0., f = 0.;
    double x = n0/nC(m, T);  // For GaAs, 300 K -> x = 4.6
    //  Looking for eta
    if (x > fermiInt(.5, 5.)) {
        eta = pow(x*std::tgamma(5./2.), 2./3.);
    }
    else if (x < fermiInt(.5, -2.))
        eta = log(x);  // *pow(2,-3/2.)*x
    else {
        if (x > fermiInt(.5, 0.)) {
            while (f < x) {
                eta += 1e-3;
                f = fermiInt(.5, eta);
            }
        }
        else {
            while (f > x) {
                eta -= 1e-3;
                f = fermiInt(.5, eta);
            }
        }
    }
    eta = eta*KB*T - M_PI*M_PI*KB*T/12./eta;
    return eta/AU_eV;
}


// --------------
// MISC
// --------------


/** 
 * Fermi dirac distribution
 * to consider: camelCase -> snake_case, or vice versa, for all functions
 * @param k wave vector
 */
double WignerSolver::fermiDirac(double k){
    double mu = k > 0 ? uL_ + uBias_ : uR_;
    double m = m_;
    double ex = (k*k/m/2.-mu)/(KB/AU_eV*temp_);     // [au]
    return 1./(exp(ex)+1);
}


/** 
 * Maxwell-Boltzmann distribution
 * @param nE electron density
 * @param k wave vector
 * @return Maxwell-Boltzmann distribution   value
 */
double WignerSolver::maxwellBoltzmann(double nE,double k){
    return nE*pow(2*M_PI*m_*KB/AU_eV*temp_, -3/2.) * 
           exp(-k*k/2/2./(KB/AU_eV*temp_));
}

