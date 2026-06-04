#include "lib.hpp"
#include "WignerSolver.hpp"
#include <armadillo>


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
    arma::rowvec outer = arma::trapz(x_, inner % arma::square(x_), 0);
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
    arma::vec outer = arma::trapz(k_, inner % arma::square(k_).as_row(), 1);
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
void WignerSolver::addGaussBarr(double u0 = 0.3/AU::eV, double x0 = 1000, double sX = 100) {
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
void WignerSolver::addRectBarr(double u0 = 0.3/AU::eV, double x0 = 1000, double wB = 100, double p = 10) {
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
    std::cout << "# Setting up GWP with parameters:\n"
    << "# x0 = " << x0*AU::nm << " nm = " << x0 << " au, delta_x = " << delta_x*AU::nm << " nm = " << delta_x <<" au\n"
    << "# k0 = " << k0 << " au, delta_k = " << delta_k << " au\n" <<std::endl;
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
 * Supply function as a function of wave vector
 * @param k wave vector
 * @return supply function value
 */
double WignerSolver::supplyFunction(double k){
    double mu = k > 0 ? uL_: uR_;
    double m = m_;
    double c = m/M_PI*AU::KB/AU::eV*temp_, ex = -(k*k/m/2.-mu)/(AU::KB/AU::eV*temp_);     // [au]
    if (ex < 700)
        return c * log(exp(ex)+1);
    else
        return ex > 0 ? c * ex : 0;
}


/**
 * Lorentz profile
 */
double WignerSolver::lorentz(double x) {
    double g = scG_/2.;
    return g/(x*x+g*g)/M_PI;
}


/**
 * Gauss profile
 */
double WignerSolver::gauss(double x) {
    double g = scG_/2.;
    return 1/g/sqrt(2*M_PI)*exp(-x*x/2./g/g);
}


/**
 * (pseudo-)Voigt profile
 * @see https://en.wikipedia.org/wiki/Voigt_profile 
 */
double WignerSolver::voigt(double x) {
    double g = scG_/2.;
    // Calculating eta - mixing parameter in pseudo-Voigt profile
    double gammaL = 2*scG_/2., gammaG = 2*sqrt(2*log(2))*scG_/2.;
    double gamma = pow(pow(gammaG,5) +
        2.69269*pow(gammaG,4)*gammaL +
        2.42843*pow(gammaG,3)*pow(gammaL,2) +
        4.47163*pow(gammaG,2)*pow(gammaL,3) +
        0.07842*gammaG*pow(gammaL,4) + pow(gammaL,5), 0.2);
    double eta = 1.36603*(gammaL/gamma) - 0.47719*pow(gammaL/gamma,2) + 0.11116*pow(gammaL/gamma, 3);
    return eta*( g/(x*x+g*g)/M_PI ) +                       // lorentz
           (1-eta)*( 1/g/sqrt(2*M_PI)*exp(-x*x/2./g/g) );   // gauss
}


/**
 * Contact equilibrium function (Supply function)
 * used for convolution with Lorentzian and Gauss profiles
 * @param mu chemical potential
 * @param energy energy at which equilibrium function is calculated
 * @return equilibrium function value
 */
inline double WignerSolver::eqFun(double mu, double energy) {
    double m = m_;
    double c = m/M_PI*AU::KB/AU::eV*temp_, ex = -(energy-mu)/(AU::KB/AU::eV*temp_);     // [au]
    if (ex < 700)
        return c * log(exp(ex)+1);
    else
        return ex > 0 ? c * ex : 0;
}


/**
 * Convolution
 */
double WignerSolver::conv(WignerProfile f, double k) {
    double mu = k > 0 ? uL_ : uR_;
    double u = k*k/m_/2.;
    int N = 1e5, i;
    double beta = 1/(AU::KB/AU::eV*temp_);
    double h = (40/beta+mu)/float(N);  // 40 from exp(x) -> 0 in SF
    double x0, x1, xm1, xm2, xm3;
    double fb = 0;
    for (i=1; i<N-1; i++){
        x0 = i*h, x1 = (i+1)*h;
        xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
        fb += 7 * (this->*f)(x0-u)  * eqFun(mu, x0)  +
             32 * (this->*f)(xm1-u) * eqFun(mu, xm1) +
             12 * (this->*f)(xm2-u) * eqFun(mu, xm2) +
             32 * (this->*f)(xm3-u) * eqFun(mu, xm3) +
              7 * (this->*f)(x1-u)  * eqFun(mu, x1);
    }
    return fb*2*h/4./45.;
}


/**
 * Boundary conditions
 * @param bcType boundary condition type
 * @todo implement armadillo convolution
 */
void WignerSolver::setBoundaryConditions(int bcType, double scG){
    scG_ = scG;
    switch (bcType) {
        case 0: // No electron inflow, closed system
            bc_.zeros();
            return;
        case 1: // Supply function
            for (size_t j=0; j<nk_; ++j)
                bc_(j) = supplyFunction(k_(j));
            return;
        case 2: // Splot with Lorentz profile
            for (size_t j=0; j<nk_; ++j)
                bc_(j) = conv(&WignerSolver::lorentz, k_(j));
            return;
        case 3: // Splot with  Gauss profile
            for (size_t j=0; j<nk_; ++j)
                bc_(j) = conv(&WignerSolver::gauss, k_(j));
            return;
        case 4: // Splot with (pseudo-)Voigt profile
            for (size_t j=0; j<nk_; ++j)
                bc_(j) = conv(&WignerSolver::voigt, k_(j));
            return;
        default:
            throw (bcType);
    }
}


// --------------
// FERMI INTEGRAL
// --------------


/** 
 * Density of states
 * @param m electron effective mass
 * @param T temperature
 */
inline double nC(double m, double T) { return 2.*pow(m*AU::KB*T/AU::eV/2./M_PI, 3./2.); }  // [au]


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
    eta = eta*AU::KB*T - M_PI*M_PI*AU::KB*T/12./eta;
    return eta/AU::eV;
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
    double mu = k > 0 ? uL_ : uR_;
    double m = m_;
    double ex = (k*k/m/2.-mu)/(AU::KB/AU::eV*temp_);     // [au]
    return 1./(exp(ex)+1);
}


/** 
 * Maxwell-Boltzmann distribution
 * @param nE electron density
 * @param k wave vector
 * @return Maxwell-Boltzmann distribution   value
 */
double WignerSolver::maxwellBoltzmann(double nE,double k){
    return nE*pow(2*M_PI*m_*AU::KB/AU::eV*temp_, -3/2.) * 
           exp(-k*k/2/2./(AU::KB/AU::eV*temp_));
}


/** 
 * @brief Normal distribution centered at (x_min + x_max) / 2
 * @param x_min Minimum value of x
 * @param x_max Maximum value of x
 * @param sig   Standard deviation
 * @param n     Number of points
 * @param A     Amplitude
 * @return arma::vec normal distribution
 */
arma::vec normalDistribution(double A, double sig, double x_min, double x_max, size_t n){
    arma::vec gauss(n, arma::fill::zeros);
    double x = 0;
    double mu = (x_min + x_max)*0.5;
    for (size_t i = 0; i < n; ++i){
        x = x_min + i*(x_max-x_min)/(n-1);
        gauss(i) = exp(-(x-mu)*(x-mu)/2./sig/sig) * A;
    }
    return gauss;
}


/**
 * Integral with step h using trapezoidal rule
 * @param f function
 * @param h (integration step
 * @return integral
 * @deprecated arma::trapz function is used
 */
double calcInt(arma::vec f, double h){
    size_t n = f.size();
    double ig = 0;
    for (size_t i=1; i<n/2; ++i)
        ig += (f(2*i-2)+4*f(2*i-1)+f(2*i))*h/3.;
    return ig;
}


/**
 * First derivative, second order accuracy 
 * @param f function
 * @param h differentiation step
 * @return first derivative
 */
 arma::vec calcFirstDer(arma::vec f, double h){
    size_t n = f.size();
    arma::vec df(n, arma::fill::zeros);
    for (size_t i=1; i<n-1; ++i)
        df(i) = (-f(i-1)+f(i+1))/2./h;
    df(0) = (-f(0)+f(1))/h;
    df(n-1) = (f(n-1)-f(n-2))/h;
    return df;
}


/**
 * Second derivative, second order accuracy 
 * @param f function
 * @param h differentiation step
 * @return second derivative
 */
arma::vec calcSecondDer(arma::vec f, double h){
    size_t n = f.size();
    arma::vec df(n, arma::fill::zeros);
    for (size_t i=1; i<n-1; ++i)
        df(i) = (f(i-1)-2.*f(i)+f(i+1))/h/h;
    df(0) = (2*f(0)-5.*f(1)+4*f(2)-f(3))/h/h;
    df(n-1) = (2*f(n-1)-5.*f(n-2)+4*f(n-3)-f(n-4))/h/h;
    return df;
}


/**
 * Third derivative, second order accuracy
 * @param f function
 * @param h differentiation step
 * @return third derivative
 */
arma::vec calcThirdDer(arma::vec f, double h){
    size_t n = f.size();
    arma::vec df(n, arma::fill::zeros);
    for (size_t i=2; i<n-2; ++i)
        df(i) = (-f(i-2)+2*f(i-1)-2*f(i+1)+f(i+2))/2./h/h/h;
    df(0) = (-f(0)+3.*f(1)-3*f(2)+f(3))/h/h/h;
    df(1) = (-f(1)+3.*f(2)-3*f(3)+f(4))/h/h/h;
    df(n-1) = (f(n-1)-3.*f(n-2)+3.*f(n-3)-f(n-4))/h/h/h;
    df(n-2) = (f(n-2)-3.*f(n-3)+3.*f(n-4)-f(n-5))/h/h/h;
    return df;
}


/**
 * Moving average
 */
arma::vec movAverage(arma::vec x) {
    size_t n = x.size();
    arma::vec x_av(n);
    for (size_t i=2; i<n-2; ++i)
        x_av(i) = (x(i-2)+x(i-1)+x(i)+x(i+1)+x(i+2))/5.;
    x_av(0) = x(0);
    x_av(1) = (x(0)+x(1)+x(2))/3.;
    x_av(n-2) = (x(n-1)+x(n-2)+x(n-3))/3.;
    x_av(n-1) = x(n-1);
    return x_av;
}

/**
 * Moving average
 */
arma::mat movAverage2D(arma::mat x) {
    size_t nr = x.n_rows, nc = x.n_cols;
    arma::mat x_av = arma::mat(x);
    for (size_t r=1; r<nr-1; ++r)
        for (size_t c=1; c<nc-1; ++c)       
            x_av(r, c) = (
                x(r+1, c-1) + x(r+1, c) + x(r+1, c+1) +
                x(r, c-1)   + x(r, c)   + x(r, c+1) +
                x(r-1, c-1) + x(r-1, c) + x(r-1, c+1)
            )/9.;
    return x_av;
}