#include "lib.hpp"
#include "WignerFunction.hpp"

using namespace AtomicUnits;


/** 
 * Sets up equilibrium function
 * if input_file is not empty: reads equilibrium from input_file
 * solves WTE for 0 V bias otherwise
 * @param input_file input file, empty string as default
 */
void WignerFunction::setEquilibriumFunction(std::string input_file = ""){
    if (input_file.empty())  {
        double uBias = uBias_, rR = rR_;
        uBias_ = 0, rR_ = 0.;
        solveWignerEq();
        for (size_t i=0; i<nx_; ++i)
            for (size_t j=0; j<nk_; ++j)
                fEq_(i,j) = f_(i,j);
        f_.zeros();
        uBias_ = uBias, rR_ = rR;
    }
    else {
        fEq_.load(input_file);
        if (fEq_.size() != f_.size()) {
            cout<<"ERROR IN setEquilibriumFunction: fEq_.size() != f_.size()"<<endl;
            exit(0);
        }
    }
}


/**
 * Calculates current density
 * @todo choose one differentiation scheme, review calculation method
 */ 
double WignerFunction::calcCurr() {
    currD_.zeros();  // Current density  // array<double>
    double a = 2., b = 1.;
    //  cur1 - current density in node i-1/2
    //  cur2 - current density in node i+1/2
    // Current density is calculated in node i+1/2
    if (diffSch_J_ == "CDS1") {
        for (size_t i=1; i<nx_; ++i) {
            for (size_t j=0; j<nk2_; ++j)  // k < 0
                currD_(i) += k_(j)*f_(i, j);
                // currD_(i) += k_(j)*f_(i+1, j);  // j(i-1/2)
            for (size_t j=nk2_; j<nk_; ++j)  // k > 0
                currD_(i) += k_(j)*f_(i-1, j);
                // currD_(i) += k_(j)*f_(i, j);  // j(i-1/2)
        }
        currD_(0) = currD_(1);
        currD_(nx_-1) = currD_(nx_-2);
    }
    else if (diffSch_J_ == "UDS1") {
        for (size_t i=1; i<nx_; ++i) {
            for (size_t j=0; j<nk2_; ++j)  // k < 0
                currD_(i) += k_(j)*f_(i, j);
                // currD_(i) += k_(j)*f_(i+1, j);  // j(i-1/2)
            for (size_t j=nk2_; j<nk_; ++j)  // k > 0
                currD_(i) += k_(j)*f_(i-1, j);
                // currD_(i) += k_(j)*f_(i, j);  // j(i-1/2)
            // currD_(i) /= 2.;
            // currD_(i) -= rR_*(nc(i)-nc_eq(i));
        }
        currD_(0) = currD_(1);
        currD_(nx_-1) = currD_(nx_-2);
    }
    else if (diffSch_J_ == "UDS2") {
        for (size_t i=1; i<nx_-2; ++i) {
            for (size_t j=0; j<nk2_; ++j)  // k < 0
                currD_(i) += k_(j)*(3*f_(i+1,j)-f_(i+2,j))*dk_/2./m_;
                // currD_(i) += k_(j)*(3.0*f_(i,j)-f_(i+1,j))*dk_/2./m_;  // j(i-1/2)
            for (size_t j=nk2_; j<nk_; ++j)  // k > 0
                currD_(i) += k_(j)*(3*f_(i,j)-f_(i-1,j))*dk_/2./m_;
                // currD_(i) += k_(j)*(3.0*f_(i-1,j)-f_(i-2,j))*dk_/2./m_;  // j(i-1/2)
        }
        currD_(0) = currD_(1);
        currD_(nx_-2) = currD_(nx_-3);
        currD_(nx_-1) = currD_(nx_-2);
    }
    else if (diffSch_J_ == "UDS3") {
        for (size_t i=2; i<nx_-3; ++i) {
            for (size_t j=0; j<nk2_; ++j)  // k < 0
                currD_(i) += k_(j)*(2*f_(i+3,j)-7*f_(i+2,j)+11*f_(i+1,j))*dk_/6./m_;
            for (size_t j=nk2_; j<nk_; ++j)  // k > 0
                currD_(i) += k_(j)*(2*f_(i-2,j)-7*f_(i-1,j)+11*f_(i,j))*dk_/6./m_;
        }
        currD_(1) = currD_(2);
        currD_(0) = currD_(1);
        currD_(nx_-3) = currD_(nx_-4);
        currD_(nx_-2) = currD_(nx_-3);
        currD_(nx_-1) = currD_(nx_-2);
    }
    else if (diffSch_J_ == "HDS22") {
        for (size_t i=1; i<nx_-2; ++i) {
            for (size_t j=0; j<nk2_; ++j)  // k < 0
                currD_(i) += k_(j)*dk_/m_*(a*f_(i,j)+(a+3*b)*f_(i+1,j)-b*f_(i+2,j))/2./(a+b);
            for (size_t j=nk2_; j<nk_; ++j)  // k > 0
                currD_(i) += k_(j)*dk_/m_*(a*f_(i+1,j)+(a+3*b)*f_(i,j)-b*f_(i-1,j))/2./(a+b);
        }
        currD_(0) = currD_(1);
        currD_(nx_-2) = currD_(nx_-3);
        currD_(nx_-1) = currD_(nx_-2);
    }
    else {
        cout<<"ERROR: WRONG DIFFERENTIATION SCHEME IN CURRENT DENSITY. "
        "PLEASE CHOOSE: UDS1, UDS2, UDS3 OR HDS22"<<endl;
        exit(0);
    }
    if (bcType_ > 0) currD_ /= 2.*M_PI;
    return sum(currD_)*dx_/l_;
}


/**
 * Calculates carrier density in x space
 * @todo decide on trapezoid or simpson method
 */
arma::vec WignerFunction::calcCD_X(){
    cdX_.zeros();
    for (size_t i=nx_; i--;) {
        for (size_t j=1; j<nk_/2; ++j)
            // cdX_(i) += (f_(i,j-1) + f_(i,j))*dk_/2.;  //  / 2./M_PI  // trapezoid
            cdX_(i) += (f_(i,2*j-2)+4*f_(i,2*j-1)+f_(i,2*j))*dk_/3.;  // simpson
    }
    if (bcType_ > 0) cdX_ /= 2.*M_PI;
    return cdX_;
}


/**
 * Calculates carrier density in k space
 */
arma::vec WignerFunction::calcCD_K(){
    cdK_.zeros();
    for (size_t j=nk_; j--;) {
        for (size_t i=1; i<nx_/2; ++i)
            // cd(j) += (f_(i-1,j) + f_(i,j))*dx_/2.;  // trapezoid
            cdK_(j) += (f_(2*i-2,j)+4*f_(2*i-1,j)+f_(2*i,j))*dx_/3.;  // simpson
    }
    return cdK_;
}


/**
 * Calculates density function norm (integral over whole space)
 */
double WignerFunction::calcNorm(){
    double sum_x = 0;
    double sum_k = 0;
    for (size_t i=nx_; i--;) {
        sum_k = 0;
        for (size_t j=1; j<nk_/2; ++j)
            // sum_k += ( f_(i,j) + f_(i,j-1) ) * dk_/2.;
            sum_k += (f_(i,2*j-2)+4*f_(i,2*j-1)+f_(i,2*j))*dk_/3.;  // simpson
        // sum_k += ( f_(i,nk_-1) + f_(i,nk_-2) ) * dk_/2.;
        if (i != 0 && i != nx_-1)
            sum_x += sum_k * dx_;
        else
            sum_x += sum_k * dx_/2.;
    }
    if (bcType_ > 0) sum_x /= 2.*M_PI;
    return sum_x;  // /2./M_PI
}


/**
 * Expected value in x-space
 */
double WignerFunction::calcEX(){
    double sum_x = 0;
    double sum_k = 0;
    for (size_t i=nx_; i--;) {
        sum_k = 0;
        for (size_t j=nk_-1; j--;)
            sum_k += ( f_(i,j) + f_(i,j+1) ) * x_(i)  * dk_/2.;
        sum_k += ( f_(i,nk_-1) + f_(i,nk_-2) ) * x_(i) * dk_/2.;
        if (i != 0 && i != nx_-1)
            sum_x += sum_k * dx_;
        else
            sum_x += sum_k * dx_/2.;
    }
    if (bcType_ > 0) sum_x /= 2.*M_PI;
    return sum_x / calcNorm();
}


/**
 * Expected value in k-space
 */
double WignerFunction::calcEK(){
    double sum_x = 0;
    double sum_k = 0;
    for (size_t i=nx_; i--;) {
        sum_k = 0;
        for (size_t j=1; j<nk2_; ++j)
            // sum_k += (f_(i,j-1) + f_(i,j)) * k_(j) * dk_/2.;  //  / 2./M_PI  // trapezoid
            sum_k += (f_(i,2*j-2)+4*f_(i,2*j-1)+f_(i,2*j)) * k_(j) * dk_/3.;  // simpson
        // sum_k += ( f_(0) + f_(1) ) * k_(0) * dk_/2.;
        if (i != 0 && i != nx_-1)
            sum_x += sum_k * dx_;
        else
            sum_x += sum_k * dx_/2.;
    }
    if (bcType_ > 0) sum_x /= 2.*M_PI;
    return sum_x / calcNorm();
}


/**
 * Expected value in k-space, squared
 */
double WignerFunction::calcEK2(){
    double sum_x = 0;
    double sum_k = 0;
    for (size_t i=nx_; i--;) {
        sum_k = 0;
        for (size_t j=nk_-1; j--;)
            sum_k += ( f_(i,j) + f_(i,j+1) ) * k_(j)*k_(j) * dk_/2.;
        sum_k += ( f_(i,nk_-1) + f_(i,nk_-2) ) * k_(nk_-1) * k_(nk_-1) * dk_/2.;
        if (i != 0 && i != nx_-1)
            sum_x += sum_k * dx_;
        else
            sum_x += sum_k * dx_/2.;
    }
    if (bcType_ > 0) sum_x /= 2.*M_PI;
    return sum_x / calcNorm();
}


/**
 * Standard deviation in k-space
 */
double WignerFunction::calcSDK(){
    double ev = 0, ev2 = 0, s_ev, s_ev2;
    for (size_t i=nx_; i--;) {
        s_ev = 0, s_ev2 = 0;
        for (size_t j=nk_-1; j--;) {
            s_ev += ( f_(i,j) + f_(i,j+1) ) * k_(j) * dk_/2.;
            s_ev2 += ( f_(i,j) + f_(i,j+1) ) * k_(j) * k_(j) * dk_/2.;
        }
        if (i == 0 || i == nx_-1) {
            ev += s_ev * dx_/2.;
            ev2 += s_ev2 * dx_/2.;
        }
        else {
            ev += s_ev * dx_;
            ev2 += s_ev2 * dx_;
        }
    }
    if (bcType_ > 0) ev /= 2.*M_PI;
    if (bcType_ > 0) ev2 /= 2.*M_PI;
    ev = ev / calcNorm();
    ev2 = ev2 / calcNorm();
    return sqrt(ev2-ev*ev);
}


/**
 * Calculates standard deviation in x
 */
double WignerFunction::calcSDX(){
    double ev = 0, ev2 = 0, s_ev, s_ev2;
    for (size_t i=0; i<nx_; ++i) {
        s_ev = 0, s_ev2 = 0;
        for (size_t j=1; j<nk_; ++j) {
            s_ev += ( f_(i,j) + f_(i,j-1) ) * x_(i) * dk_/2.;
            s_ev2 += ( f_(i,j) + f_(i,j-1) ) * x_(i) * x_(i) * dk_/2.;
        }
        if (i == 0 || i == nx_-1) {
            ev += s_ev * dx_/2.;
            ev2 += s_ev2 * dx_/2.;
        }
        else {
            ev += s_ev * dx_;
            ev2 += s_ev2 * dx_;
        }
    }
    if (bcType_ > 0) ev /= 2.*M_PI;
    if (bcType_ > 0) ev2 /= 2.*M_PI;
    ev = ev / calcNorm();
    ev2 = ev2 / calcNorm();
    return sqrt(ev2-ev*ev);
}


/** 
 * Sets system potential as linear drop from uBias/2 to -uBias/2
 * @param uBias - potential bias
 */ 
void WignerFunction::setPotBias(double uBias) {
    uBias_ = uBias;
    uC_.zeros();
    double x;
    for (size_t i = 0; i < nx_; ++i) {
        x = x_(i);
        uC_(i) = uBias_*(0.5-x/l_);
    }
}

/** 
 * Adds gaussian barrier
 * @param u0 height
 * @param x0 position
 * @param sX width
 */
void WignerFunction::addGaussBarr(double u0 = 0.3/AU_eV, double x0 = 1000, double sX = 100) {
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
void WignerFunction::addRectBarr(double u0 = 0.3/AU_eV, double x0 = 1000, double wB = 100, double p = 10) {
    double sig = wB/2.;
    for (size_t i=0; i<nx_; ++i)
        uB_(i) += exp(-pow(x_(i)-x0, 2*p)/2./pow(sig, 2*p))*u0;
}


/** 
 * Gaussian wave packet initial conditions
 * @param gwp_x0, gwp_k0 GWP position
 * @param gwp_dx, gwp_dk GWP size
 */
void WignerFunction::addWavePacket(double gwp_x0, double gwp_dx,
    double gwp_k0, double gwp_dk) {
    double A = cD_*lC_ / (gwp_dx*gwp_dk*2*M_PI);
    cout<<"# Setting up GWP with parameters:\n"
    <<"# gwp_x0 = "<<gwp_x0*AU_nm<<" nm, "<<gwp_x0<<" a.u.\n"
    <<"# gwp_dx = "<<gwp_dx*AU_nm<<" nm, "<<gwp_dx<<" a.u.\n"
    <<"# gwp_p0 = "<<gwp_k0<<" a.u.\n"
    <<"# gwp_dp = "<<gwp_dk<<" a.u.\n"
    <<"# gwp_A = "<<A/AU_cm2<<" cm^-2, "<<A<<" a.u.\n"
    <<"# cD_*lC_ = "<<cD_*lC_/AU_cm2<<" cm^-2\n"<<endl;
    double sx = 2*gwp_dx*gwp_dx, sk = 2*gwp_dk*gwp_dk;
    for (size_t i=0; i<nx_; ++i) {
        for (size_t j=0; j<nk_; ++j)
            f_(i,j) += exp(
                - (k_(j)-gwp_k0)*(k_(j)-gwp_k0)/sk
                - (x_(i)-gwp_x0)*(x_(i)-gwp_x0)/sx  ) * A; // * gwp_A_, / M_PI
    }
}


/** 
 * Analitical solution of wave packet time evolution (no potential)
 */
double WignerFunction::wavePacket_TEV(double gwp_x0, double gwp_dx,
    double gwp_k0, double gwp_dk, double x, double k)
    { return exp( -(k-gwp_k0)*(k-gwp_k0)/2./gwp_dk/gwp_dk
        -(x-gwp_x0)*(x-gwp_x0)/2./gwp_dx/gwp_dx ) / M_PI; }


/**
 * Boundary conditions
 * @todo implement armadillo convolution
 */
void WignerFunction::setBoundCond(){
    uL_ = uBias_BC_ ? uL_ + uBias_ : uL_;
    Gamma_ = rG_*.5;
    if (bcType_ == 1)
        for (size_t j=0; j<nk_; ++j)
            bc_(j) = supplyFunction(k_(j));
    else if (bcType_ == 0)
        bc_.zeros();
    else if (bcType_ == -1)
        for (size_t j=0; j<nk_; ++j)
            bc_(j) = gaussian_bc(k_(j));
    else {
        if (rG_ > 0) {
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
            cout << "# rG_ = " << rG_
                << " SCATTERING RATE SHOULD BE GREATER THAN 0" << endl;
            exit(0);
        }
    }
    std::ofstream file;
    file.open("OutData/BC.out", std::ios::out);
    for (size_t j=0; j<nk_; ++j) {
        file<<k_(j)<<' '<<bc_(j)<<'\n';
    }
    file<<"# "<<calcInt(bc_, dk_)/2./M_PI/AU_cm3;
    file.close();
}


/// Supply function (E(k))
double WignerFunction::supplyFunction(double k){
    double mu = k > 0 ? uL_ : uR_;
    double m = m_;
    double c = m/M_PI*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);     // [au]
    if (ex < 700)
        return c * log(exp(ex)+1);
    else{
        if (ex > 0)
            return c * ex;
        else return 0;
    }
}

/// Fermi dirac distribution
double WignerFunction::fermiDirac(double k){
    double mu = k > 0 ? uL_ : uR_;
    double m = m_;
    double ex = (k*k/m/2.-mu)/(KB/AU_eV*temp_);     // [au]
    return 1./(exp(ex)+1);
}


/// Maxwell-Boltzmann distribution
double WignerFunction::maxwell_boltzmann(double k){;
    double m = m_;
    double c = cD_*pow(2*M_PI*m*KB/AU_eV*temp_, -3/2.);
    double ex = exp(-k*k/2/2./(KB/AU_eV*temp_));
    return c * ex;
}


/// Gaussian with sigma parameter
double WignerFunction::gaussian_bc(double k){
    double m = m_;
    double kBT = KB/AU_eV*temp_;
    double sigma = 1;
    double c = pow((2.*M_PI*m*kBT)*sigma*sigma,-1/2.)*cD_;
    double ex = exp(-(k*k/2./m)/kBT/sigma/sigma);  // (k-pF)*(k-pF)
    return c * ex;
}


/// Supply function (x)
inline double WignerFunction::sf_x(double mu, double x) {
    double m = m_;
    double c = m/M_PI*KB/AU_eV*temp_, ex = -(x-mu)/(KB/AU_eV*temp_);     // [au]
    if (ex < 700)
        return c * log(exp(ex)+1);
    else{
        if (ex > 0)
            return c * ex;
        else return 0;
    }
}


/// Equilibrium function for electrons in contact
double WignerFunction::eqFun_x(double mu, double x) {
    double m = m_;
    double kBT = KB/AU_eV*temp_;
    double c, ex;
    if (bcType_ > 0) {
        c = m/M_PI*KB/AU_eV*temp_, ex = -(x-mu)/(KB/AU_eV*temp_);     // [au]
        if (ex < 700)
            return c * log(exp(ex)+1);
        else{
            if (ex > 0)
                return c * ex;
            else return 0;
        }
    }
    else {
        double sigma = 1;
        c = pow((2.*M_PI*m*kBT)*sigma*sigma,-1/2.)*cD_;
        ex = exp(-x/kBT/sigma/sigma);     // [au]
        return c * ex;
    }
}


/// Lorentz * SF convolution
double WignerFunction::lorentz(double k){
    double mu = k > 0 ? uL_ : uR_;
    double m = m_;
    double u = k*k/m/2., g = Gamma_;
    double beta = 1/(KB/AU_eV*temp_);
    // Lorentz profile
    auto f = [g](double x) { return g/(x*x+g*g)/M_PI; };
    // convolution
    size_t N = 1e4, i;
    double h = (40/beta+mu)/float(N);  // 40 from exp(x) -> 0 in SF
    double fb = 0;
    double x0, x1, xm1, xm2, xm3;
    for (i=1; i<N-1; i++){
        x0 = i*h, x1 = (i+1)*h;
        xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
        fb += 7*f(x0-u) * eqFun_x(mu, x0) +
            32*f(xm1-u) * eqFun_x(mu, xm1) +
            12*f(xm2-u) * eqFun_x(mu, xm2) +
            32*f(xm3-u) * eqFun_x(mu, xm3) +
            7*f(x1-u) * eqFun_x(mu, x1);
    }
    fb *= 2*h/4./45.;
    //  a(N, arma::fill::zeros), b(N, arma::fill::zeros), c(N, arma::fill::zeros);
    // for (i=0; i<N; ++i) {
    //     b(i) = f(i*h-u), a(i) = sf_x(mu, i*h);
    // }
    // c = conv(a, b, "same");
    // fb = sum(c);
    return fb;
}


/// Gauss * SF convolution
double WignerFunction::gauss(double k){
    double mu = k > 0 ? uL_ : uR_;
    double m = m_;
    double u = k*k/m/2., g = Gamma_;
    double beta = 1/(KB/AU_eV*temp_);
    // -------------
    // Gauss profile
    // -------------
    auto f = [g](double x) { return 1/g/sqrt(2*M_PI)*exp(-x*x/2./g/g); };
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
        fb += 7*f(x0-u) * eqFun_x(mu, x0) +
            32*f(xm1-u) * eqFun_x(mu, xm1) +
            12*f(xm2-u) * eqFun_x(mu, xm2) +
            32*f(xm3-u) * eqFun_x(mu, xm3) +
            7*f(x1-u) * eqFun_x(mu, x1);
    }
    return fb*2*h/4./45.;
}


/// Voigt * SF convolution
double WignerFunction::voigt(double k) {
    double mu = k > 0 ? uL_ : uR_;
    double m = m_;
    double u = k*k/m/2., g = Gamma_;
    double beta = 1/(KB/AU_eV*temp_);
    // Calculating eta - mixing parameter in pseudo-Voigt profile
    double gammaL = 2*Gamma_, gammaG = 2*sqrt(2*log(2))*Gamma_;
    double gamma = pow(pow(gammaG,5) +
        2.69269*pow(gammaG,4)*gammaL +
        2.42843*pow(gammaG,3)*pow(gammaL,2) +
        4.47163*pow(gammaG,2)*pow(gammaL,3) +
        0.07842*gammaG*pow(gammaL,4) + pow(gammaL,5), 0.2);
    double eta = 1.36603*(gammaL/gamma) - 0.47719*pow(gammaL/gamma,2) + 0.11116*pow(gammaL/gamma, 3);
    // ----------------------
    // (pseudo-)Voigt profile
    // ----------------------
    auto f = [eta, g](double x)->double{ 
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
        fb += 7*f(x0-u) * eqFun_x(mu, x0) +
            32*f(xm1-u) * eqFun_x(mu, xm1) +
            12*f(xm2-u) * eqFun_x(mu, xm2) +
            32*f(xm3-u) * eqFun_x(mu, xm3) +
            7*f(x1-u) * eqFun_x(mu, x1);
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
 * Calculates Fermi energy from el. density using Fermi integral
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


// ---------------
// Other functions
// ---------------


void WignerFunction::calc_IVchar(double v_min, double v_max, size_t nv){
    double dv = (v_max-v_min)/(nv-1), v = 0, curr, currD_range, carrNum;
    double uBias = uBias_;
    iv_v_.resize(nv), iv_i_.resize(nv);
    iv_iRange_.resize(nv), iv_n_.resize(nv);
    // double E;
    std::ofstream ivChar, vpMap;
    ivChar.open("OutData/ivChar.out", std::ios::out);
    vpMap.open("OutData/vpMap.out", std::ios::out);
    vpMap<<"# u_B [eV]  p [a.u]  1/4  1/2  3/4\n";
    ivChar<<"# u_B [eV]  J [Acm^{-2}]  |max(J)-min(J)|  N\n";
    // cout<<"# u_B [eV]  J [Acm^{-2}]  |max(J)-min(J)|  N\n";
    for (size_t i = 0; i < nv; ++i) {
        v = v_min+i*dv;
        // set_uBias(v);
        // solveWignerEq();
        solveWignerPoisson(v, 1e-5, 1, 1000, false);
        curr = calcCurr();
        currD_range = range(currD_);
        carrNum = calcNorm();
        calcCD_K();
        // E = i*dv/lD_ * AU_eV/1e3 / AU_cm;
        ivChar<<v*AU_eV<<' '<<curr*AU_Acm2<<' '<<currD_range*AU_Acm2<<' '<<carrNum/AU_cm2<<endl;
        // cout<<v*AU_eV<<' '<<curr*AU_Acm2<<' '<<currD_range*AU_Acm2<<' '<<carrNum/AU_cm2<<endl;  // <<' '<<calcNorm()/AU_cm2<<endl;
        iv_v_(i) = v, iv_i_(i) = curr, iv_iRange_(i) = currD_range, iv_n_(i) = carrNum;
        for (size_t j=0; j<nk_; ++j) {
                vpMap<<v*AU_eV<<' '<<k_(j)
                    <<' '<<cdK_(j)
                    // <<' '<<f_(size_t(nx_/4.),j)  // col. 3
                    // <<' '<<f_(size_t(nx_/2.),j)  // col. 4
                    // <<' '<<f_(size_t(nx_*3/4.),j)  // col. 5
                    <<'\n';
        }
        vpMap<<'\n';
    }
    ivChar.close();
    vpMap.close();
    arma::vec fit = arma::polyfit(iv_v_, iv_i_, 1);
    fit.print();
    uBias_ = uBias;
}


void WignerFunction::calcMobility() {

    arma::vec el_f(nx_, arma::fill::zeros);
    arma::vec mob(nx_, arma::fill::zeros);

    calcCD_X();
    calcCurr();
    arma::vec dnedx = calcFirstDer(cdX_, dx_);

    for (size_t i = 0; i < nx_; ++i)
        el_f(i) = -du_(i);

    for (size_t i = 0; i < nx_; ++i) {
        double m = cdX_(i)*el_f(i)+KB/AU_eV*temp_*dnedx(i);
        mob(i) = currD_(i)/m * AU_cm2/AU_eV/AU_s;
    }

    /*
    for (size_t i = 0; i < nx_; ++i)
        cout<<x_(i)<<' '<<mob(i)<<' '<<ne(i)<<' '<<dnedx(i) \
            <<' '<<u_(i)*AU_eV<<' '<<el_f(i)<<' '<<jn(i)<<' '<<mob(i)*KB*temp_<<endl;
    */

    // Average mobility
    double mob_av = 0;
    for (size_t i = 1; i < nx_; ++i)
        mob_av += (mob(i-1) + mob(i))*dx_/2./lD_;
    cout<<m_<<' '<<mob_av<<endl;

}
