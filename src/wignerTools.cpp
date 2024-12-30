#include "lib.hpp"
#include "WignerFunction.hpp"

/// Sets up equilibrium function
/// If read_data = true, reads equilibrium from input_file
void WignerFunction::setEquilibriumFunction(std::string input_file = "", bool read_data = false){
    // solveWignerEq();
    // calcCD_X();
    //  nx = cdX_;
    if (read_data) {
        fEq_.load(input_file);
        if (fEq_.size() != f_.size()) {
            cout<<"ERROR IN setEquilibriumFunction: fEq_.size() != f_.size()"<<endl;
            exit(0);
        }
    }
    else {
        double uBias = uBias_, rR = rR_;
        uBias_ = 0, rR_ = 0.;
        solveWignerEq();
        // calcCD_X();
        //  nx_0 = cdX_;
        for (size_t i=0; i<nx_; ++i)
            for (size_t j=0; j<nk_; ++j)
                fEq_(i,j) = f_(i,j);  // *nx(i)/nx_0(i)
        f_.zeros();
        uBias_ = uBias, rR_ = rR;
    }
}


/// Calculates current density
double WignerFunction::calcCurr() {
    currD_.zeros();  // Current density  // array<double>
    double a = 2., b = 1.;
    //  cur1 - current density in node i-1/2
    //  cur2 - current density in node i+1/2
    //  nc(nx_, arma::fill::zeros), nc_eq(nx_, arma::fill::zeros);
    // for (size_t i=0; i<nx_; ++i) {
    //     for (size_t j=1; j<nk_/2; ++j) {
    //         nc(i) += (f_(i,2*j-2)+4*f_(i,2*j-1)+f_(i,2*j))*dk_/3.;  // simpson
    //         nc_eq(i) += (fEq_(i,2*j-2)+4*fEq_(i,2*j-1)+fEq_(i,2*j))*dk_/3.;  // simpson
    //     }
    //     // if (bcType_ > 0) nc(i) /= 2.*M_PI, nc_eq(i) /= 2.*M_PI;
    // }
    // Current density is calculated in node i+1/2
    if (diffSch_J_ == "CDS1") {
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
    //
    // Metoda prostokątów
    // for (size_t i=0; i<nx_; ++i) {
    //     for (size_t j=0; j<nk_; ++j) {
    //         currD_(i) += k_(j)/m_*f_(i,j)*dk_;
    //         if (bcType_ > 0) currD_(i) /= 2.*M_PI;
    //     }
    // }
    //
    // Metoda trapezów / simpsona
    // for (size_t i=0; i<nx_; ++i) {
    //     for (size_t j=1; j<nk_/2; ++j)
    //         // cdX_(i) += (f_(i,j-1) + f_(i,j))*dk_/2.;  //  / 2./M_PI  // trapezoid
    //         currD_(i) +=
    //             ( k_(2*j-2)*f_(i,2*j-2) +
    //             4*k_(2*j-1)*f_(i,2*j-1) +
    //             k_(2*j)*f_(i,2*j) )*dk_/3.;  // simpson
    //     if (bcType_ > 0) currD_(i) /= 2.*M_PI;
    // }
    //
    // Armadillo function
    // arma::mat f = f_;
    // f.each_row() %= k_.t()/m_;
    // currD_ = trapz(k_, f, 1);
    if (bcType_ > 0) currD_ /= 2.*M_PI;
    return sum(currD_)*dx_/l_;
}


arma::vec WignerFunction::calcCD_X(){
    // Calculates carrier density in x space
    cdX_.zeros();
    for (size_t i=nx_; i--;) {
        for (size_t j=1; j<nk_/2; ++j)
            // cdX_(i) += (f_(i,j-1) + f_(i,j))*dk_/2.;  //  / 2./M_PI  // trapezoid
            cdX_(i) += (f_(i,2*j-2)+4*f_(i,2*j-1)+f_(i,2*j))*dk_/3.;  // simpson
    }
    // arma::mat f = f_.t();
    // for (size_t i=nx_; i--;) {
    //     for(arma::mat::col_iterator j = f.begin_col(i)+1; j != f.end_col(i); ++j)
    //         cdX_(i) += ( *j + *(j-1) ) * dk_/2.;
    //     // cdX_(i) += ( f_(i,0) + f_(i,1) ) * dk_/2.;
    // }
    if (bcType_ > 0) cdX_ /= 2.*M_PI;
    return cdX_;
}


arma::vec WignerFunction::calcCD_K(){
    // Calculates carrier density in k space
    //  cd(nk_, arma::fill::zeros);  // array<double>
    cdK_.zeros();
    for (size_t j=nk_; j--;) {
        for (size_t i=1; i<nx_/2; ++i)
            // cd(j) += (f_(i-1,j) + f_(i,j))*dx_/2.;  // trapezoid
            cdK_(j) += (f_(2*i-2,j)+4*f_(2*i-1,j)+f_(2*i,j))*dx_/3.;  // simpson
    }
    // for (size_t j=nk_; j--;) {
    //     for(arma::mat::col_iterator i = f_.begin_col(j)+1; i != f_.end_col(j); ++i)
    //         cdK_(j) += ( *i + *(i-1) ) * dk_/2.;
    //     // cdK_(j) += ( f_(0,j) + f_(1,j) ) * dk_/2.;
    // }
    return cdK_;
}


double WignerFunction::calcNorm(){
    // Calculates density function norm (integral over whole space)
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


double WignerFunction::calcSDK(){
    double ev = 0, ev2 = 0, s_ev, s_ev2;
    /*
    double sum_x = 0, sum_k;
    // E[k]
    sum_x = 0;
    for (size_t i=0; i<nx_; ++i) {
        sum_k = 0;
        for (size_t j=1; j<nk_; ++j)
            sum_k += ( f_(i,j) + f_(i,j-1) ) * k_(j) * dk_/2.;
        if (i == 0 || i == nx_-1)
            sum_x += sum_k * dx_/2.;
        else
            sum_x += sum_k * dx_;
    }
    ev = sum_x/2./M_PI;  //  / calcNorm()
    // E[k^2]
    sum_x = 0;
    for (size_t i=0; i<nx_; ++i) {
        sum_k = 0;
        for (size_t j=1; j<nk_; ++j)
            sum_k += ( f_(i,j) + f_(i,j-1) ) * k_(j) * k_(j) * dk_/2.;
        if (i == 0 || i == nx_-1)
            sum_x += sum_k * dx_/2.;
        else
            sum_x += sum_k * dx_;
    }
    ev2 = sum_x/2./M_PI;  //  / calcNorm()
    */
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


/// Sets system potential as linear function
void WignerFunction::setPotBias(double uBias) {
    uBias_ = uBias;
    uC_.zeros();
    double x;
    // for (size_t i = 0; i < nx_; ++i) {
    //     x = x_(i)-lC_;
    //     if (x<0)
    //         uC_(i) = 0;
    //     else if (x>=0 && x<=lD_)
    //         uC_(i) = -uBias_*x/lD_;
    //     else if (x>lD_)
    //         uC_(i) = -uBias_;
    // }
    for (size_t i = 0; i < nx_; ++i) {
        x = x_(i);
        uC_(i) = uBias_*(0.5-x/l_);
    }
}  // End of setPotBias

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
    // double s2 = 2*gwp_dx_*gwp_dx_;
    double sx = 2*gwp_dx*gwp_dx, sk = 2*gwp_dk*gwp_dk;
    for (size_t i=0; i<nx_; ++i) {
        for (size_t j=0; j<nk_; ++j)
            f_(i,j) += exp(
                - (k_(j)-gwp_k0)*(k_(j)-gwp_k0)/sk
                - (x_(i)-gwp_x0)*(x_(i)-gwp_x0)/sx  ) * A; // * gwp_A_, / M_PI
    }
}


//
//  WP evolution (no potential)
//

double WignerFunction::wavePacket_TEV(double gwp_x0, double gwp_dx,
    double gwp_k0, double gwp_dk, double x, double k)
    { return exp( -(k-gwp_k0)*(k-gwp_k0)/2./gwp_dk/gwp_dk
        -(x-gwp_x0)*(x-gwp_x0)/2./gwp_dx/gwp_dx ) / M_PI; }


//
// Boundary conditions
//

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
            /// TODO: Armadillo convolution
            // else if (bcType_ == 4)
            //      a = (nk_, arma::fill::zeros), b = (nk_, arma::fill::zeros), c;
            //     for (size_t j=0; j<nk_; ++j) {
            //         a(j) = supplyFunction(k_(j));
            //         b(j) = lorentz(k_(j));
            //     }
            //     c = conv(a, b);
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
    // double c = 4*M_PI*m*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);     // [au]
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
    // double mu = k > 0 ? uL_ : uR_;
    double m = m_;
    double c = cD_*pow(2*M_PI*m*KB/AU_eV*temp_, -3./2.);
    double ex = exp(-k*k/2/2./(KB/AU_eV*temp_));
    // double v = k/m;
    // double c = pow(2.*M_PI*m*(KB/AU_eV*temp_), 0.5);  // * cD_
    // double ex = exp(-(m*v*v/2.-mu)/(KB/AU_eV*temp_));     // [au]
    // double c = 4*M_PI*m*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);     // [au]
    return c * ex;
}

/// Gaussian with sigma parameter
double WignerFunction::gaussian_bc(double k){
    double m = m_;
    // double pF = k > 0 ? sqrt(uL_*2.*m) : sqrt(uR_*2.*m);
    double kBT = KB/AU_eV*temp_;
    double sigma = 1;
    // double c = 1./sqrt(2.*M_PI*m*(KB/AU_eV*temp_)*sigma*sigma);
    // double ex = exp(-(k*k/m/2.-mu)/(KB/AU_eV*temp_)/sigma/sigma);     // [au]
    // double c = 4*M_PI*m*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);     // [au]
    // double c = 1./sqrt(2.*M_PI*m*kBT)/sigma*cD_;
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


double WignerFunction::eqFun_x(double mu, double x) {
    // Equilibrium function for electrons in contact
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
        // double c = 4*M_PI*m*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);     // [au]
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
    // Gauss profile
    auto f = [g](double x) { return 1/g/sqrt(2*M_PI)*exp(-x*x/2./g/g); };
    // Equilibrium function for electrons in contact
    // auto f_eq = [mu, bcType, gaussian_x, sf_x](double x) {
    //     if (bcType < 0)
    //         return gaussian_x(mu, x);
    //     else
    //         return sf_x(mu, x);
    // };
    // convolution
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
    // (pseudo-)Voigt profile
    auto f = [eta, g](double x)
        { return eta*( g/(x*x+g*g)/M_PI ) +  // lorentz
            (1-eta)*( 1/g/sqrt(2*M_PI)*exp(-x*x/2./g/g) ); };  // gauss
    // Equilibrium function for electrons in contact
    // auto f_eq = [mu, bcType, gaussian_x, sf_x](double x) {
    //     if (bcType < 0)
    //         return gaussian_x(mu, x);
    //     else
    //         return sf_x(mu, x);
    // };
    // convolution
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


//
// FERMI INTEGRAL
//


/** 
 * Density of states
 * @param m electron effective mass
 * @param T temperature
*/
inline double nC(double m, double T)
    { return 2.*pow(m*KB*T/AU_eV/2./M_PI, 3./2.); }  // [au]


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
 *Calc Fermi integral
*/
double calcFermiEn(double n0, double m, double T) {
    // Calculates Fermi energy from el. density using Fermi integral
    double eta = 0., f = 0.;
    double x = n0/nC(m, T);  // For GaAs, 300 K -> x = 4.6
    //  Looking for eta
    if (x > fermiInt(.5, 5.)) {
        eta = pow(x*std::tgamma(5./2.), 2./3.);
        /*
        eta = log(x) + 0.25*pow(2., .5)*x \
            + (3/16.-1/9.*pow(3., .5))*x*x \
            + (1/8.+5/48.*pow(2., .5)-1/9.*pow(6., .5))*x*x*x \
            + (1585/6912+5/32.*pow(2., .5)-5/24.*pow(3., .5)-1/25.*pow(5., .5))*x*x*x*x; */
        // eta = log(x) + 3.53553E-1*x - 4.95009E-3*pow(x, 2.) + 1.48386E-4*pow(x, 3.) - 4.42563E-6*pow(x, 4.);
        /*
        while (f < x) {
            eta += 0.01;
            f = pow(eta, 3./2.)/std::tgamma(5/2);
            }
        */
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
    // cout<<"eta = "<<eta/(KB*T)<<", x = "<<x<<endl;
    return eta/AU_eV;
}


//
// Other functions
//


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
