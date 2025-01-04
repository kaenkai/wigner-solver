#ifndef LIB_H
#define LIB_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <ctime>
#include <map>
#include <iomanip>  // std::setw

#include <omp.h>

#define ARMA_USE_SUPERLU 1
// #define ARMA_OPENMP_THREADS 2
// #define ARMA_PRINT_ERRORS 1
// #define ARMA_WARN_LEVEL 3
// #define OMP_NUM_THREADS 2

#include <armadillo>

// using namespace std;
using std::cout;
using std::endl;
using std::setw;
// using namespace arma;

// ############################## GLOBAL VARIABLES ##############################
// #define N_THREADS 4

// ############################## CONSTANTS ##############################
/// TODO: put into a namespace (maybe)

double const PI {M_PI};
double const KB {8.617333262145179E-05};  // Boltzmann constant [eV/K]
double const KB_J {1.380649E-23};  // Boltzmann constant [J/K]
double const E0 {1.602176634E-19};  // Elementary charge [C]
double const HBAR_J {1.0545718176461565E-34};  // Reduced Planck constant [Js]
double const HBAR_eV {6.582119569509067E-16};  // Reduced Planck constant [eVs]
double const M0 {9.1093837015E-31};  // Electron mass [kg]
double const EPS0 {8.8541878128E-12};  // Electric constant

// ############################## SYSTEM CONSTANTS ##############################

double const A_GaAs {0.565};  // GaAs lattice constant [nm]
double const EPS_GaAs {13.1};  // GaAs relative permittivity
double const M_GaAs {0.067};  // GaAs effective mass
double const TEMP {300};
double const DIFF {1};  // Diffusion coefficient [cm^2/s]
double const ND {2E18};  // Donor concentration [cm^-3]

// ############################## atomic units ##############################

double const AU_eV {27.211386245988};  // Hartree energy [eV]
double const AU_nm {0.0529177210903};  // Bohr radius [nm]
double const AU_m {AU_nm*1E-9};  // Bohr radius [m]
double const AU_m2 {AU_m*AU_m};  // [cm**2]  (AU_nm*1e-7)**2
double const AU_m3 {AU_m*AU_m*AU_m};  // [cm**3]  (AU_nm*1e-7)**3
double const AU_cm {AU_m*1e2};  // [cm]
double const AU_cm2 {AU_cm*AU_cm};  // [cm**2]  (AU_nm*1e-7)**2
double const AU_cm3 {AU_cm*AU_cm*AU_cm};  // [cm**3]  (AU_nm*1e-7)**3
double const AU_s {HBAR_eV/AU_eV};  // Time [s]
double const AU_A {E0/AU_s};  // [A]  _e0/_tau0
double const AU_Acm2 {AU_A/AU_cm2};  // [A/cm**2]  AU_A/AU_cm/AU_cm

// double const AU_m {4.*M_PI*EPS0*HBAR_J*HBAR_J/M0/E0/E0};  // Bohr radius [m] also 1./(4.*M_PI*EPS0*AU_eV)
// double const AU_eV {HBAR_J*HBAR_J/M0/AU_m/AU_m*1/E0 };  // Hartree energy [eV] (*1/E0: J -> eV)
// double const AU_s {HBAR_eV/AU_eV};  // Time [s]
// double const AU_nm {AU_m*1E9};  // Bohr radius [nm]
// double const AU_cm {AU_m*1e2};  // [cm]
// double const AU_cm2 {AU_cm*AU_cm};  // [cm**2]  (AU_nm*1e-7)**2
// double const AU_cm3 {AU_cm*AU_cm*AU_cm};  // [cm**3]  (AU_nm*1e-7)**3
// double const AU_A {E0/AU_s};  // [A]  _e0/_tau0
// double const AU_Acm2 {AU_A/AU_cm2};  // [A/cm**2]  AU_A/AU_cm/AU_cm

// ############################## atomic units - numberic ##############################

/*
double const AU_eV {27.211386245988};  // Hartree energy [eV]
double const AU_nm {0.0529177210903};  // Bohr radius [nm]
double const AU_m {AU_nm*1E-9};  // Bohr radius [m]
double const AU_cm {AU_m*1e2};  // [cm]
double const AU_cm2 {AU_cm*AU_cm};  // [cm**2]  (AU_nm*1e-7)**2
double const AU_cm3 {AU_cm*AU_cm*AU_cm};  // [cm**3]  (AU_nm*1e-7)**3
double const AU_s {2.4188843265857225E-17}; // Time [s]
double const AU_A {6.623618237509881E-3};  // [A]
double const AU_Acm2 {2.365337010944052E14};  // [A/cm**2]
*/

// ############################## atomic units - simplified ##############################

// double const AU_nm {0.053};  // Bohr radius [nm]
// double const AU_m {AU_nm*1E-9};  // Bohr radius [m]
// double const AU_eV {27.21};  // Hartree energy [eV]
// double const AU_s {2.42E-17};  // Time [s]
// double const AU_cm {AU_nm*1e-7};  // [cm]  AU_nm*1e-7
// double const AU_cm2 {AU_cm*AU_cm};  // [cm**2]  (AU_nm*1e-7)**2
// double const AU_cm3 {AU_cm*AU_cm*AU_cm};  // [cm**3]  (AU_nm*1e-7)**3
// double const AU_A {6.62E-3};  // [A]  _e0/_tau0
// double const AU_Acm2 {AU_A/AU_cm2};  // [A/cm**2]  AU_A/AU_cm/AU_cm

// ############################## scaling factors ##############################

// double const AU_m {sqrt(EPS0*EPS_GaAs*KB_J*TEMP/E0/E0/(ND*1E6))};  // Debye length [m] (*1e6: cm->m)
// double const AU_eV {KB*TEMP};  // Thermal energy [eV]
// double const AU_nm {AU_m*1E9};  // Bohr radius [nm]
// double const AU_cm {AU_m*1e2};  // [cm]
// double const AU_cm2 {AU_cm*AU_cm};  // [cm**2]  (AU_nm*1e-7)**2
// double const AU_cm3 {AU_cm*AU_cm*AU_cm};  // [cm**3]  (AU_nm*1e-7)**3
// double const AU_s {AU_cm*AU_cm/DIFF};  // Time [s]
// double const AU_A {E0/AU_s};  // [A]  _e0/_tau0
// double const AU_Acm2 {AU_A/AU_cm2};  // [A/cm**2]  AU_A/AU_cm/AU_cm

// #################### klasa array ####################
template <class T>
class array {
    private:
        std::vector<T> a;
        size_t n_row;
    public:
        array() {}
        array(size_t n) : a (std::vector<T>(n)), n_row (n) {}
        array(size_t n, T v) : a (std::vector<T>(n, v)), n_row (n) {}
        ~array() { a.clear(); }
        T& operator () (size_t i) { return a.at(i); }
        // Vector addition
        array<T> operator-(array<T> ai) {
            array<T> af(n_row);
            for (size_t i=0; i<n_row; ++i) af.at(i) = a.at(i)-ai.at(i);
            return af;
        }
        // Vector substraction
        array<T> operator+(array<T> ai) {
            array<T> af(n_row);
            for (size_t i=0; i<n_row; ++i) af.at(i) = a.at(i)+ai.at(i);
            return af;
        }
        // Vector multiplication
        array<T> operator*(array<T> ai) {
            array<T> af(n_row);
            for (size_t i=0; i<n_row; ++i) af.at(i) = a.at(i)*ai.at(i);
            return af;
        }
        // Vector division
        array<T> operator/(array<T> ai) {
            array<T> af(n_row);
            for (size_t i=0; i<n_row; ++i) af.at(i) = a.at(i)/ai.at(i);
            return af;
        }
        // Multiplying vector by scalar
        // array<T> operator*(Q s) {
        //   array<T> af(n_row);
        //   for (size_t i=0; i<n_row; ++i) af.at(i) = s*a.at(i);
        //   return af;
        // }
        T& at(size_t i) { return a.at(i); }
        size_t size() { return a.size(); }
        void add(T x) { a.push_back(x); }
        void copy(array<T> ac) { for (size_t i=0; i<n_row; ++i) a.at(i) = ac(i); }
        void zeros() { for (size_t i=0; i<n_row; ++i) a.at(i) = 0.; }
        double max() { return *max_element(a.begin(), a.end()); }
        double min() { return *min_element(a.begin(), a.end()); }
        double sum() {
            double sum = 0;
            for (size_t i=0; i<n_row; ++i) sum += a.at(i);
            return sum;
        }
};

// #################### klasa matrix ####################
template <class T>
class matrix {
    private:
        array< array<T> > m;
        size_t n_row, n_col;
    public:
        matrix() {}
        matrix(size_t n1, size_t n2) : n_row (n1), n_col (n2)
                { for (size_t i=0; i<n_row; ++i) m.add( array<T>(n_col) ); }
        ~matrix() {}
        array<T>& operator () (size_t i) { return m.at(i); }
        T& operator () (size_t i, size_t j) { return m.at(i).at(j); }
        array<T>& at(size_t i) { return m.at(i); }
        matrix<T> operator-(matrix<T> mi) {
            matrix<T> mf(n_row, n_col);
            for (size_t i=0; i<n_row; ++i)
                for (size_t j=0; j<n_col; ++j)
                    mf.at(i) = m.at(i)-mi.at(i);
            return mf;
        }
        matrix<T> operator+(matrix<T> mi) {
            matrix<T> mf(n_row, n_col);
            for (size_t i=0; i<n_row; ++i)
                for (size_t j=0; j<n_col; ++j)
                    mf.at(i) = m.at(i)+mi.at(i);
            return mf;
        }
        T& at(size_t i, size_t j) { return m.at(i).at(j); }
        size_t size() { return n_col*n_row; }
        void copy(matrix<T> mc) { for (size_t i=0; i<n_row; ++i) m(i).copy(mc(i)); }
        // void add(T x){ m.push_back(x); }
        void zeros() { for (size_t i=0; i<n_row; ++i) m(i).zeros(); }
};


//
// Funkcje do liczenia całek i pochodnych
// Definicje tych funkcji znajdują się w pliku wignerTools.cpp
//


/// Integral (Simpson rule)
template <class T>
double calcInt(arma::vec f, T h){
    size_t n = f.size();
    T ig = 0;
    for (size_t i=1; i<n/2; ++i)
        // ig += (f(i-1) + f(i))*h/2.;  // trapezoid
        ig += (f(2*i-2)+4*f(2*i-1)+f(2*i))*h/3.;  // simpson
    return ig;
}


/// First derivative
template <class T>
arma::vec calcFirstDer(arma::vec f, T h){
    size_t n = f.size();
    arma::vec df(n, arma::fill::zeros);
    for (size_t i=2; i<n-2; ++i)
        df(i) = (1/12.*f(i-2)-2/3.*f(i-1)+2/3.*f(i+1)-1/12.*f(i+2))/h;
        // df(i) = (-f(i-1)+f(i+1))/h/2.;
    df(0) = (-3.*f(0)+4.*f(1)-f(2))/h/2.;
    df(n-1) = (3.*f(n-1)-4.*f(n-2)+f(n-3))/h/2.;
    df(1) = (-f(0)+f(2))/h/2.;
    df(n-2) = (-f(n-3)+f(n-1))/h/2.;
    // df(0) = (-f(0)+f(1))/h;
    // df(n-1) = (f(n-1)-f(n-2))/h;
    // df(0) = df(1);
    // df(n-1) = df(n-2);
    return df;
}


/// Second derivative
template <class T>
arma::vec calcSecondDer(arma::vec f, T h){
    size_t n = f.size();
    arma::vec df(n, arma::fill::zeros);
    for (size_t i=1; i<n-1; ++i)
        df(i) = (f(i-1)-2.*f(i)+f(i+1))/h/h;
    df(0) = (2*f(0)-5.*f(1)+4*f(2)-f(3))/h/h;
    df(n-1) = (2*f(n-1)-5.*f(n-2)+4*f(n-3)-f(n-4))/h/h;
    return df;
}


/// Third derivative
template <class T>
arma::vec calcThirdDer(arma::vec f, T h){
    size_t n = f.size();
    arma::vec df(n, arma::fill::zeros);
    for (size_t i=2; i<n-2; ++i)
        df(i) = (-f(i-2)/2.+f(i-1)-f(i+1)+f(i+2)/2.)/h/h/h;
    df(0) = (-f(0)+3.*f(1)-3*f(2)+f(3))/h/h/h;
    df(1) = (-f(1)+3.*f(2)-3*f(3)+f(4))/h/h/h;
    df(n-1) = (f(n-1)-3.*f(n-2)+3.*f(n-3)-f(n-4))/h/h/h;
    df(n-2) = (f(n-2)-3.*f(n-3)+3.*f(n-4)-f(n-5))/h/h/h;
    return df;
}

/** 
 * @brief Calculates gaussian distribution centered at (x_min + x_max) / 2
 * @param x_min Minimum value of x
 * @param x_max Maximum value of x
 * @param sig   Standard deviation
 * @param n     Number of points
 * @param A     Amplitude
 * @return arma::vec Gaussian distribution
*/
template <class T>
arma::vec gaussian_dist(T x_min = 0., T x_max = 1., T sig = 1., T A = 1., size_t n = 100){
    arma::vec gauss(n, arma::fill::zeros);
    double x = 0;
    double mu = (x_min + x_max)*0.5;
    for (size_t i = 0; i < n; ++i){
        x = x_min + i*(x_max-x_min)/(n-1);
        gauss(i) = exp(-(x-mu)*(x-mu)/2./sig/sig) * A;
    }
    return gauss;
}

double calcFermiEn(double, double, double);
std::map<std::string, double> readParam(std::string);


#endif
