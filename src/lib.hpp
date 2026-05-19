#ifndef LIB_H
#define LIB_H

#include <vector>
#include <string>
#include <omp.h>

#define ARMA_USE_SUPERLU 1
#include <armadillo>


namespace AU {
    // ---------
    // constants
    // ---------
    double const KB {8.617333262145179E-05};  // Boltzmann constant [eV/K]
    double const KB_J {1.380649E-23};  // Boltzmann constant [J/K]
    double const E0 {1.602176634E-19};  // Elementary charge [C]
    double const HBAR_J {1.0545718176461565E-34};  // Reduced Planck constant [Js]
    double const HBAR_eV {6.582119569509067E-16};  // Reduced Planck constant [eVs]
    double const M0 {9.1093837015E-31};  // Electron mass [kg]
    double const EPS0 {8.8541878128E-12};  // Electric constant

    // --------------------------------------------------
    // atomic units
    // source: https://en.wikipedia.org/wiki/Atomic_units
    // --------------------------------------------------
    double const eV {27.211386245981};  // Hartree energy [eV]
    double const nm {0.0529177210544};  // Bohr radius [nm]
    double const m {nm*1E-9};  // Bohr radius [m]
    double const m2 {m*m};  // [cm**2]  (nm*1e-7)**2
    double const m3 {m*m*m};  // [cm**3]  (nm*1e-7)**3
    double const cm {m*1e2};  // [cm]
    double const cm2 {cm*cm};  // [cm**2]  (nm*1e-7)**2
    double const cm3 {cm*cm*cm};  // [cm**3]  (nm*1e-7)**3
    double const s {HBAR_eV/eV};  // Time [s]
    double const A {E0/s};  // [A]  _e0/_tau0
    double const Acm2 {A/cm2};  // [A/cm**2]  A/cm/cm
}


// ------------------
// default parameters
// ------------------
double const A_GaAs {0.565};  // GaAs lattice constant [nm]
double const EPS_GaAs {13.1};  // GaAs relative permittivity
double const M_GaAs {0.067};  // GaAs effective mass
double const TEMP {300};
double const DIFF {1};  // Diffusion coefficient [cm^2/s]
double const ND {2E18};  // Donor concentration [cm^-3]


/**
 * array class
 * @deprecated now armadillo vec is used instead
 */
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


/**
 * matrix class
 * @deprecated now armadillo mat is used instead
 */
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




// ----------------------
// Functions declarations
// ----------------------

double calcInt(arma::vec, double);
arma::vec calcFirstDer(arma::vec, double);
arma::vec calcSecondDer(arma::vec, double);
arma::vec calcThirdDer(arma::vec, double);
arma::vec normalDistribution(double, double, double, double, size_t);
double calcFermiEn(double, double, double);
std::map<std::string, double> readParameters(std::string);
void saveMatGP(arma::mat, std::string);

void testDerivatives();
void testBTE();


#endif
