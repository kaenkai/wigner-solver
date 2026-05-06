#ifndef Poisson1D_HPP
#define Poisson1D_HPP

#include <armadillo>

using namespace arma;

class Poisson1D {

    public:

        Poisson1D(size_t nx, double h) : nx_ (nx), h_ (h) {
            rho_ = arma::vec(nx_, arma::fill::zeros);
            nE_ = arma::vec(nx_, arma::fill::zeros);
            uOld_ = arma::vec(nx_, arma::fill::zeros);
            uNew_ = arma::vec(nx_, arma::fill::zeros);
            du_ = arma::vec(nx_, arma::fill::zeros);
            dPu_ = arma::sp_mat(nx_, nx_);
            pFun_ = arma::vec(nx_, arma::fill::zeros);
            epsilonR_ = 1, temp_ = 300;
            dirichletL_ = 0, dirichletR_ = 0;
        }
        Poisson1D() : nx_ (100), h_ (1.) {
            rho_ = arma::vec(nx_, arma::fill::zeros);
            nE_ = arma::vec(nx_, arma::fill::zeros);
            uOld_ = arma::vec(nx_, arma::fill::zeros);
            uNew_ = arma::vec(nx_, arma::fill::zeros);
            du_ = arma::vec(nx_, arma::fill::zeros);
            dPu_ = arma::sp_mat(nx_, nx_);
            pFun_ = arma::vec(nx_, arma::fill::zeros);
            epsilonR_ = 1, temp_ = 300;
            dirichletL_ = 0, dirichletR_ = 0;
        }
        ~Poisson1D(){};

        void solve();
        void solve_tridiag();

        void set_boundary_conditions(double dirichletL, double dirichletR) {
            this -> dirichletL_ = dirichletL;
            this -> dirichletR_ = dirichletR;
        };
        void set_epsilonR(double epsilonR) {epsilonR_ = epsilonR;};
        void set_temp(double temp) {temp_ = temp;};

        double get_h() {return h_;};
        size_t get_nx() {return nx_;};
        double get_epsilonR() {return epsilonR_;};
        double get_temp() {return temp_;};
        double get_dirichletL() {return dirichletL_;};
        double get_dirichletR() {return dirichletR_;};


        /*
        Poisson equasion tests
        */
        static void testUniformCharge();
        static void testExponentCharge();
        static void testSineCharge();

        double testSine();
        static void testGrid();
        static void testSelfConsistency();

        static void testChargedPlane();

        arma::vec rho_;     // Charge density
        arma::vec nE_;      // Electron density
        arma::vec uOld_;    // Potential from previous iteration
        arma::vec uNew_;    // Potential from current iteration
        arma::vec du_;      // Potential difference between iterations

        arma::vec pFun_;    // Poisson function (right-hand side of the linearized equation)
        arma::sp_mat dPu_;  // Derivative of Poisson function with respect to potential

    private:
        size_t nx_;                         // number of grid points
        double h_;                          // grid step
        double dirichletL_, dirichletR_;    // Dirichlet boundary conditions
        double epsilonR_, temp_;            // relative permitivitty, temperature
};

# endif
