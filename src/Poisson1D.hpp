#ifndef Poisson1D_HPP
#define Poisson1D_HPP

#include "lib.hpp"


class Poisson1D {

    public:

        Poisson1D(size_t nx, double h) : nx_ (nx), h_ (h) {
            rho_ = arma::vec(nx_, arma::fill::zeros);
            uOld_ = arma::vec(nx_, arma::fill::zeros);
            uNew_ = arma::vec(nx_, arma::fill::zeros);
            du_ = arma::vec(nx_, arma::fill::zeros);
            epsilonR_ = 1, temp_ = 300;
            dirichletL_ = 0, dirichletR_ = 0;
        }
        Poisson1D() : nx_ (100), h_ (1.) {
            rho_ = arma::vec(nx_, arma::fill::zeros);
            uOld_ = arma::vec(nx_, arma::fill::zeros);
            uNew_ = arma::vec(nx_, arma::fill::zeros);
            du_ = arma::vec(nx_, arma::fill::zeros);
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
        void set_uOld(arma::vec uOld) {uOld_ = uOld;};  // todo: size check
        void set_rho(arma::vec rho) {rho_ = rho;};

        double get_h() {return h_;};
        size_t get_nx() {return nx_;};
        double get_epsilonR() {return epsilonR_;};
        double get_temp() {return temp_;};
        double get_dirichletL() {return dirichletL_;};
        double get_dirichletR() {return dirichletR_;};
        arma::vec get_uNew() {return uNew_;};
        arma::vec get_uOld() {return uOld_;};
        arma::vec get_du() {return du_;};
        arma::vec get_rho() {return rho_;};

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

    private:
        size_t nx_;                         // number of grid points
        double h_;                          // grid step
        double dirichletL_, dirichletR_;    // Dirichlet boundary conditions
        double epsilonR_, temp_;            // relative permitivitty, temperature
        arma::vec rho_;                     // Charge density
        arma::vec uOld_;                    // Potential from previous iteration
        arma::vec uNew_;                    // Potential from current iteration
        arma::vec du_;                      // Potential difference between iterations
};

# endif
