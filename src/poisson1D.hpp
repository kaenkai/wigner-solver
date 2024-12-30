#ifndef POISSON1D_HPP
#define POISSON1D_HPP

#include "lib.hpp"

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
            epsilonR_ = 1, temp_ = 300, beta_ = 1;
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
            epsilonR_ = 1, temp_ = 300, beta_ = 1;
            dirichletL_ = 0, dirichletR_ = 0;
        }
        ~Poisson1D(){};

        void solve();
        void solve_gummel();
        void solve_tridiag();

        void set_boundary_conditions(double dirichletL, double dirichletR) {
            dirichletL_ = dirichletL;
            dirichletR_ = dirichletR;
        };
        void set_epsilonR(double epsilonR) {epsilonR_ = epsilonR;};
        void set_temp(double temp) {temp_ = temp;};
        void set_beta(double beta) {beta_ = beta;};

        double get_h() {return h_;};
        size_t get_nx() {return nx_;};
        double get_epsilonR() {return epsilonR_;};
        double get_temp() {return temp_;};
        double get_beta() {return beta_;};
        double get_dirichletL() {return dirichletL_;};
        double get_dirichletR() {return dirichletR_;};

        static void testPoisson();

        arma::vec rho_;
        arma::vec nE_;
        arma::vec uOld_;
        arma::vec uNew_;
        arma::vec du_;

        arma::vec pFun_;
        arma::sp_mat dPu_;

    private:
        size_t nx_;
        double h_;
        double dirichletL_, dirichletR_;
        double epsilonR_, temp_;
        double beta_;  // potential mixing
        
};

# endif
