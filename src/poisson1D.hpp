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
            neumannL_ = 0, neumannR_ = 0;
            pBC_D_ = true, pBC_vN_ = false;
        }
        ~Poisson1D(){};

        void solve();
        void solve_gummel();
        void solve_tridiag();

        void testPoisson();

        size_t nx_;
        double h_;
        double dirichletL_, dirichletR_;
        double neumannL_, neumannR_;
        bool pBC_D_, pBC_vN_;
        double epsilonR_, temp_;
        double beta_;  // potential mixing

        arma::vec rho_;
        arma::vec nE_;
        arma::vec uOld_;
        arma::vec uNew_;
        arma::vec du_;

        arma::vec pFun_;
        arma::sp_mat dPu_;

};

# endif
