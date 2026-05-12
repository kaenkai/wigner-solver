#include "lib.hpp"


/**
 * Derivative test
 * @todo better differentiation scheme for third order derivative
 */
void testDerivatives() {
    const size_t nx = 101; // grid size
    double len = 2*M_PI;
    double h = len/(nx-1); // step size
    int k=2;  // number of oscillations

    arma::vec x = arma::linspace(0, len, nx);
    arma::vec f = arma::linspace(0, len, nx).transform(
        [k](double xx) -> double {return std::sin(k*xx);}
    ) + arma::vec(nx, arma::fill::randn)*1E-4;

    arma::vec f1_err = calcFirstDer(f, h)/k -
        arma::linspace(0, len, nx).transform(
            [k](double xx) -> double {
                return std::cos(k*xx);
            });

    arma::mat out_data;
    out_data.insert_cols(0, x);
    out_data.insert_cols(1, f);
    out_data.insert_cols(2, calcFirstDer(f, h)/k);
    out_data.insert_cols(3, -calcSecondDer(f, h)/k/k);
    out_data.insert_cols(4, -calcThirdDer(f, h)/k/k/k);
    out_data.insert_cols(5, f1_err);
    out_data.print();
}