#include "lib.hpp"
#include "WignerFunction.hpp"
#include "Poisson1D.hpp"

using namespace AtomicUnits;


/**
 * Solves Wigner-Poisson transport equation system
 * @todo put into separate class
 */
void WignerFunction::solveWignerPoisson
	(double u_bias, double i_alpha, double beta, size_t i_n_max, bool timeDependent){

    arma::vec j0(nx_, arma::fill::zeros), j1(nx_, arma::fill::zeros);
    arma::vec dj(nx_, arma::fill::zeros);
    arma::vec dRho(nx_, arma::fill::zeros);
    arma::vec u_der(nx_, arma::fill::zeros);
    arma::vec rho_old(nx_, arma::fill::zeros);
    arma::vec currD_a, nc_a, q_a, dJ_x_a, dU_x_a, dRho_x_a;
    arma::vec rho_test(nx_, arma::fill::zeros);

	//
	// Convergence criteria
    size_t n_max = i_n_max, n_min = 2;
    size_t n_dJ = 0, n_dU = 0, n_conv = 2;
    double max_dJ = 1e-6, max_dU = 1e-6;
    bool conv = false;
    double curr = 0, nc = 0, nd = 0, q = 0;  // Current, carrier nr, dopant nr, total charge
    double dJ_x = 0, dU_x = 0, dRho_x;  // Maximum and minimum values
    double alpha = i_alpha;  // Density mixing parameter

	//
	// Setting up Poisson solver
    Poisson1D p(nx_, dx_);  // Setting up Poisson solver
    p.set_boundary_conditions(u_bias/2., -u_bias/2.);  // Dirichlet BC
    p.set_epsilonR(epsilonR_), p.set_temp(temp_);   // Permittivity and temperature
	// p.uNew_ = uStart_, p.uOld_ = uStart_;  // Starting potential values
    p.uNew_.zeros(), p.uOld_.zeros(); 
    p.set_beta(beta);  // Potential mixing parameter

	//
	// Doping profile
    set_doping_profile(0.01);

	//
	// Initial conditions
    setBoundCond();
    p.rho_.zeros();
    p.solve();
	// Start electron concentration
	// uC_ = uStart_;
    uC_ = p.uNew_;
    solveWignerEq();
    calcCD_X();
    rho_old = nD_ - cdX_;
	// f_.zeros();
	// // Creating starting density function with condition that int(f_*dj) = nD_(i)
	// for (size_t i=0; i<nx_; ++i) {
	// 	// cout<<i*dx_*AU_nm*1e-3<<'\t'<<lC_*AU_nm*1e-3<<endl;
	//     for (size_t j=0; j<nk_; ++j) {
	// 		    f_(i,j) = bc_(j); // *nD_(i)/cD_
	// 		}
	// }

	/// TODO: Pointers to functions

    std::ofstream poisson_step("OutData/poisson_step.out");
    poisson_step<<"it,x [nm],rho [cm^{-3}],uNew [eV],{/Symbol d}u [eV],du/dx [au],J [Acm^{-2}]\n";

    std::ofstream tr_char("OutData/tr_char.csv");
    if ( timeDependent ) {
	    tr_char<<"t[fs],J [au],n_E [au],q [au],max({/Symbol d}J/J),max({/Symbol d}U/U),max({/Symbol d}{/Symbol r}/{/Symbol r})"<<endl;
	}
    else {
	    tr_char<<"it.,J [au],n_E [au],q [au],max({/Symbol d}J/J),max({/Symbol d}U/U),max({/Symbol d}{/Symbol r}/{/Symbol r})"<<endl;
	}
    if ( timeDependent ) {
	    cout<<"# it.\tt[fs]\tCurr. [Acm^-2]\tnE [au]\tnD [au]\tq [au]\tdU(L/2) [meV]\tmax(U)-min(U)"<<endl;
	}
    else {
	    cout<<setw(10)<<"# it."<<setw(20)<<"Curr. [Acm^-2]"<<setw(15)<<"nE [au]"<<setw(15)<<"nD [au]"<<setw(15)<<"q [au]"<<setw(15)<<"dU(L/2) [meV]"<<setw(15)<<"max(U)-min(U)"<<endl;
	}
    for ( size_t n_it = 0; n_it < n_max; ++n_it ) {

		/// TODO: Change order (BTE -> PE)

		//
		// Solve Poisson EQ
		// If time dependent: alpha = 0
	    p.rho_ = (1.-alpha)*rho_old + alpha*(nD_ - cdX_);
	    p.nE_ = cdX_;
	    p.solve();

		// cout<<"rho_old: "<<rho_old(nx_/2)/AU_cm3<<", rho: "<<p.rho_(nx_/2)/AU_cm3<<", nE: "<<cdX_(nx_/2)/AU_cm3<<endl;

		//
		// Solve Wigner/Boltzmann eq
	    uC_ = p.uNew_;
	    timeDependent ? solveTimeEv() : solveWignerEq();
	    calcCD_X(); // calcCD_K();
	    curr = calcCurr();
	    j0 = j1, j1 = currD_;

		//
		// Maximum values
	    dj = j0 - j1;  // Vectors
	    dRho = rho_old - p.rho_;
	    dJ_x = arma::max(arma::abs(dj/j1));
	    dU_x = arma::max(arma::abs(p.du_/p.uNew_));
		// pFun_x = arma::max(arma::abs(p.pFun_));
	    dRho_x = arma::max(arma::abs(dRho/p.rho_));

	    nc = calcNorm(), nd = calcInt(nD_, dx_), q = nd - nc;
		// u_der = calcFirstDer(p.uNew_, dx_);
		// rho_test = calcSecondDer(p.uNew_, dx_)*epsilonR_/4./M_PI;

	    cout<<setw(10)<<n_it;
	    if ( timeDependent ) {
			// cout.setf( ios::fixed ),
		    cout.precision( 2 );
		    cout<<setw(10)<<std::fixed<<n_it*dt_*AU_s*1e15;
		}
		// cout.setf( ios::scientific );
	    cout.precision( 3 );
	    cout<<std::scientific
			<<setw(20)<<curr*AU_Acm2
			<<setw(15)<<nc
			<<setw(15)<<nd
			<<setw(15)<<q
			<<setw(15)<<p.du_(nx_/2)*AU_eV*1e3
			<<setw(15)<<arma::max(p.uNew_)-arma::min(p.uNew_)<<endl;

	    if (timeDependent) { tr_char<<std::fixed<<n_it*dt_*AU_s*1e15; }
	    else { tr_char<<n_it; }
	    tr_char<<std::scientific<<','<<curr
		<<','<<nc
		<<','<<q
		<<','<<dJ_x
		<<','<<dU_x
		<<','<<dRho_x<<endl;

		//
		// Saving data

		// field<std::string> header(9);
		// arma::mat out_data;
		// out_data.insert_cols(0, x_*AU_nm), header(0) = "x [nm]";
		// out_data.insert_cols(1, p.uNew_*AU_eV), header(1) = "U_{new} [eV]";
		// out_data.insert_cols(2, p.uOld_*AU_eV), header(2) = "U_{old} [eV]";
		// out_data.insert_cols(3, p.du_*AU_eV), header(3) = "{/Symbol d}U [eV]";
		// out_data.insert_cols(4, p.rho_/AU_cm3), header(4) = "{/Symbol r}_{new} [cm^{-3}]";
		// out_data.insert_cols(5, rho_old/AU_cm3), header(5) = "{/Symbol r}_{old} [cm^{-3}]";
		// out_data.insert_cols(6, cdX_/AU_cm3), header(6) = "n_{E} [cm^{-3}]";
		// out_data.insert_cols(7, nD_/AU_cm3), header(7) = "n_{D} [cm^{-3}]";
		// out_data.insert_cols(8, j1*AU_Acm2), header(8) = "J [Acm^{-2}]";
		// // out_data.insert_cols(9, rho_test/AU_cm3), header(9) = "{/Symbol r}_{test} [cm^{-3}]";
		// out_data.save(csv_name("out_data/poisson_test.csv", header));

	    poisson_step<<"## it  j [au]  n_el [au]  n_D [au]  q [au]  max(du)"<<endl;
	    poisson_step<<"# "
			// <<n_it
			<<n_it*dt_*AU_s*1e15
			<<setw(15)<<curr
			<<setw(15)<<nc
			<<setw(15)<<nd
			<<setw(15)<<q
			<<setw(15)<<dU_x<<'\n';
	    poisson_step<<"## it,x [au],rho [au],uNew [au],du [au],u_der [au],j [au],|du/u|\n";
	    for (size_t i=0; i<nx_; ++i)
		    poisson_step
				<<n_it  // col. 1
				<<','<<x_(i)  // col. 2
				<<','<<p.rho_(i)  // col. 3
				<<','<<p.uNew_(i)  // col. 4
				<<','<<p.du_(i)  // col. 5
				<<','<<u_der(i)  // col. 6
				<<','<<j1(i)  // col. 7
				<<','<<p.du_(i)/p.uNew_(i)<<'\n';  // col. 8
	    poisson_step<<"\n";
		// currD_a.resize(n_it), nc_a.resize(n_it), q_a.resize(n_it);
		// dJ_x_a.resize(n_it), dU_x_a.resize(n_it), dRho_x_a.resize(n_it);
		// currD_a(n_it-1) = curr, nc_a(n_it-1) = nc, q_a(n_it-1) = q;
		// dJ_x_a(n_it-1) = dJ_x, dU_x_a(n_it-1) = dU_x, dRho_x_a(n_it-1) = dRho_x;

		//
		// Check convergance
	    conv = dJ_x < max_dJ ? true : false;
	    n_dJ = conv && n_it > n_min ? n_dJ+1 : 0;  // How many times the condition has been met
		//
	    conv = dU_x < max_dU ? true : false;
	    n_dU = conv && n_it > n_min  ? n_dU + 1 : 0;
		//
		// if (n_dU > n_conv) break;  //  && n_dU > n_conv
	    if (n_dU > n_conv) continue;  //  && n_dU > n_conv

		//
		// Values from previous iteration are set to *_old
	    rho_old = p.rho_;
	    p.uOld_ = p.uNew_;
	}
    poisson_step.close();

    std::ofstream test("OutData/poisson_test.csv");
    test << "x [nm],U_{new} [eV],U_{old} [eV],{/Symbol d}U [eV],"
	"{/Symbol r}_{new} [C/cm^{3}],{/Symbol r}_{old} [C/cm^{3}],"
	"n_{E} [cm^{-3}],n_{D} [cm^{-3}],J [Acm^{-2}]"<<endl;
    for (size_t i=0; i<nx_; ++i)
    test<<x_(i)*AU_nm<<','<<p.uNew_(i)*AU_eV<<','<<p.uOld_(i)*AU_eV<<','
	<<p.du_(i)*AU_eV<<','<<p.rho_(i)*E0/AU_cm3<<','<<rho_old(i)*E0/AU_cm3<<','
	<<cdX_(i)/AU_cm3<<','<<nD_(i)/AU_cm3<<','<<j1(i)*AU_Acm2<<endl;
	    test.close();
    saveWignerFun();

    p.uNew_.save("OutData/poisson_pot.bin");
    p.rho_.save("OutData/poisson_rho.bin");

    uStart_ = p.uNew_;
	// std::ofstream pot_out;
	// pot_out.open("OutData/poisson_pot/poisson_pot.out", std::ios::out);
	// pot_out<<"# i u(i)\n";
	// for (size_t i=0; i<nx_; ++i)
	// 	    pot_out<<x_(i)<<' '<<uStart_(i)<<'\n';
	// pot_out.close();

}
