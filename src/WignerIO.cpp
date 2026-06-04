#include "lib.hpp"
#include "WignerSolver.hpp"


// -------------------------
// Read parameters from file
// -------------------------

/**
 * Reads single parameter value
 * @param str input file line
 * @param val_name parameter to search
 * @param ln input file line number
 * @return parameter value
 */
template<typename T>
T ReadPar(std::string str, std::string val_name, std::size_t ln){
    size_t pp, eol;
    T val;
    eol = str.find(';');
    if(eol==std::string::npos) eol = str.size()-1;
    pp = str.substr(0, eol).find(val_name);
    if (pp!=std::string::npos){
		    val = std::stod( str.substr(pp+val_name.size(), eol) );
		    std::cout<<"found "+val_name+" at ("<<ln<<','<<pp<<") equal to "<<val<<std::endl;
			// pp = str.find(val_name, pp+1);
	}
    else val = -1;
    return val;
}


/**
 * Reads parameters from a file
 * @deprecated parameters are set through main function 
 */
std::map<std::string, double> readParameters(std::string filename){
    std::ifstream file (filename);
    std::map<std::string, double> val = {
		{"calc_mode", 0},
		{"contact_temp", 77},
		{"fermi_energy", .00316},  // 0.086 eV
		{"dop_con_left", 0},
		{"dop_con_right", 0},
		{"effective_mass", 0.067},
		{"device_lenght", 3780.72},  // 200 nm
		{"contact_lenght", 0},
		{"max_k", -1},
		{"xspace_step_nr", 100},
		{"kspace_step_nr", 100},
		{"max_voltage", .0147},  // approx. 0.4 V
		{"voltage_step_nr", 40},
		{"courant_num", 1},
		{"gwp_x0", -1},
		{"gwp_dx", -1},
		{"gwp_p0", -1},
		{"gwp_dp", -1},
		{"pot_type", 0},
		{"part_num",1000},
		{"inelastic_sc", 0},
		{"elastic_sc", 0},
		{"spatial_decoh", 0},
		{"contact_diss", 0},
		{"contact_dist", 0},
		{"dconc_left", 0},
		{"dconc_right", 0},
		{"v_bias", 0}
	};
    double tmp;
    if (file.is_open()){
	    std::string line;
	    for (std::pair<std::string, double> ival : val){
		    std::size_t ln=0;
		    while ( getline (file,line) ){
			    ln++;
			    if (line[0] != '#'){
					// std::cout << line << ' ' << line[0] << std::endl;
				    tmp = ReadPar<double>(line, ival.first, ln);
				    if (tmp != -1)
						    val[ival.first] = tmp;
				}
			}
		    file.clear();
		    file.seekg (0, std::ios::beg);
		    if (val[ival.first] == -1)
				    std::cout<<ival.first+" NOT FOUND, set to default value: "<<val[ival.first]<<std::endl;
		}
	    file.close();
	}

    if (val["max_k"]==-1){
	    val["max_k"] = M_PI/2./(val["device_lenght"]/val["xspace_step_nr"]);
	    std::cout<<"max_k NOT FOUND or EQUAL TO -1, is set to "<<val["max_k"]<<" nm^-1"<<std::endl;
	}
    return val;
}


// -----------------------------
// Print and save data to a file
// -----------------------------


/**
 * Save matrix in a format readable by gnuplot
 */
void saveMatGP(arma::mat data, std::string filename) {
    std::ofstream file(filename);
    file<<"# x y z\n";
    for (size_t i=0; i<data.n_rows; ++i){
	    for (size_t j=0; j<data.n_cols; ++j)
		    file<<i<<' '<<j<<' '<<data(i,j)<<'\n';
	    file<<"\n";
	}
    file.close();
}


/**
 * Saves Distribution function to wf.out, wf.bin, and wf.z (GLE format) files
 */
void WignerSolver::saveDistFun() {
    // --------------------------------
    // Gnuplot format (3 columns x y z)
    // --------------------------------
    std::ofstream wf_out("output/wf.out");
    wf_out<<"# x [nm] k [a.u.] f [a.u.]\n";
    for (size_t i=0; i<nx_; ++i){
	    for (size_t j=0; j<nk_; ++j)
		    wf_out<<x_(i)<<' '<<k_(j)<<' '<<f_(i,j)<<'\n';
	    wf_out<<"\n";
	}
    wf_out.close();
    // -----------
    // Binary file
    // -----------
    f_.save("output/wf.bin");
    // -------------------
    // GLE *.z file format
    // -------------------
    wf_out.open("output/wf.z", std::ios::out);
    wf_out<<"! nx "<<nx_<<" ny "<<nk_<<" xmin "<<0<<" xmax "<<l_*AU::nm<<" ymin "<<-kmax_<<" ymax "<<kmax_<<'\n';
    for (size_t j=0; j<nk_; ++j){
	    for (size_t i=0; i<nx_; ++i)
		    wf_out<<f_(i,j)<<' ';
	    wf_out<<'\n';
	}
    wf_out.close();
}


/// Print parameters to console
void WignerSolver::printParam()
{
    int cw_n = 25, cw_v = 25;

    std::cout << std::left;

    std::cout<<std::endl;
    std::cout.fill('=');
    std::cout.width(cw_n+2*cw_v);
    std::cout<<'#'<<'#'<<std::endl;

    std::cout.fill(' ');

    std::cout.width(cw_n+2*cw_v); std::cout<<"# SET PARAMETERS"<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# Variable name";
    std::cout.width(cw_v); std::cout<<"variable's value (a.u.)";
    std::cout.width(cw_v); std::cout<<"variable's value (SI)"<<'#'<<std::endl;
	// ////////// Effective mass //////////
    std::cout.width(cw_n); std::cout<<"# m*";
    std::cout.width(cw_v); std::cout<<m_;
    std::cout.width(cw_v); std::cout<<'-'<<'#'<<std::endl;
	// ////////// Temperature //////////
    std::cout.width(cw_n); std::cout<<"# contact_temp_";
    std::cout.width(cw_v); std::cout<<temp_;
    std::cout.width(cw_v); std::cout<<'-'<<'#'<<std::endl;
	// ////////// Lenght //////////
    std::cout.width(cw_n); std::cout<<"# L";
    std::cout.width(cw_v); std::cout<<l_;
    std::cout.width(cw_v); std::cout<<l_*AU::nm<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# kmax";
    std::cout.width(cw_v); std::cout<<kmax_;
    std::cout.width(cw_v); std::cout<<kmax_/AU::nm<<'#'<<std::endl;
	// ////////// Numerical grid parameters //////////
    std::cout.width(cw_n); std::cout<<"# Nx";
    std::cout.width(cw_v); std::cout<<nx_;
    std::cout.width(cw_v); std::cout<<'-'<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# dx";
    std::cout.width(cw_v); std::cout<<dx_;
    std::cout.width(cw_v); std::cout<<dx_*AU::nm<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# Nk";
    std::cout.width(cw_v); std::cout<<nk_;
    std::cout.width(cw_v); std::cout<<'-'<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# dk";
    std::cout.width(cw_v); std::cout<<dk_;
    std::cout.width(cw_v); std::cout<<dk_/AU::nm<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# dx*dk";
    std::cout.width(cw_v); std::cout<<dx_*dk_;
    std::cout.width(cw_v); std::cout<<'-'<<'#'<<std::endl;
	// ////////// Dissipation //////////
    std::cout.width(cw_n); std::cout<<"# scR";
    std::cout.width(cw_v); std::cout<<scR_;
    std::cout.width(cw_v); std::cout<<scR_/AU::s<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# scM";
    std::cout.width(cw_v); std::cout<<scM_;
    std::cout.width(cw_v); std::cout<<scM_/AU::s<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# lambda";
    std::cout.width(cw_v); std::cout<<lambda_;
    std::cout.width(cw_v); std::cout<<lambda_*AU::nm*AU::nm*AU::s<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# scG";
    std::cout.width(cw_v); std::cout<<scG_;
    std::cout.width(cw_v); std::cout<<scG_/AU::s<<'#'<<std::endl;
	// ////////// Boundary condition //////////
    std::cout.width(cw_n); std::cout<<"# fermi_energy (left)";
    std::cout.width(cw_v); std::cout<<uL_;
    std::cout.width(cw_v); std::cout<<uL_*AU::eV<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# fermi_energy (right)";
    std::cout.width(cw_v); std::cout<<uR_;
    std::cout.width(cw_v); std::cout<<uR_*AU::eV<<'#'<<std::endl;

    std::cout.fill('=');
    std::cout.width(cw_n+2*cw_v);
    std::cout<<'#'<<'#'<<std::endl;
    std::cout.fill(' ');
    std::cout<<std::endl;

    std::cout.fill('=');
    std::cout.width(cw_n+2*cw_v);
    std::cout<<'#'<<'#'<<std::endl;
    std::cout.fill(' ');

    std::cout.width(cw_n+2*cw_v); std::cout<<"# SCALING PARAMETERS"<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# Variable name";
    std::cout.width(cw_v); std::cout<<"variable's value (a.u.)";
    std::cout.width(cw_v); std::cout<<"variable's value (SI)"<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# AU::nm";
    std::cout.width(cw_v); std::cout<<1;
    std::cout.width(cw_v); std::cout<<AU::nm<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# AU::eV";
    std::cout.width(cw_v); std::cout<<1;
    std::cout.width(cw_v); std::cout<<AU::eV<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# AU::s";
    std::cout.width(cw_v); std::cout<<1;
    std::cout.width(cw_v); std::cout<<AU::s<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# AU::A";
    std::cout.width(cw_v); std::cout<<1;
    std::cout.width(cw_v); std::cout<<AU::A<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# AU::cm";
    std::cout.width(cw_v); std::cout<<1;
    std::cout.width(cw_v); std::cout<<AU::cm<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# AU::cm^2";
    std::cout.width(cw_v); std::cout<<1;
    std::cout.width(cw_v); std::cout<<AU::cm2<<'#'<<std::endl;
    std::cout.width(cw_n); std::cout<<"# AU::cm^3";
    std::cout.width(cw_v); std::cout<<1;
    std::cout.width(cw_v); std::cout<<AU::cm3<<'#'<<std::endl;

    std::cout.fill('=');
    std::cout.width(cw_n+2*cw_v);
    std::cout<<'#'<<'#'<<std::endl;
    std::cout.fill(' ');
    std::cout<<std::endl;
}
