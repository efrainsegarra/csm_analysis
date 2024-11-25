#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#include <iterator>
#include <algorithm>
#include "cminpack.h"
#include <stdlib.h>
#include <iomanip>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
int Chi2( void *p, int m, int n, const double *pars, double *fvec, int iflag);

int nDataPoints;
int num_pars = -1;
int grad_order = -1;
int bootstrap_opt = -1;
int rand_offset = -1;
struct DataPoint{
	double x;
	double y;
	double z;
	double B;
	double err;
};
std::vector<DataPoint> data;
std::vector<DataPoint> bootstrap_data;

// Main
int main(int argc, char ** argv){

	if (argc<4){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [grad order] [bootstrap opt] [rand offset]\n";
		cerr << "**********************************\n";
		cerr << "\t\t[grad order==0]: G00\n";
		cerr << "\t\t[grad order==1]: G00, G1m1, G10, G11\n";
		cerr << "\t\t[grad order==2]: G00, G1m1, G10, G11, G2m2, ...\n";
		cerr << "**********************************\n";
		cerr << "\t\t[bootstrap opt == 0]: Don't bootstrap\n";
		cerr << "\t\t[bootstrap opt == 1]: Bootstrap\n";
		cerr << "**********************************\n";
		cerr << "\t\t[rand offset == 0]: No offset\n";
		cerr << "\t\t[rand offset == 1]: Randomize an offset\n";
		cerr << "**********************************\n";
		return -1;
	}

	// Random seed for offsets
	srand (time(NULL));
	//srand(1);	

	///////////////////////////////////////////////////////////////////////////
	// Load the data structure that we will use:
	rand_offset = atoi(argv[3]);
	if( rand_offset != 0 && rand_offset != 1 ){ cerr << "unexpected input!\n"; exit(-1); }
	std::ifstream f;
	std::string line;
	f.open("../test_data.txt");
	f.ignore(1000, '\n'); // skip 1 line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			double x, y, z, b, err;
			int i;
			double offset = 0;
			if( rand_offset){
				offset = rand() % 20 - 10; // in pT
				offset *= 1E-6;
			}
			ss >> i >> x >> y >> z >> b >> err;
			DataPoint dat;
			dat.x = x;
			dat.y = y;
			dat.z = z;
			dat.B = b + offset;
			dat.err = err;
			data.push_back( dat );
		}
	}
	f.close();
	nDataPoints = data.size();

	bootstrap_opt = atoi(argv[2]);
	if( bootstrap_opt != 0 && bootstrap_opt != 1 ){ cerr << "unexpected input!\n"; exit(-1); }
	if( bootstrap_opt){
		/////////////////////////////////////////////////////////////////////////////
		// Create RANDOM samples of data sets for boostrapping error estimation
		std::random_device rd; // obtain a random number from hardware
		std::mt19937 gen(rd()); // seed the generator

		std::uniform_int_distribution<> distr_data(0, data.size()-1);
		for( int i = 0 ; i < data.size() ; ++i ){ 	
			int rand_elem = distr_data(gen);
			bootstrap_data.push_back( data.at(rand_elem) );
		}
		data.clear();
		data = bootstrap_data;
	}

	///////////////////////////////////////////////////////////////////////////
	// Now we can setup our Chi2 function and minimize it!
	grad_order = atoi(argv[1]);
	int nParams = -1;
	if( grad_order == 0 ) nParams = 1;
	else if( grad_order == 1 ) nParams = 4;
	else if( grad_order == 2 ) nParams = 5;
	else if( grad_order == 3 ) nParams = 8;
	else{ cerr << "unexpected input!\n"; exit(-1); }


	double *params = new double[nParams];
	// Set initial guess for params
	if( grad_order == 0 ){
		params[0] = 0.94195;	// G00 (muT)
	}
	else if( grad_order == 1 ){
		params[0] = 0.94195;	// G00 (muT)
		params[1] = 0;		// G11 
		params[2] = 1;		// G10
		params[3] = 0;		// G1m1
	}
	else if( grad_order == 2 ){
		params[0] = 0.94195;	// G00
		params[1] = 0;		// G11
		params[2] = 1;		// G10
		params[3] = 0;		// G1m1
		params[4] = 0;		// G20
	}
	else if( grad_order == 3 ){
		params[0] = 0.94195;	// G00
		params[1] = 1;		// G10
		params[2] = 0;		// G1m1
		params[3] = 0;		// G20
		params[4] = 0;		// G2m2
		params[5] = 0;		// G2m1
		params[6] = 0;		// G21
		params[7] = 0;		// G22
	}

	// Define some working memory for the minimization algorithm
	nDataPoints+=1;
	double 	*fvec	=	new double[nDataPoints];
	int 	*iwork	=	new int[nDataPoints];
	int 	sdwork	=	nDataPoints*nParams + 5*nParams+nDataPoints;
	double 	*dwork	=	new double[sdwork];

	///////////////////////////////////////////////////////////////////////////
	// Minimize!
	num_pars = nParams;
	//cout << "...minimizing...\n";
	int r	=	lmdif1(Chi2,	NULL,
			nDataPoints,	nParams,	params,
			fvec,	1e-7,
			iwork,	dwork,	sdwork);
	//cout << "...done minimizing\n";
	//cout << "r value of minimization: " << r << "\n";
	// Print out the final parameters found:
	for( int i = 0 ; i < nParams ; ++i ){
		//cout << "Par " << i << ":\t" << params[i] << "\n";
		cout << std::setprecision(15); 
		cout << params[i] << "\t";
	}
	cout << "\n";



	///////////////////////////////////////////////////////////////////////////
	// Cleanup
	delete []params;
	delete []fvec;
	delete []iwork;
	delete []dwork;

	return 1;
}

int Chi2( void *p, int m, int n, const double *pars, double *fvec, int iflag){
	// p not used??
	// m = num data points
	// n = 0, not used??
	// x = parameters
	// fvec = residuals
	// iflag not used??

	// Want to return chi2 which is calculated given the parameters to be minimized
	std::cerr << "************ Fitting in progress *************\n";
	int dat_point = 0;

	///////////////////////////////////////////////////////////////////////////
	// Perform loop over exp data
	double chi2 = 0.;
	for( int i = 0 ; i < data.size(); ++i ){
		double theo = 0;

		double xd = data.at(i).x;
		double yd = data.at(i).y;
		double zd = data.at(i).z;
			
		if( grad_order == 0 ){ // just G00 (muT)
			theo = pars[0];
		}
		if( grad_order == 1 ){ // G00 (muT), G11 & G10 & G1m1 in pT/cm
			theo = pars[0] + pars[1]/1E6*xd+ pars[2]/1E6*zd + pars[3]/1E6*yd;
		}
		if( grad_order == 2 ){ // G00 (muT), G11 & G10 & G1m1 & G20 in pT/cm
			theo = pars[0] + pars[1]/1E6*xd+ pars[2]/1E6*zd + pars[3]/1E6*yd + pars[4]/1E6*(zd*zd - 0.5*(xd*xd + yd*yd));
		}
		if( grad_order == 3 ){ // G00 (muT), G11 & G10 & G1m1 & G20 & G2m2 & G2m1 & G21 & G22
			theo = pars[0]  // G00
				+ pars[1]/1E6*xd+ pars[2]/1E6*zd + pars[3]/1E6*yd 
				+ pars[4]/1E6*(zd*zd - 0.5*(xd*xd + yd*yd)) +  pars[5]/1E6*(2*yd*zd) + pars[6]/1E6*(2*xd*zd) + pars[7]/1E6*(xd*xd-yd*yd);
		}

		std::cerr << "\tData-pt: " << i << " " << data.at(i).x << " " << data.at(i).y << " " << data.at(i).z << " " << data.at(i).B << " " << theo << "\n";

		double this_chi2 = (data.at(i).B - theo) / data.at(i).err;
		fvec[dat_point] = this_chi2;
		++dat_point;
		chi2 += pow(this_chi2,2);
	} // end loop over data


	fvec[dat_point] = 0;
	++dat_point;

	std::cerr << "------------Finished calculations!------------\n";
	std::cerr << "\tOptions used:\n";
	std::cerr << "\t\tgrad_order = " << grad_order << "\n";
	std::cerr << "\tCurrent parameters: " << num_pars << "\n";
	for( int i = 0 ; i < num_pars ; ++i ) std::cerr << "\t\t" << pars[i] << "\n";
	std::cerr << "\tCurrent chi2: " << chi2 << "\n";
	std::cerr << "\tReduced chi2: " << chi2 / (nDataPoints-num_pars) << "\n";
	std::cerr << "\tCurrent residuals: \n";
	for(int i = 0; i < m ; ++i) std::cerr << "\t\t" << fvec[i] << "\n";
	std::cerr << "**********************************************\n\n";

	return 0.;
}
