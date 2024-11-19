#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#include <iterator>
#include <algorithm>
#include "cminpack.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
int Chi2( void *p, int m, int n, const double *pars, double *fvec, int iflag);

int nDataPoints;
int num_pars = -1;
int grad_order = -1;

struct DataPoint{
	double x;
	double y;
	double z;
	double B;
	double err;
};
std::vector<DataPoint> data;

// Main
int main(int argc, char ** argv){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [grad order]\n";
		cerr << "**********************************\n";
		cerr << "\t\t[grad order==0]: G00\n";
		cerr << "\t\t[grad order==1]: G00, G1m1, G10, G11\n";
		cerr << "\t\t[grad order==2]: G00, G1m1, G10, G11, G2m2, ...\n";
		cerr << "**********************************\n";
		return -1;
	}

	

	///////////////////////////////////////////////////////////////////////////
	// Load the data structure that we will use:
	std::ifstream f;
	std::string line;
	f.open("../test_data.txt");
	f.ignore(1000, '\n'); // skip 1 line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			double x, y, z, b, err;
			int i;
			ss >> i >> x >> y >> z >> b >> err;
			DataPoint dat;
			dat.x = x;
			dat.y = y;
			dat.z = z;
			dat.B = b;
			dat.err = err;
			data.push_back( dat );
		}
	}
	f.close();
	nDataPoints = data.size();

	///////////////////////////////////////////////////////////////////////////
	// Now we can setup our Chi2 function and minimize it!
	grad_order = atoi(argv[1]);
	int nParams = -1;
	if( grad_order == 0 ) nParams = 1;
	else if( grad_order == 1 ) nParams = 4;
	else if( grad_order == 2 ) nParams = 9;
	else{ cerr << "unexpected input!\n"; exit(-1); }


	double *params = new double[nParams];
	// Set initial guess for params
	if( grad_order == 0 ){
		params[0] = 0.94195;
	}
	else if( grad_order == 1 ){
		params[0] = 0.94195;
		params[1] = 0;
		params[2] = 1;
		params[3] = 0;
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
	cout << "...minimizing...\n";
	int r	=	lmdif1(Chi2,	NULL,
			nDataPoints,	nParams,	params,
			fvec,	1e-7,
			iwork,	dwork,	sdwork);
	cout << "...done minimizing\n";
	cout << "r value of minimization: " << r << "\n";
	// Print out the final parameters found:
	for( int i = 0 ; i < nParams ; ++i ){
		cout << "Par " << i << ":\t" << params[i] << "\n";
	}



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

		if( grad_order == 0 ){
			theo = pars[0];
		}
		if( grad_order == 1 ){
			theo = pars[0] + pars[1]/1E6*data.at(i).x + pars[2]/1E6*data.at(i).z + pars[3]/1E6*data.at(i).y;
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
