#include <iostream>
#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>
#include <cstdint>

#include "susyqm.hpp"
#include "timing.hpp"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
	TIMING_INIT;
	/* This is the lightweight code that will read operator matrix elements (from *.tmp files) 
	 and eigenvalues of H (from evals*.dat files) and only calculate and save the OTOCs.
	 This is meant to work in the big M - small nev regime that appears more advantageous for analysis
	 of OTOCs */
	
	int      nthreads = 0;      //If 0, the number of threads will be set automatically
	//XXZ model parameters
	int      M        = 4;  //Maximal value of k+l in state indexing
	double   L        = -1.0;    //Width parameter in the oscillator basis
	int      nev      = 0;      //Number of lowest eigenstates ...
	//Trotter decomposition parameters
	double   tmin     =  0.0;   //Min evolution time
	double   tmax     = 10.0;   //Max evolution time
	double   beta     = 1.0;    //Inverse temperature
	double   dt       = 0.01;   //Time step for OTOC calculation
	string   datadir  = "./data/";
	bool     susy     = false;
	bool     free_system  = false;
	bool     vevs_only    = false;
	bool     rewrite_data = false; //If true, overwrite output rather than append 

	//Option parsing
	po::options_description desc("Simulation of chaotic quantum mechanics with x2*y2 potential");
	desc.add_options()
	     ("help,h", "produce this help message")
	     ("nthreads",	   po::value<int>(       &(nthreads))->default_value(      0  ), "Number of OpenMP threads to use, 0 = automatic choice")
	     ("M",             po::value<int>(              &(M))->default_value(      4  ), "Maximal value of k+l in state indexing"               )
	     ("L",   		   po::value<double>(           &(L))->default_value(   -1.0  ), "Width parameter for the oscillator basis"             )
	     ("nev",           po::value<int>(            &(nev))->default_value(      0  ), "Number of lowest eigenstates"                         )  
		 ("tmax",          po::value<double>(        &(tmax))->default_value(   10.0  ), "Max. evolution time"								    )
		 ("tmin",          po::value<double>(        &(tmin))->default_value(    0.0  ), "Min. evolution time"								    )
		 ("beta",          po::value<double>(        &(beta))->default_value(    1.0  ), "Inverse temperature"								    )
		 ("dt",            po::value<double>(          &(dt))->default_value(    0.1  ), "Trotter decomposition step"							)
		 ("datadir", 	   po::value<string>(      &datadir                           ), "Directory for data output"                            )
		 ("susy",          po::bool_switch( &susy                                     ), "Use supersymmetric Hamiltonian"                       )
		 ("free",          po::bool_switch( &free_system                              ), "For consistency checks - use Hamiltonian without any potential" )
		 ("vevs-only",     po::bool_switch( &vevs_only                                ), "Only print out the thermal vevs and stop, do not calculate/write OTOCs")
		 ("rewrite-data",  po::bool_switch( &rewrite_data                             ), "Rewrite data files on each run");
		 
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	// Check if help is needed
	if(vm.count( "help" )){cout<<desc<<endl; return 1;};
		
	//Automatically setting the value of L
	bool L_auto_set = false;
	if(L<0.0)
	{
		double A = 0.961624; double B = 0.0409431; double C = 0.377368;
		L = (susy? A + B*pow((double)(2*M-1), C) : 1.12246);
		L_auto_set = true;
	};
	
	nthreads = (nthreads==0? omp_get_max_threads() : nthreads);
	omp_set_num_threads(nthreads);
	openblas_set_num_threads(nthreads);
	
	cout << ansi::white << "Using " << ansi::magenta << nthreads << ansi::white << " OpenMP threads, openblas_num_threads = " << ansi::magenta << openblas_get_num_threads() << ansi::reset << endl;
	
	cout << ansi::green << "M  =   " << ansi::magenta << M    << ansi::reset << endl; //", total number of states is " << N << endl;
	cout << ansi::green << "L  =   " << ansi::magenta << L    << ansi::yellow << (L_auto_set? "[Auto set]" : "") << ansi::reset << endl;
	cout << ansi::green << "tmin = " << ansi::magenta << tmin << ansi::reset << endl;
	cout << ansi::green << "tmax = " << ansi::magenta << tmax << ansi::reset << endl;
	cout << ansi::green << "beta = " << ansi::magenta << beta << ansi::reset << endl;
	cout << ansi::green << "dt =   " << ansi::magenta << dt   << ansi::reset << endl;
	
	HamiltonianSystem* HS = NULL;
	if(susy)
	{
		SupersymmetricSystem SS(M, L);
		HS = &SS;
	}
	else
	{
		BosonicSystem BS(M, L, free_system);
		HS = &BS;
	};
	cout << ansi::green << "Using " << ansi::magenta << HS->model_name << ansi::green << " system!" << ansi::reset << flush << endl;
	
	cout << ansi::green << "N  =        " << ansi::magenta << HS->N          << ansi::reset << endl;
	cout << ansi::green << "N2 =        " << ansi::magenta << HS->N2         << ansi::reset << endl;
	cout << ansi::green << "Model name: " << ansi::magenta << HS->model_name << ansi::reset << endl << endl;
	
	uint Nt = (uint)std::ceil((tmax - tmin)/dt);
	double* pr1   = new double[2*Nt]; std::fill(pr1, pr1 + 2*Nt, 0.0);
	double* pr2   = new double[2*Nt]; std::fill(pr2, pr2 + 2*Nt, 0.0);
	double  pZ[2]  = {0.0, 0.0};
	double  pA[2]  = {0.0, 0.0};
	double  pB[2]  = {0.0, 0.0};
	double  pZr[2] = {0.0, 0.0};
	double  pAr[2] = {0.0, 0.0};
	double  pBr[2] = {0.0, 0.0};
	
	for(uint parity_sector=0; parity_sector<2; parity_sector++)
	{
		cout << ansi::blue << "\t\t X1 PARITY IS " << ansi::cyan << (parity_sector==0? "EVEN" : "ODD") << ansi::reset << endl << flush;
	
		//First trying to read the eigensystem from the file and check whether we have a sufficient number of eigenvectors
		uint fnev = HS->read_eigensystem(datadir, parity_sector, false);
		if(nev > fnev)
		{
			cerr << ansi::red << "Only data files with " << ansi::yellow << fnev << ansi::red << " eigenvalues found in the folder " << ansi::yellow << datadir << ansi::red << ", but " << ansi::yellow << nev << ansi::red << " required!!!"  << ansi::reset << flush << endl; 
			continue;
		};
	
		if(HS->read_OTOC_data(datadir, parity_sector, nev)<=0)
		{
			cerr << ansi::red << "Failed to read in the OTOC data!!! " << ansi::reset << flush << endl; 
			continue;
		};
	
		cout << ansi::red << "nev = " << nev << ansi::reset << endl << flush;
		HS->resize_OTOC_data(nev);
		cout << ansi::green << "Number of eigenstates adjusted to " << ansi::magenta << HS->nev  << ansi::reset << endl << flush;
	
		//HS->A = HS->B; //SUPERKASTYL

		cout << endl << ansi::green << "\t Calculating thermal VEVs of both operators..." << ansi::reset << endl << flush;
		TIMING_START;
		HS->thermal_vevs(            beta, pA[parity_sector],  pB[parity_sector],  pZ[parity_sector]);
		HS->thermal_vevs_regularized(beta, pAr[parity_sector], pBr[parity_sector], pZr[parity_sector]);
		TIMING_STOP;
		cout << ansi::green << "... Done, " << ansi::magenta << a_time/(double)Nt << ansi::green << " sec. per step " << ansi::reset << endl << flush;

		if(vevs_only)
		{
			cout << endl;
			cout << ansi::green << "\t\t pA  = " << ansi::magenta << pA[parity_sector]  << ansi::green << "\t, pB  = " << ansi::magenta << pB[parity_sector]  << flush << ansi::reset << endl;
			cout << ansi::green << "\t\t pAr = " << ansi::magenta << pAr[parity_sector] << ansi::green << "\t, pBr = " << ansi::magenta << pBr[parity_sector] << flush << ansi::reset << endl;
			cout << endl;
			continue;
		};
		
		cout << endl << ansi::green << "\t Calculating the OTOCs for " << ansi::magenta << Nt << ansi::green << " time steps..." << ansi::reset << endl << flush;
		TIMING_START;
		for(uint it=0; it<Nt; it++)
		{
			HS->OTOC_regularized(beta, tmin + (double)it*dt, pr1 + Nt*parity_sector + it, pr2 + Nt*parity_sector + it);
			if(it%1==0) //(it%(Nt/10))
			{
				fprintf(stdout, "%3.4lf\t|\t%+1.4E\t%+1.4E\t|\t\t%+1.4E\t%+1.4E\t%+1.4E\n", tmin + (double)it*dt, pr1[Nt*parity_sector + it], pr2[Nt*parity_sector + it], pA[parity_sector], pB[parity_sector], pZ[parity_sector]);
				fflush(stdout);
			};
		};
		TIMING_STOP;
		cout << ansi::green << "... Done, " << ansi::magenta << a_time/(double)Nt << ansi::green << " sec. per step " << ansi::reset << endl << flush;
	};
		
	double w0 = pZ[0]/(pZ[0] + pZ[1]);
	double w1 = pZ[1]/(pZ[0] + pZ[1]);
	double aA = w0*pA[0] + w1*pA[1];
	double aB = w0*pB[0] + w1*pB[1];
	
	double w0r = pZr[0]/(pZr[0] + pZr[1]);
	double w1r = pZr[1]/(pZr[0] + pZr[1]);
	double aAr = w0r*pAr[0] + w1r*pAr[1];
	double aBr = w0r*pBr[0] + w1r*pBr[1];
	
	cout << endl << ansi::green << "Iterations over parity sectors finished, " << ansi::magenta << " w0 = " << w0 << ", w1 = " << w1 << ansi::reset << endl << flush;
	
	if(vevs_only)
	{
		cout << endl;
		cout << ansi::green << "\t\t aA  = " << ansi::magenta << aA  << ansi::green << "\t, aB  = " << ansi::magenta << aB  << flush << ansi::reset << endl;
		cout << ansi::green << "\t\t aAr = " << ansi::magenta << aAr << ansi::green << "\t, aBr = " << ansi::magenta << aBr << flush << ansi::reset << endl;
		cout << endl;
		return EXIT_SUCCESS;
	};	
	
	double* r1    = new double[Nt]; double* r2    = new double[Nt]; 
	double* otoc  = new double[Nt]; double* dotoc = new double[Nt];
	
	for(uint it=0; it<Nt; it++)
	{
		r1[it] = w0*pr1[Nt*0 + it] + w1*pr1[Nt*1 + it];
		r2[it] = w0*pr2[Nt*0 + it] + w1*pr2[Nt*1 + it];
		otoc[it] = (((r2[it] - r1[it])>0.0)? log(2.0*(r2[it] - r1[it])) : -20.0);
	};
	
	char otocs_fname[512];
	sprintf(otocs_fname, "%s/otoc_to_plot_%s_M%u_L%1.4lf_n%u_beta%1.4lf.dat", datadir.c_str(), HS->model_name, HS->M, HS->L, HS->nev, beta);
	FILE* otocs_file = (rewrite_data? fopen(otocs_fname, "w") : fopen(otocs_fname, "a"));
	if(otocs_file!=NULL)
		cout << ansi::cyan << "Saving OTOCs to the file " << ansi::magenta << otocs_fname << ansi::reset << endl << endl << flush;
	else
		cerr << ansi::red  << "The file " << otocs_fname << " could not be opened for writing! " << ansi::reset << endl << flush;	
		
	dotoc[0] = (-otoc[2] + 4*otoc[1] - 3*otoc[0])/(2*dt); //Second order forward derivative discretization
	for(uint it=1; it<Nt-1; it++)
		dotoc[it] = (otoc[it+1] - otoc[it-1])/(2*dt);
	dotoc[Nt-1] = (otoc[Nt-3] - 4*otoc[Nt-2] + 3*otoc[Nt-1])/(2*dt); //Second order backward derivative discretization	
	
	if(otocs_file!=NULL)
	{
		for(uint it=0; it<Nt; it++)
			fprintf(otocs_file, "%3.6lf\t%+1.8E\t%+1.8E\t%+1.8E\t%+1.8E\t%+1.8E\t%+1.8E\t%+1.8E\t%+1.8E\t%+1.8E\t%+1.8E\t%+1.8E\t%+1.8E\n", tmin + (double)it*dt, r1[it], r2[it], otoc[it], dotoc[it], aA*aB, aA, aB, w0, w1, aAr*aBr, aAr, aBr);
		fflush(otocs_file);
		fclose(otocs_file);
	};

	return EXIT_SUCCESS;
}

