#include <iostream>
#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>
#include <cstdint>

#include "susyqm.hpp"
#include "timing.hpp"

using namespace std;
namespace po = boost::program_options;

double blas_performance_test(uint64 vec_size);

int main(int argc, char **argv)
{
	int      nthreads      = 0;      //If 0, the number of threads will be set automatically
	//XXZ model parameters
	int      M             = 4;       //Maximal value of k+l in state indexing
	double   L             = -1.0;    //Width parameter in the oscillator basis
	uint     parity_sector = 0;       //parity sector for the block H
	int      nev      	   = 0;      //Number of lowest eigenstates ...
	int      ncv      	   = 0;      //ncv parameter in Arnoldi
	double   prec     	   = 1.0E-10; //Target precision for Arnoldi iterations
	//Trotter decomposition parameters
	double   tmin          =  0.0;   //Min evolution time
	double   tmax          = 10.0;   //Max evolution time
	double   beta          = 1.0;    //Inverse temperature
	double   dt            = 0.01;   //Time step for OTOC calculation
	string   datadir       = "./data/";
	bool     susy          = false;
	bool     free_system   = false;
	bool     dont_read_eigensystem  = false;
	bool     dont_write_eigensystem = false;
	bool	 dont_write_eigenvecs   = false; 
	bool     dont_calculate_otocs   = false;
	bool     check_eigensystem      = false;
	bool     blas_test              = false;
	bool     dont_write_otoc_data   = false;
	bool     dont_read_otoc_data    = false;
	bool     print_out_eigenvalues  = false;
	bool     dont_classify_irreps   = false;
	bool     print_irrep_info       = false;
	
	//Option parsing
	po::options_description desc("Simulation of chaotic quantum mechanics with x2*y2 potential");
	desc.add_options()
	     ("help,h", "produce this help message")
	     ("nthreads",	   po::value<int>(       &(nthreads))->default_value(      0  ), "Number of OpenMP threads to use, 0 = automatic choice")
	     ("M",             po::value<int>(              &(M))->default_value(      4  ), "Maximal value of k+l in state indexing"               )
	     ("L",   		   po::value<double>(           &(L))->default_value(   -1.0  ), "Width parameter for the oscillator basis"             )
	     ("parity-sector", po::value<uint>( &(parity_sector))->default_value(      0  ), "Parity sector, 0 for x1-even states, 1, for x1-odd  " )
	     ("nev",           po::value<int>(            &(nev))->default_value(      0  ), "Number of lowest eigenstates"                         )  
		 ("ncv",           po::value<int>(            &(ncv))->default_value(      0  ), "ncv parameter in ARPACK"                              )    
		 ("tmax",          po::value<double>(        &(tmax))->default_value(   10.0  ), "Max. evolution time"								    )
		 ("tmin",          po::value<double>(        &(tmin))->default_value(    0.0  ), "Min. evolution time"								    )
		 ("beta",          po::value<double>(        &(beta))->default_value(    1.0  ), "Inverse temperature"								    )
		 ("dt",            po::value<double>(          &(dt))->default_value(    0.1  ), "Trotter decomposition step"							)
		 ("precision",     po::value<double>(        &(prec))->default_value(1.0E-10  ), "Arnoldi algorithm precision"						    )	
		 ("datadir", 	   po::value<string>(                    &datadir        ), "Directory for data output"                            )
		 ("susy",          po::bool_switch( &susy                                ), "Use supersymmetric Hamiltonian"                       )
		 ("free",          po::bool_switch( &free_system                         ), "For consistency checks - use Hamiltonian without any potential" )
		 ("dont-read-eigsys",       po::bool_switch( &dont_read_eigensystem   ), "Do not attempt to read eigensystem from a file"       )
		 ("dont-write-eigsys",      po::bool_switch( &dont_write_eigensystem  ), "Do not attempt to write eigensystem to a file "       )
		 ("dont-write-eigvecs",     po::bool_switch( &dont_write_eigenvecs    ), "Do not attempt to write eigenvectors to a file (useful to save only eigenvalues if disk space is limited) "    )
		 ("no-otocs",               po::bool_switch( &dont_calculate_otocs    ), "Do not calculate OTOCS (only save eigensystem) "      )
		 ("blas-test",              po::bool_switch( &blas_test               ), "Just test the BLAS performance (RAM bandwidth)"       )
		 ("check-eigensystem",      po::bool_switch( &check_eigensystem       ), "Each time an eigensystem is loaded/calculated, check it for orthonormality/diagonality ")
		 ("dont-write-otoc-data",   po::bool_switch( &dont_write_otoc_data    ), "Don't save operator matrix elements to .tmp files"    )
		 ("dont-read-otoc-data",    po::bool_switch( &dont_read_otoc_data     ), "Don't read operator matrix elements from .tmp files"  )
		 ("print-out-eigenvalues",  po::bool_switch( &print_out_eigenvalues   ), "Print out the eigenvalues that were read/found "      )
		 ("dont-classify-irreps",   po::bool_switch( &dont_classify_irreps    ), "Do not create a file for irrep classification "       )
		 ("print-irrep-info",       po::bool_switch( &print_irrep_info        ), "Print out eigenvalues and their irreps "              );
	     
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	// Check if help is needed
	if(vm.count( "help" )){cout<<desc<<endl; return 1;};
	
	if(parity_sector>1){cerr << ansi::red << "Wrong value! parity_sector = " << parity_sector << ansi::reset << endl << flush; return -1;};
		
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
	
	if(blas_test)
	{
		uint64 ps = (uint64)8*(uint64)1024*(uint64)1024*(uint64)1024;
		double at = blas_performance_test(ps);
		FILE* fr = fopen("./blas_test.dat", "a");
		fprintf(fr, "%03i %2.4E\n", nthreads, openblas_get_num_threads(), at);
		fclose(fr);
		return EXIT_SUCCESS;
	};
	
	cout << ansi::green << "\t\t X1 PARITY IS " << ansi::cyan << (parity_sector==0? "EVEN" : "ODD") << ansi::reset << endl << flush;
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
	cout << ansi::green << "Using " << ansi::magenta << HS->model_name  << ansi::green << " system!" << ansi::reset << flush << endl;

	
	cout << ansi::green << "N  =        " << ansi::magenta << HS->N          << ansi::reset << endl;
	cout << ansi::green << "N2 =        " << ansi::magenta << HS->N2         << ansi::reset << endl;
	cout << ansi::green << "Model name: " << ansi::magenta << HS->model_name << ansi::reset << endl << endl;
		
	HS->CheckSymmetryGroupImplementation();
		
	TIMING_INIT;
	//First trying to read the eigensystem from the file and check whether we have a sufficient number of eigenvectors
	uint fnev = (dont_read_eigensystem? 0 : HS->read_eigensystem(datadir, parity_sector));
	if(nev==0)//If nev was not explicitly specified
	{
		if(fnev==0) //We have not found eigensystem file
		{
			cout << ansi::cyan << "Finding all " << ansi::magenta << HS->N << ansi::cyan << " eigenstates of the Hamiltonian ... " << ansi::reset << endl << flush;
			TIMING_START;
			HS->FindAllEigenstates(parity_sector);
			TIMING_STOP;
			cout << ansi::green << "... Done in " << ansi::magenta << a_time << ansi::green << " sec." << ansi::reset << endl << endl << flush;
			check_eigensystem = true;
			
			if(!dont_write_eigensystem)
			{
				cout << ansi::green << "Saving the eigensystem to the file... " << ansi::reset << flush << endl;  
				HS->write_eigensystem(datadir, parity_sector, !dont_write_eigenvecs);
				cout << ansi::green << "...Done!" << ansi::reset << flush << endl << endl;
			};
		};
		//Otherwise we just go on using fnev eigenstates
	}
	else
	{
		if(fnev<nev)
		{
			cout << ansi::cyan << "Finding the lowest " << ansi::magenta << nev << ansi::cyan << " eigenstates of the Hamiltonian, ncv = " << ncv << " ... " << ansi::reset << endl << flush;
			TIMING_START;
			HS->FindLowestEigenstates(nev, parity_sector, prec, ncv);
			TIMING_STOP;
			cout << ansi::green << "... Done in " << ansi::magenta << a_time << ansi::green << " sec." << ansi::reset << endl << endl << flush;
			check_eigensystem = true;
			
			if(!dont_write_eigensystem)
			{
				cout << ansi::green << "Saving the eigensystem to the file... " << ansi::reset << flush << endl;  
				HS->write_eigensystem(datadir, parity_sector, !dont_write_eigenvecs);
				cout << ansi::green << "...Done!" << ansi::reset << flush << endl << endl;
			};
		}
		else
			HS->nev = nev;
	};
	
	if(check_eigensystem)
	{
        cout << ansi::cyan << "Checking the eigensystem with " << ansi::magenta << HS->nev << ansi::cyan << " eigenstates ... E = " << HS->E << ", psi =" << HS->psi << endl << ansi::reset << flush;
	    TIMING_START;
	    HS->CheckEigensystem(parity_sector);
	    TIMING_STOP;
	    cout << ansi::green << "... Done in " << ansi::magenta << a_time << ansi::green << " sec." << ansi::reset << endl << endl << flush;
    };
    
    if(print_out_eigenvalues)
    {
   		cout << ansi::green << "Eigenvalues in parity sector " << ansi::magenta << parity_sector << ansi::reset << endl << flush;
   		for(uint ie=0; ie<HS->nev; ie++)
   		{
   			cout << std::fixed << std::setprecision(3) << std::showpos << "\t" << ansi::green << ie << ansi::cyan << "\t" << HS->E[nev*parity_sector + ie] << " | " << ansi::yellow; // <<  ansi::reset << endl << flush;
   			for(uint i=0; i<5; i++)
   				cout << "\t" << ansi::yellow << HS->psi[HS->N*nev*parity_sector + ie*HS->N + i];
   			cout << ansi::reset << endl << flush;
   		};
   	};
   	
	if(!dont_classify_irreps)
	{
		cout << ansi::cyan << "Classifying the eigenstates according to irreps of D4d... HS->nev = " << ansi::magenta << HS->nev << endl << ansi::reset << flush;
		TIMING_START;
		HS->ClassifyIrreps(datadir, parity_sector, print_irrep_info);
		TIMING_STOP;
		cout << ansi::green << "... Done in " << ansi::magenta << a_time << ansi::green << " sec." << ansi::reset << endl << endl << flush;
	};
		
	      		
	if(dont_calculate_otocs) return EXIT_SUCCESS;
			
	cout << ansi::cyan << "Initializing the OTOCs... " << endl << ansi::reset << flush;
	TIMING_START;
	bool otoc_data_available = false;
	if(!dont_read_otoc_data)  otoc_data_available = (HS->read_OTOC_data(datadir, parity_sector, nev) > 0);
	if(!otoc_data_available)  HS->init_OTOCS(parity_sector);
	if(!dont_write_otoc_data) HS->write_OTOC_data(datadir, parity_sector);
	TIMING_STOP;
	cout << ansi::green << "... Done in " << ansi::magenta << a_time << ansi::green << " sec." << ansi::reset << endl << endl << flush;
	
	uint Nt = (uint)std::ceil((tmax - tmin)/dt);

	char otocs_fname[512];
	sprintf(otocs_fname, "%s/otoc_%s_M%u_L%1.4lf_p%u_n%u_beta%1.4lf.dat", datadir.c_str(), HS->model_name, HS->M, HS->L, parity_sector, HS->nev, beta);
	FILE* otocs_file = fopen(otocs_fname, "a");
	cout << ansi::cyan << "Calculating OTOCs for " << ansi::magenta << Nt << ansi::cyan << " time steps " << ansi::reset << endl << flush;
	if(otocs_file!=NULL)
		cout << ansi::cyan << "Saving OTOCs to the file " << ansi::magenta << otocs_fname << ansi::reset << endl << endl << flush;
	else
		cerr << ansi::red  << "The file " << otocs_fname << " could not be opened for writing! " << ansi::reset << endl << flush;	
	
	double aA = 0.0, aB = 0.0, Z = 0.0;
	HS->thermal_vevs_regularized(beta, aA, aB, Z);
	
	TIMING_START;
	for(uint it=0; it<Nt; it++)
	{
		double r1, r2; double t = tmin + (double)it*dt;
		HS->OTOC_regularized(beta, t, &r1, &r2);
		fprintf(stdout, "%02.4lf\t|\t%+1.4E\t%+1.4E\t|\t\t%+1.4E\t%+1.4E\t%+1.6E\n", t, r1, r2, aA, aB, Z);
		fflush(stdout);
		if(otocs_file!=NULL){ fprintf(otocs_file, "%02.4lf\t%+1.4E\t%+1.4E\t%+1.4E\t%+1.4E\t%+1.4E\n", t, r1, r2, aA, aB, Z); fflush(otocs_file); };
	};
	TIMING_STOP;
	cout << ansi::green << "... Done, " << ansi::magenta << a_time/(double)Nt << ansi::green << " sec. per step " << ansi::reset << endl << flush;
	if(otocs_file!=NULL) fclose(otocs_file);

	return EXIT_SUCCESS;
}

double blas_performance_test(uint64 vec_size)
{
	double vsize = (double)vec_size*(double)sizeof(t_complex)/(double)(1024*1024*1024);
	cout << ansi::yellow << "JUST TESTING THE SCALING OF SIMPLE BLAS, size of vector is " << ansi::magenta << vsize << " Gb"  << ansi::reset << endl << flush;
	t_complex* v1 = new t_complex[vec_size];
	t_complex* v2 = new t_complex[vec_size];
	
	#pragma omp parallel for
	for(uint i=0; i<vec_size; i++)
	{
		v1[i] = (double)(i+1) + 0.0i;
		v2[i] = 1.0/(double)(i+1) + 0.0i;
	};
	
	TIMING_INIT;
	TIMING_START;
	for(uint iter=0; iter<20; iter++)
	{
		A_pluseq_bB(v1, 0.1 + 0.1i, v2, vec_size);
		cblas_zcopy(vec_size, v2, 1, v1, 1);
	};
	TIMING_STOP;
	cout << ansi::green << "... Done in " << ansi::magenta << a_time << ansi::green << " sec." << ansi::reset << endl << endl << flush;
	
	delete [] v1;
	delete [] v2;
	
	return a_time;
}

/*
double max_err = 0.0;
	for(uint i=0; i<BS.N; i++)
	{
		t_complex* t1 = new t_complex[BS.N];
		t_complex* c1 = new t_complex[BS.N];
		t_complex* c2 = new t_complex[BS.N];
		
		//for(uint j=0; j<BS.N; j++)
		//	t1[j] = (j==i? 1.0 + 0.0i : 0.0 + 0.0i);
		rand_vec(t1, BS.N);
		
		t_complex* HM = BS.H(); //Real-valued N x N matrix
	
		psi_eq_A_mult_chi(c1, HM, t1, BS.N);
		BS.H(t1, c2);
		
		max_err = std::max(max_err, norm_diff(c1, c2, BS.N));
		
		delete [] t1; delete [] c1; delete [] c2; delete [] HM;
	};
	cout << "Sparse vs Dense err: " << max_err << endl;
*/



/*	double max_err = 0.0;
   	for(uint ie=0; ie<HS->nev; ie++)
   	{
   		double* epsi = new double[2*N];
   		std::copy(HS->psi  + N*ie, HS->psi + N*ie + N, epsi);
   		std::fill(epsi + N, epsi + 2*N, 0.0);
   		
   		HS->H(epsi, v1); rescale(v1, 1.0/HS->E[ie], 2*N);
   		double err = norm_diff(epsi, v1, 2*N);
   		max_err = std::max(max_err, err);
   		
   		double max_gerr = 0.0;
   		for(uint ig=0; ig<16; ig++)
   		{
   			HS->SymmetryTransformation(ig, epsi, v2);
   			HS->H(v2, v1); rescale(v1, 1.0/HS->E[ie], 2*N);
   			double gerr = norm_diff(v2, v1, 2*N);
   			max_gerr = std::max(max_gerr, gerr);
   		};
   		
   		cout << "Eigenstate " << ie << ", E = " << HS->E[ie] << ", err = " << err << ", gerr = " << max_gerr << endl << flush;
   		
	   	//HS->ClassifyIrreps(HS->psi  + HS->N*ie, parity_sector);
	}; */
