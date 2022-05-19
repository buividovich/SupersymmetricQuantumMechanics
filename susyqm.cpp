#include "susyqm.hpp"

uint64       HamiltonianSystem::cart2idx(uint k1, uint k2)
{
	uint ak1 = k1/2; uint ak2 = k2/2;
	uint64 r = ak1 + ak2;
	return 2*((uint64)ak1 + (r + r*r)/2) + (k2%2);
}

uint64       HamiltonianSystem::cart2idx_full(uint k1, uint k2, uint &parity_sector)
{
	uint s1 = k1%2; uint s2=k2%2; parity_sector = s1;
	uint ak1 = k1/2; uint ak2 = k2/2;
	uint64 r = ak1 + ak2;
	return 2*((uint64)ak1 + (r + r*r)/2) + s2;
}

void         HamiltonianSystem::idx2cart(uint64 idx, uint parity_sector, uint &k1, uint &k2)
{
	uint  idx0 = idx/2; uint s2 = idx%2;
	uint s1     = parity_sector;
	uint64 r    = (uint64)(floor(sqrt(2.0*idx0 + 0.25) - 0.5));
	uint64 k1k2 = (r*(r + 1))/2;
	k1 = 2*(uint)(idx0 - k1k2)     + s1;
	k2 = 2*(uint)(r + k1k2 - idx0) + s2;
}
	
HamiltonianSystem::HamiltonianSystem(int M, double L)
{
	this->L = L;
	this->M = M;
	NB  = M + M*M; //Size of the Hilbert space given the truncation parameter M
	NB2 = (uint64)NB*(uint64)NB;
	N = 0; 
	N2 = 0;
}

void  HamiltonianSystem::H(double* in, double* out)
{
	this->H(in    ,     out, 0);
	this->H(in + N, out + N, 1);
}

void HamiltonianSystem::SymmetryTransformation(uint ig, double* in, double* out, uint parity_sector)
{
	double* fin = new double[2*N];
	std::copy(in, in + N, fin + N*parity_sector); 
	std::fill(fin + N*(1-parity_sector), fin + N*(1 - parity_sector) + N, 0.0);
	SymmetryTransformation(ig, fin, out);
	delete [] fin;
};

void	HamiltonianSystem::CheckSymmetryGroupImplementation()
{
	//Checking the multiplication table for all irreps ...
	cout << ansi::green << "Checking multiplication rules for all " << ansi::magenta << NIRREPS << ansi::green << " irreps of the symmetry group of the " << ansi::magenta << model_name << ansi::green << " model with " << ansi::magenta << NG << ansi::green << " elements..." << ansi::reset << endl << flush;
	for(uint irrep=0; irrep<NIRREPS; irrep++)
	{
		uint    D  = irrep_dims[irrep];
		t_irrep T  = irreps[irrep];
		
		double* GiGj = new double[D*D];
		double  max_mult_err  = -1.0;
		for(uint i=0; i<NG; i++)
			for(uint j=0; j<NG; j++)
			{
				A_eq_B_mult_C(GiGj, (double*)(T + D*D*i), (double*)(T + D*D*j), D);
				uint k = group_mult_table[NG*i + j];
				max_mult_err = std::max(max_mult_err, norm_diff(GiGj, (double*)(T + D*D*k), D*D));
			};
		delete [] GiGj;
		
		cout << ansi::white << "\t Irrep " << ansi::magenta << irrep << ansi::white << ": max.err = " << ansi::magenta << max_mult_err << ansi::reset << endl << flush;
	};
	cout << ansi::green << " ...Done!" << ansi::reset << endl << endl << flush;  
	
	//Checking the Schur orthogonality relations for all irreps
	cout << ansi::green << "Checking Schur orthogonality relations among all irreps ..." << ansi::reset << endl << flush;
	double max_schur_err  = -1.0;
	for(uint ir1=0; ir1<NIRREPS; ir1++)
		for(uint ir2=0; ir2<NIRREPS; ir2++)
		{
			uint    D1 = irrep_dims[ir1]; uint    D2 = irrep_dims[ir2];
			t_irrep T1 = irreps[ir1];	  t_irrep T2 = irreps[ir2];
			
			for(uint ij1=0; ij1<D1*D1; ij1++)
				for(uint ij2=0; ij2<D2*D2; ij2++)
				{
					double res = 0.0;			
					for(uint ig=0; ig<NG; ig++)
						res += T1[D1*D1*ig + ij1]*T2[D2*D2*ig + ij2];
					if(ir1==ir2 && ij1==ij2)
					{
						max_schur_err = std::max(max_schur_err, std::abs(res - (double)NG/(double)D1));
					}
					else
						max_schur_err = std::max(max_schur_err, std::abs(res));
				};
		};
	cout << ansi::green << " ...Done! Max. error in Schur orthogonality: " << ansi::magenta << max_schur_err << ansi::reset << endl << endl << flush;
  
    cout << ansi::green << "Checking the representation of the symmetry group in the basis of harmonic oscillator eigenfunctions ..." << ansi::reset << endl << flush;
	double max_func_rep_err  = -1.0;
	double* v1 = new double[2*N]; double* v2 = new double[2*N]; double* v3 = new double[2*N]; double* v4 = new double[2*N];
	rand_vec(v1, 2*N);
	for(uint ig1=0; ig1<NG; ig1++)
		for(uint ig2=0; ig2<NG; ig2++)
		{
			SymmetryTransformation(ig2, v1, v2);
			SymmetryTransformation(ig1, v2, v3); //v3 = ig1*ig2*v1;
   			
			uint ig12 = group_mult_table[ig1*NG + ig2];
			SymmetryTransformation(ig12, v1, v2); //v2 = (ig1*ig2)*v1;
   			
			max_func_rep_err = std::max(max_func_rep_err, norm_diff(v2, v3, 2*N));
		};
	cout << ansi::green << " ...Done! Max. error of the functional representation: " << ansi::magenta << max_func_rep_err << ansi::reset << endl << endl << flush;
	
	cout << ansi::green << "Checking that symmetry transformations indeed commute with the Hamiltonian ..." << ansi::reset << endl << flush;
	double max_commutativity_err  = -1.0;
	rand_vec(v1, 2*N);
	for(uint ig=0; ig<NG; ig++)
    {
		SymmetryTransformation(ig, v1, v2); //v2 = G*v1
		H(v2, v3); //v3 = H*v2 = H*G*v1;
		
		H(v1, v2); //v2 = H*v1
		SymmetryTransformation(ig, v2, v4); //v4 = G*v2 = G*H*v1
   		 			
		max_commutativity_err = std::max(max_commutativity_err, norm_diff(v3, v4, 2*N));
	};
	cout << ansi::green << " ...Done! Max. norm of [G, H]: " << ansi::magenta << max_commutativity_err << ansi::reset << endl << endl << flush;
	delete [] v1; delete [] v2; delete [] v3; delete [] v4;
}

BosonicSystem::BosonicSystem(int M, double L, bool free_system)  : HamiltonianSystem(M, L)
{
	N   = NB;
	N2  = (uint64)N*(uint64)N;
	NG  = 8;
	NIRREPS = 5;
	irreps     = new t_irrep[NIRREPS];
	irrep_dims = new    uint[NIRREPS];
	group_mult_table = C4v_Mult_Table;
	
	irreps[0] = C4v_A1; irrep_dims[0] = 1;
	irreps[1] = C4v_A2; irrep_dims[1] = 1;
	irreps[2] = C4v_B1; irrep_dims[2] = 1;
	irreps[3] = C4v_B2; irrep_dims[3] = 1;
	irreps[4] = C4v_E0; irrep_dims[4] = 2;
	
	this->free_system = free_system;
	
	sprintf(model_name, (free_system? "free" : "bosonic")); 
}

void BosonicSystem::SymmetryTransformation(uint ig, double* in, double* out)
{
	//Construct npsi after acting with a group element number ig
	//First we apply the coordinate transformation
	std::fill(out, out + 2*N, 0.0);
	//#pragma omp parallel for
	for(uint parity_sector=0; parity_sector<2; parity_sector ++)
		for(uint idx=0; idx<NB; idx++)
		{
			uint k1, k2, nk1, nk2, new_parity_sector = 0, nidx;
				
			idx2cart(idx, parity_sector, k1, k2);
			double phase = C4v_k1k2_fwd(ig, k1,  k2, nk1, nk2);
			
			nidx = cart2idx_full(nk1, nk2, new_parity_sector);
			out[new_parity_sector*N + nidx] += phase*in[parity_sector*N + idx];
		};
}

void BosonicSystem::H(double* in, uint parity_sector) //Hamiltonian matrix in the parity sector (s1, s2)
{
	#pragma omp parallel for
	for(uint64 idxk=0; idxk<NB; idxk++)
		for(uint64 idxl=0; idxl<NB; idxl++)
		{
			uint64 idxkl = idxk*NB + idxl;
			in[idxkl] = 0.0;
			uint k1, k2, l1, l2;
			idx2cart(idxk, parity_sector, k1, k2);
			idx2cart(idxl, parity_sector, l1, l2);
			//Kinetic terms in the Hamiltonian
			if(k1==l1) in[idxkl] += psq(k2, l2);
			if(k2==l2) in[idxkl] += psq(k1, l1);
			//Potential term in the Hamiltonian
			if(!free_system) in[idxkl] += xsq(k1, l1)*xsq(k2, l2);
		};	
}

double* HamiltonianSystem::H(uint parity_sector)
{
	double* res = new double[N2];
    H(res, parity_sector);
	return res;
}

void        BosonicSystem::H(double* in, double* out, uint parity_sector)
{
	uint s1 = parity_sector%2;
	uint s2 = (parity_sector/2)%2;
	//Some helpful variables
	double i2L  = 0.5/(L*L);
	double L4i4 = (free_system? 0.0 : 0.25*L*L*L*L);
	#pragma omp parallel for
	for(uint64 idxk=0; idxk<NB; idxk++)
	{
		uint k1, k2; uint64 idxl;
		idx2cart(idxk, parity_sector, k1, k2);
		//Some commonly used algebraic constants
		double a1  = (double)(2*k1 + 1);
		double a2  = (double)(2*k2 + 1);
		double s1p = sqrt((double)(k1 + 1)*(double)(k1 + 2));
		double s1m = sqrt((double)(k1    )*(double)(k1 - 1));
		double s2p = sqrt((double)(k2 + 1)*(double)(k2 + 2));
		double s2m = sqrt((double)(k2    )*(double)(k2 - 1));
		//Diagonal terms
		out[idxk]  = (i2L*(a1 + a2) + L4i4*a1*a2)*in[idxk];
		//l1 = k1 + 2, l2=k2 
		idxl = cart2idx(k1+2, k2);
		if(idxl < NB) out[idxk] += (L4i4*a2 - i2L)*s1p*in[idxl];
		//l1 = k1 - 2, l2=k2 
		if(k1>1)
		{
			idxl = cart2idx(k1-2, k2);
			out[idxk] += (L4i4*a2 - i2L)*s1m*in[idxl];
		};
		//l1=k1, l2 = k2 + 2
		idxl = cart2idx(k1, k2+2);
		if(idxl < NB) out[idxk] += (L4i4*a1 - i2L)*s2p*in[idxl];
		//l1=k1, l2 = k2 - 2
		if(k2>1)
		{
			idxl = cart2idx(k1, k2-2);
			out[idxk] += (L4i4*a1 - i2L)*s2m*in[idxl];
		};
		//l1 = k1+2, l2=k2+2
		idxl = cart2idx(k1+2, k2+2);
		if(idxl < NB) out[idxk] += L4i4*s1p*s2p*in[idxl];
		//l1 = k1+2, l2=k2-2
		if(k2>1)
		{
			idxl = cart2idx(k1+2, k2-2);
			if(idxl < NB) out[idxk] += L4i4*s1p*s2m*in[idxl];
		};
		if(k1>1)
		{
			//l1 = k1-2, l2=k2+2
			idxl = cart2idx(k1-2, k2+2);
			if(idxl < NB) out[idxk] += L4i4*s1m*s2p*in[idxl];
			//l1 = k1-2, l2=k2-2
			if(k2>1)
			{
				idxl = cart2idx(k1-2,k2-2);
				out[idxk] += L4i4*s1m*s2m*in[idxl];
			};
		};
	};		
}

void        HamiltonianSystem::FindAllEigenstates(uint parity_sector)
{
	delete [] psi; psi = new double[N2];
	delete [] E;   E   = new double[N];
	
	//Initialize storage for evecs
	H(psi, parity_sector); //Temporarily
	
	/*cout << "Matrix elements of H: " << endl;
	for(uint i=0; i<N; i++)
	{
		cout << std::fixed << std::setprecision(3) << std::showpos;
		for(uint j=0; j<N; j++)
			cout << psi[i*N + j] << " ";
		cout << endl << flush; 
	};
	cout << endl;*/
	
	cout << "Hi from FindAllEigenstates with parity_sector = " << ansi::cyan << parity_sector << ansi::reset << endl << flush;
	int res = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N, psi, N, E);
	if(res!=0) cout << endl << ansi::red << "Something went wrong, LAPACKE_dsyev returned " << res << " !!!" << endl << flush;
		cout << "Now lapacke_dsyev finished ..." << endl << flush;
	//Transposing psi to make it into an eigensystem
	transpose(psi, N);
	sort_eigensystem(E, psi, N, N);
	this->nev = N;
}

void       HamiltonianSystem::ReorthogonalizeEigenstates(double precision)
{
	uint ie = 0;
	while(ie<nev)
	{
		uint nd = 1;
		while(ie + nd < nev && std::abs(scalar_prod(psi + N*ie, psi + N*(ie + nd), N))>precision)
			nd ++;
		if(nd>1) //Re-orthogonalize eigenvectors
		{
			double* spm = new double[nd*nd];
			for(uint i=0; i<nd; i++)
				for(uint j=0; j<nd; j++)
					spm[i*nd + j] = scalar_prod(psi + N*(ie + j), psi + N*(ie + i), N); //Note the reversed indices i and j in scalar_prod
			double* spme = new double[nd];
			int res = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nd, spm, nd, spme);
			if(res!=0) cout << endl << ansi::red << "Something went wrong, LAPACKE_dsyev in orthogonalization returned " << res << " !!!" << endl << flush;
			//Reconstructing spm
			double* rwt = new double[nd*nd];
			for(uint i=0; i<nd; i++)
				for(uint j=0; j<nd; j++)
				{
					rwt[i*nd + j] = 0.0;
					for(uint k=0; k<nd; k++)
						rwt[i*nd + j] += spm[i*nd + k]*spm[j*nd + k]/sqrt(spme[k]);
				};
				//Reweighting the vectors
				double* new_evecs = new double[N*nd];
				std::fill(new_evecs, new_evecs + N*nd, 0.0);
				for(uint i=0; i<nd; i++)
					for(uint j=0; j<nd; j++)
						A_pluseq_bB(new_evecs + N*i, rwt[i*nd + j], psi + N*(ie + j), N);
				std::copy(new_evecs, new_evecs + N*nd, psi + N*ie);
				delete [] new_evecs;
				delete [] rwt;
				delete [] spm;
				delete [] spme;
		};
		ie += nd;
	};	
}

int        HamiltonianSystem::FindLowestEigenstates(uint nev, uint parity_sector, double precision, uint ncv)
{
	if(psi!=NULL) delete [] psi; psi = new double[nev*N];
	if(E!=NULL)   delete [] E;   E   = new double[nev];
	int myncv = (ncv==0? 2*nev + 1 : ncv);
	
	uint niter = 0; bool success = true;
	
	ARrcSymStdEig<double>* prob = new ARrcSymStdEig<double>(N, nev, "SA", myncv, precision, INT_MAX);
	uint64 itc = 0;
	while (!prob->ArnoldiBasisFound())
	{
	   	prob->TakeStep();
		if ((prob->GetIdo() == 1)||(prob->GetIdo() == -1)) 
		{
			H(prob->GetVector(), prob->PutVector(), parity_sector);
			itc++;
       	};
    };
    cout << ansi::white << "\t No. of iterations: " << ansi::magenta << prob->GetIter() << ansi::white << " (my own counter " << ansi::magenta << itc << ansi::white << ")" << ansi::reset << flush << endl;
    	// Finding eigenvalues and eigenvectors.
    prob->FindEigenvectors();
    niter = prob->GetIter();
    if(prob->EigenvaluesFound() && prob->ConvergedEigenvalues()==nev)
    {
    	cout << ansi::yellow << "\t Starting reorthogonalization... " << ansi::reset << endl << flush;
		prob->EigenValVectors(psi, E);
		sort_eigensystem(E, psi, N, nev);
		ReorthogonalizeEigenstates(10.0*precision);
	};
	delete prob;
	this->nev = nev;
	return niter;
}

void        HamiltonianSystem::CheckEigensystem(uint parity_sector, double* eigensys_err, double* orthogonality_err)
{
	double* tmpvec = new double[N];
	double max_eigensys_err = 0.0;
	double max_ortho_err    = 0.0;
	
	//Eigenvector error
	for(uint ie=0; ie<nev; ie++)
	{
		H(psi + N*ie, tmpvec, parity_sector);
		A_pluseq_bB(tmpvec, -E[ie], psi + N*ie, N);
		max_eigensys_err = std::max(max_eigensys_err, norm(tmpvec, N));
	};

	//Orthonormality
	for(uint ie1=0; ie1<nev; ie1++)
		for(uint ie2=0; ie2<nev; ie2++)
		{
			double sp = scalar_prod(psi + N*ie1, psi + N*ie2, N) - (ie1==ie2? 1.0 : 0.0 );
			max_ortho_err = std::max(max_ortho_err, std::abs(sp));
		};
	
	delete [] tmpvec;
	
	if(eigensys_err!=NULL)
		*eigensys_err = max_eigensys_err;
	else
		cout << ansi::white << "\t Max. eigensystem error (N = " << ansi::magenta << N << ansi::white <<"):\t " << ansi::magenta << max_eigensys_err << ansi::reset << endl; 

	if(orthogonality_err!=NULL)
		*orthogonality_err = max_ortho_err;
	else
		cout << ansi::white << "\t Max. orthogonality error(N = " << ansi::magenta << N << ansi::white <<"):\t " << ansi::magenta << max_ortho_err << ansi::reset << endl; 	
}

void        HamiltonianSystem::CheckHermiticity(uint parity_sector, double* err)
{
	double max_err  = 0.0;
	double* HM = H(parity_sector);
	for(uint i=0; i<N; i++)
		for(uint j=0; j<N; j++)
		{
			double aerr = HM[i*N + j] - HM[j*N + i];
			max_err = std::max(max_err, std::abs(aerr));	
		};
	delete [] HM;	
	if(err!=NULL)
		*err = max_err;
	else
		cout << ansi::white << "\t Max. hermiticity error (N = " << ansi::magenta << N << ansi::white <<"):\t " << ansi::magenta << max_err << ansi::reset << endl; 	
};

uint         HamiltonianSystem::read_eigensystem(string datadir, uint parity_sector, bool read_evecs) //Returns the number of eigenvalues/vectors that were read in successfully
{
	uint fnev = 0;
	
	char evals_fname[512];
	sprintf(evals_fname, "%s/eigensystems/evals_%s_M%u_L%1.4lf_p%u.dat", datadir.c_str(), model_name, M, L, parity_sector);
	FILE* evals_file = fopen(evals_fname, "rb"); //Binary read
	if(evals_file!=NULL)
	{
		//Check the file size and determine how many eigenvalues it contains
		fseek(evals_file, 0L, SEEK_END);
		size_t fs = ftell(evals_file);
		rewind(evals_file);
		fnev = fs/sizeof(double);
		delete [] E;
		E = new double[fnev];
		if(fread(E, sizeof(double), fnev, evals_file)!=fnev)
		{
			cout << ansi::red << "Some error, Could not read " << fnev << " eigenvalues from the file " << evals_fname << "!!!" << ansi::reset << endl << flush;
			fclose(evals_file);
			return 0;
		};
		cout << ansi::white << "\t Loaded " << ansi::magenta << fnev << ansi::white << " eigenvalues from the file " << ansi::magenta << evals_fname << ansi::reset << endl << flush;
		fclose(evals_file);
	}
	else
	{
		cout << ansi::yellow << "The file " << ansi::magenta << evals_fname << ansi::yellow << " could not be opened!!!" << ansi::reset << endl << flush;
		return 0;
	};
	
	if(read_evecs)
	{
		char evecs_fname[512];
		sprintf(evecs_fname, "%s/eigensystems/evecs_%s_M%u_L%1.4lf_p%u.dat", datadir.c_str(), model_name, M, L, parity_sector);
		FILE* evecs_file = fopen(evecs_fname, "rb"); //Binary write
		if(evecs_file!=NULL)
		{
			delete [] psi;
			psi = new double[fnev*N];
			cout << ansi::white << "\t Reading " << ansi::magenta << fnev << ansi::white << " eigenvectors from the file " << ansi::magenta << evecs_fname << ansi::white << " ..."<< ansi::reset << endl << flush;
			TIMING_INIT;
			TIMING_START;
			if(fread(psi, sizeof(double), N*fnev, evecs_file)!=N*fnev)
			{
				cerr << ansi::red << "Failed to read " << fnev << " eigenvectors from the file " << evecs_fname << "!!!" << ansi::reset << endl << flush;
				return 0;
			};
			fclose(evecs_file);
			TIMING_STOP;
			cout << ansi::white << "\t ... Done in " << ansi::magenta << a_time << ansi::white << " sec." << ansi::reset << endl << flush;
		}
		else
		{
			cerr << ansi::red << "Failed to open the file " << evals_fname << " to read the eigenvectors!" << ansi::reset << endl << flush;
			return 0;
		};
	};
	
	this->nev = fnev;
	return this->nev;
}

void        HamiltonianSystem::write_eigensystem(string datadir, uint parity_sector, bool write_evecs)
{
	char evals_fname[512];
	sprintf(evals_fname, "%s/eigensystems/evals_%s_M%u_L%1.4lf_p%u.dat", datadir.c_str(), model_name, M, L, parity_sector);
	cout << ansi::white << "\t File to save eigenvalues: " << ansi::magenta << evals_fname << ansi::reset << flush << endl;
	FILE* evals_file = fopen(evals_fname, "wb"); //Binary write
	if(evals_file!=NULL)
	{
		size_t res = fwrite(E, sizeof(double), nev, evals_file);
		if(res!=nev)
			cerr << ansi::red << "Only " << res << " doubles out of " << nev << " could be written to the file " << evals_fname << "!!!" << ansi::reset << endl << flush;
		fclose(evals_file);
	}
	else
		cerr << ansi::red << "Failed to open the file " << evals_fname << " for binary writing!" << ansi::reset << endl << flush;
	
	if(!write_evecs) return;
		
	char evecs_fname[512];
	sprintf(evecs_fname, "%s/eigensystems/evecs_%s_M%u_L%1.4lf_p%u.dat", datadir.c_str(), model_name, M, L, parity_sector);
	cout << ansi::white << "\t File to save eigenvectors: " << ansi::magenta << evecs_fname << ansi::reset << flush << endl;
	FILE* evecs_file = fopen(evecs_fname, "wb"); //Binary write
	if(evecs_file!=NULL)
	{
		if(fwrite(psi, sizeof(double), nev*N, evecs_file)!=nev*N)
			cerr << ansi::red << "Failed to write " << nev*N << " doubles to the file " << evecs_fname << "!!!" << ansi::reset << endl << flush;
		fclose(evecs_file);
	}
	else
		cerr << ansi::red << "Failed to open the file " << evals_fname << " for binary writing!" << ansi::reset << endl << flush;	 
}

double      HamiltonianSystem::x(uint k, uint l)
{
	if(l==k+1) return  sqrt(0.5*(double)(k+1))*L;
	if(l==k-1) return  sqrt(0.5*(double)(k  ))*L;
	return 0.0;
}

double      HamiltonianSystem::ip(uint k, uint l)
{
	if(l==k+1) return  -sqrt(0.5*(double)(k+1))/L;
	if(l==k-1) return   sqrt(0.5*(double)(k  ))/L;
	return 0.0;
}

double      HamiltonianSystem::xsq(uint k, uint l)
{
	if(l==k)   return  0.5*(double)(2*k + 1)*L*L;
	if(l==k+2) return  0.5*sqrt((double)((k+1)*(k+2)))*L*L;
	if(l==k-2) return  0.5*sqrt((double)((k  )*(k-1)))*L*L;
	return 0.0;
}

double      HamiltonianSystem::psq(uint k, uint l)
{
	if(l==k)   return   0.5*(double)(2*k + 1)/(L*L);
	if(l==k+2) return  -0.5*sqrt((double)((k+1)*(k+2)))/(L*L);
	if(l==k-2) return  -0.5*sqrt((double)((k  )*(k-1)))/(L*L);
	return 0.0;
}

void BosonicSystem::pluseq_x1(double* in, double* out, uint parity_sector_in, double c)
{
	double aL = L/sqrt(2.0);
	#pragma omp parallel for
	for(uint64 idxk=0; idxk<NB; idxk++)
	{
		uint k1, k2; uint64 idxl;
		idx2cart(idxk, 1-parity_sector_in, k1, k2);
		//l1 = k1 + 1, l2=k2
		idxl = cart2idx(k1 + 1, k2);
		if(idxl<NB) out[idxk] += c*aL*sqrt((double)(k1 + 1))*in[idxl];
		//l1 = k1 - 1, l2=k2 
		if(k1>0)
		{
			idxl = cart2idx(k1 - 1, k2);
			out[idxk] += c*aL*sqrt((double)k1)*in[idxl];
		};
	};
}

void BosonicSystem::pluseq_x2(double* in, double* out, uint parity_sector, double c)
{
	double aL = L/sqrt(2.0);
	#pragma omp parallel for
	for(uint64 idxk=0; idxk<NB; idxk++)
	{
		uint k1, k2; uint64 idxl;
		idx2cart(idxk, parity_sector, k1, k2);
		//l1 = k1, l2=k2+1
		idxl = cart2idx(k1, k2+1);
		if(idxl<NB) out[idxk] += c*aL*sqrt((double)(k2 + 1))*in[idxl];
		//l1 = k1, l2=k2-1 
		if(k2>0)
		{
			idxl = cart2idx(k1, k2-1);
			out[idxk] += c*aL*sqrt((double)k2)*in[idxl];
		};
	};
}

void BosonicSystem::pluseq_ip2(double* in, double* out, uint parity_sector, double c)
{
	double aL = 1.0/(sqrt(2.0)*L);
	#pragma omp parallel for
	for(uint64 idxk=0; idxk<NB; idxk++)
	{
		uint k1, k2; uint64 idxl;
		idx2cart(idxk, parity_sector, k1, k2);
		//l1 = k1, l2=k2+1
		idxl = cart2idx(k1, k2+1);
		if(idxl<NB) out[idxk] += c*aL*sqrt((double)(k2 + 1))*in[idxl];
		//l1 = k1, l2=k2-1 
		if(k2>0)
		{
			idxl = cart2idx(k1, k2-1);
			out[idxk] -= c*aL*sqrt((double)k2)*in[idxl];
		};
	};
}

void BosonicSystem::x2(double* in, double* out, uint parity_sector)
{
	std::fill(out, out + NB, 0.0);
	pluseq_x2(in,  out, parity_sector);
}

void BosonicSystem::ip2(double* in, double* out, uint parity_sector)
{
	std::fill(out, out + NB, 0.0);
	pluseq_ip2(in, out, parity_sector);
}

SupersymmetricSystem::SupersymmetricSystem(int M, double L)  : BosonicSystem(M, L)
{
	N   = 2*NB;
	N2  = (uint64)N*(uint64)N;
	NG      = 16;
	NIRREPS = 7; //Should be 7, but for debugging purposes let us only consider E1 and E2
	irreps     = new t_irrep[NIRREPS];
	irrep_dims = new    uint[NIRREPS];
	group_mult_table = D4d_Mult_Table;
	
	irreps[0] = D4d_A1; irrep_dims[0] = 1;
	irreps[1] = D4d_A2; irrep_dims[1] = 1;
	irreps[2] = D4d_B1; irrep_dims[2] = 1;
	irreps[3] = D4d_B2; irrep_dims[3] = 1;
	irreps[4] = D4d_E0; irrep_dims[4] = 2;
	irreps[5] = D4d_E1; irrep_dims[5] = 2;
	irreps[6] = D4d_E2; irrep_dims[6] = 2;
	
	sprintf(model_name, "supersymmetric");	
}

void        SupersymmetricSystem::H(double* in, uint parity_sector)	//The matrix in contains N*N Hamiltonian matrix
{
	#pragma omp parallel for
	for(uint64 idxk=0; idxk<NB; idxk++)
		for(uint64 idxl=0; idxl<NB; idxl++)
		{
			uint64 idxkl11 = (idxk     )*N + (idxl     );
			uint64 idxkl12 = (idxk     )*N + (idxl + NB);
			uint64 idxkl21 = (idxk + NB)*N + (idxl     );
			uint64 idxkl22 = (idxk + NB)*N + (idxl + NB);
			
			in[idxkl11] = 0.0; in[idxkl12] = 0.0; in[idxkl21] = 0.0; in[idxkl22] = 0.0;
			
			uint k1, k2, l1, l2;
			
			//H_11 - the first diagonal block 
			idx2cart(idxk, parity_sector, k1, k2);
			idx2cart(idxl, parity_sector, l1, l2);
			if(k1==l1) in[idxkl11] += psq(k2, l2) + x(k2, l2); //Kinetic term p2^2 + x2
			if(k2==l2) in[idxkl11] += psq(k1, l1);             //Kinetic term p1^2
			in[idxkl11] += xsq(k1, l1)*xsq(k2, l2);
			
			//H_22 - the second diagonal block 
			idx2cart(idxk, 1-parity_sector, k1, k2);
			idx2cart(idxl, 1-parity_sector, l1, l2);
			if(k1==l1) in[idxkl22] += psq(k2, l2) - x(k2, l2); //Kinetic term p2^2 - x2
			if(k2==l2) in[idxkl22] += psq(k1, l1);             //Kinetic term p1^2
			in[idxkl22] += xsq(k1, l1)*xsq(k2, l2);
						
			//H12 = x1 - the first off-diagonal block
			idx2cart(idxk,   parity_sector, k1, k2);
			idx2cart(idxl, 1-parity_sector, l1, l2);
			if(k2==l2) in[idxkl12] += x(k1, l1);
			
			//H21 = x1 - the second off-diagonal block
			idx2cart(idxk, 1-parity_sector, k1, k2);
			idx2cart(idxl,   parity_sector, l1, l2);
			if(k2==l2) in[idxkl21] += x(k1, l1); 
		};
}

void        SupersymmetricSystem::H(double* in, double* out, uint parity_sector) //Hamiltonian, in and out are real since the Hamiltonian is real
{
	//First diagonal block H_11
	BosonicSystem::H(in     , out     ,     parity_sector); 
	//First diagonal block - first NB components of in and out
	       pluseq_x2(in     , out     ,     parity_sector, +1.0);
	//Second diagonal block H_12
	BosonicSystem::H(in + NB, out + NB, 1 - parity_sector); //Second diagonal block - last NB components of in and out
	       pluseq_x2(in + NB, out + NB, 1 - parity_sector, -1.0);
	//Off-diagonal blocks H_12 = H_21 = x1
	       pluseq_x1(in + NB, out     , 1 - parity_sector, 1.0);
	       pluseq_x1(in     , out + NB,     parity_sector, 1.0);
}

void SupersymmetricSystem::x2(double* in, double* out, uint parity_sector)
{
	BosonicSystem::x2(in     , out     ,   parity_sector);
	BosonicSystem::x2(in + NB, out + NB, 1-parity_sector);
}

void SupersymmetricSystem::ip2(double* in, double* out, uint parity_sector)
{
	BosonicSystem::ip2(in     , out     ,   parity_sector);
	BosonicSystem::ip2(in + NB, out + NB, 1-parity_sector);
}

void     HamiltonianSystem::write_OTOC_data(string datadir, uint parity_sector)
{
	char fname[512];
	sprintf(fname, "%s/otoc_data_%s_M%u_L%1.4lf_p%u_n%u.tmp", datadir.c_str(), model_name, M, L, parity_sector, nev);
	cout << ansi::white << ansi::white << "\t File to save operator matrix elements for OTOCs: " << ansi::magenta << fname << ansi::reset << flush << endl;
	FILE* ofile = fopen(fname, "wb"); //Binary write
	if(ofile!=NULL)
	{
		size_t resA    = fwrite(A, sizeof(double), nev*nev, ofile);
		if(resA != nev*nev)
			cerr << ansi::red << "Only " << resA << " matrix elements of A out of " << nev*nev << " could be written to the file " << fname << "!!!" << ansi::reset << endl << flush;
		size_t resB    = fwrite(B, sizeof(double), nev*nev, ofile); 	
		if(resB != nev*nev)
			cerr << ansi::red << "Only " << resB << " matrix elements of B out of " << nev*nev << " could be written to the file " << fname << "!!!" << ansi::reset << endl << flush;
		fclose(ofile);
		cout << ansi::white << "\t\t ... Data saved! " << ansi::reset << flush << endl;
	}
	else
		cerr << ansi::red << "Failed to open the file " << ofile << " for binary writing!" << ansi::reset << endl << flush;
}

int     HamiltonianSystem::read_OTOC_data(string datadir, uint parity_sector, uint fnev)
{
	FILE* ifile = NULL; uint an = 2000; char fname[512];
	//First try exactly fnev
	sprintf(fname, "%s/otoc_data_%s_M%u_L%1.4lf_p%u_n%u.tmp", datadir.c_str(), model_name, M, L, parity_sector, fnev);	
	ifile = fopen(fname, "rb"); //Binary read
	cout << ansi::yellow << "Attempted to open the file " << ansi::magenta << fname << ansi::yellow << ", fopen returned " << ansi::magenta << ifile << ansi::reset << endl << flush;
	//Second, try exactly nev
	if(ifile==NULL)
	{
		sprintf(fname, "%s/otoc_data_%s_M%u_L%1.4lf_p%u_n%u.tmp", datadir.c_str(), model_name, M, L, parity_sector, this->nev);	
		ifile = fopen(fname, "rb"); //Binary read
		cout << ansi::yellow << "Attempted to open the file " << ansi::magenta << fname << ansi::yellow << ", fopen returned " << ansi::magenta << ifile << ansi::reset << endl << flush;
	};
	
    //If the file with exactly nev eigenvalues not found - try other values of nev...
	while(ifile==NULL && an>0)
	{
		sprintf(fname, "%s/otoc_data_%s_M%u_L%1.4lf_p%u_n%u.tmp", datadir.c_str(), model_name, M, L, parity_sector, an);
		ifile = fopen(fname, "rb"); //Binary read
		an -= 50;
	};
	
	if(ifile!=NULL)
	{
		//First determine the data size
		fseek(ifile, 0L, SEEK_END);
		size_t fs = ftell(ifile);
		rewind(ifile);
		uint fnev2 = fs/(2*sizeof(double));
		uint fnev  = (int)round(sqrt((double)fnev2));
		if(fnev2 != fnev*fnev)
		{
			cout << ansi::red << "The file " << ansi::yellow << fname << ansi::red << " is inconsistent, fnev2 = " << ansi::yellow << fnev2 << ansi::red << "is not a square of an integer number!!!" << ansi::reset << endl << flush;
			fclose(ifile);
			return -1;
		};
		
		delete [] A; A = new double[fnev*fnev];
		if(fread(A, sizeof(double), fnev*fnev, ifile) != fnev*fnev)
		{
			cout << ansi::red << "Some error, Could not read " << fnev*fnev << " matrix elements of A from the file " << fname << "!!!" << ansi::reset << endl << flush;
			fclose(ifile);
			return -2;
		};
		
		delete [] B; B = new double[fnev*fnev];
		if(fread(B, sizeof(double), fnev*fnev, ifile) != fnev*fnev)
		{
			cout << ansi::red << "Some error, Could not read " << fnev*fnev << " matrix elements of B from the file " << fname << "!!!" << ansi::reset << endl << flush;
			fclose(ifile);
			return -3;
		};
		fclose(ifile);
		cout << ansi::white << "\t Matrix elements for the first " << ansi::magenta << fnev << ansi::white << " eigenstates loaded from the file " << ansi::magenta << fname << ansi::reset << flush << endl << endl;
		
		this->nev = fnev;
	}
	else
	{
		cout << ansi::yellow << "The file " << ansi::magenta << fname << ansi::yellow << " could not be opened!!!" << ansi::reset << endl << flush;
		return -4;
	};
	return this->nev;
}

void    HamiltonianSystem::init_OTOCS(uint parity_sector)
{
	cout << ansi::white << "\t Initializing OTOC with " << ansi::magenta << this->nev << ansi::white << " eigenstates" << ansi::reset << flush << endl; 
	
	delete [] A; delete [] B;
	A = new double[nev*nev]; B = new double[nev*nev];
	
	TIMING_INIT;
	TIMING_START;
	double* tmpvec = new double[N];
	
	for(uint ie2=0; ie2<nev; ie2++)
	{
		ip2(psi + N*ie2, tmpvec, parity_sector);
		for(uint ie1=0; ie1<nev; ie1++)
			A[ie1*nev + ie2] = scalar_prod(psi + N*ie1, tmpvec, N);
	
		x2(psi + N*ie2, tmpvec, parity_sector);
		for(uint ie1=0; ie1<nev; ie1++)	
			B[ie1*nev + ie2] = scalar_prod(psi + N*ie1, tmpvec, N);
	};
	delete [] tmpvec;
	TIMING_STOP;
}

void  HamiltonianSystem::OTOC_regularized(double beta, double t, double* r1, double* r2, double* r3)
{
	t_complex* Atr = new t_complex[nev*nev];
	t_complex* Br  = new t_complex[nev*nev];
	t_complex* tmp = new t_complex[nev*nev];
	
	if(r1!=NULL) *r1 = 0.0; 
	if(r2!=NULL) *r2 = 0.0; 
	if(r3!=NULL) *r3 = 0.0;

	double Z = 0.0;	
	#pragma omp parallel for reduction(+ : Z)
	for(uint ie=0; ie<nev; ie++)
		Z += exp(-beta*E[ie]);
	
	#pragma omp parallel for
	for(uint ie1=0; ie1<nev; ie1++)
		for(uint ie2=0; ie2<nev; ie2++)
		{
			Atr[ie1*nev + ie2] = A[ie1*nev + ie2]*exp(1.0i*t*(E[ie1] - E[ie2]) - 0.125*beta*(E[ie1] + E[ie2]));
			Br[ ie1*nev + ie2] = B[ie1*nev + ie2]*exp(                         - 0.125*beta*(E[ie1] + E[ie2]));
		};
		
	A_eq_B_mult_C(tmp, Atr, Br, nev);

	if(r1!=NULL)
	{
		double ar1 = 0.0;  //r1 is F(t) in the notation of arXiv:2109.08693
		#pragma omp parallel for reduction(+ : ar1)
		for(uint ie=0; ie<nev; ie++)
			for(uint k=0; k<nev; k++)
				ar1 += real(tmp[ie*nev + k]*tmp[k*nev + ie]);
		*r1 += ar1;
	};
	
	if(r2!=NULL)
	{
		double ar2 = 0.0; //r2 is eq (2.5) in 2109.08693, but with "y" factors between all operators
		A_eq_B_mult_C(tmp, Atr, Atr, nev);
		A_eq_B_mult_C(Atr,  Br,  Br, nev); //Temporary using Atr for Br*Br
		#pragma omp parallel for reduction(+ : ar2)
		for(uint ie=0; ie<nev; ie++)
			for(uint k=0; k<nev; k++)
				ar2 += real(tmp[ie*nev + k]*Atr[k*nev + ie]);
		*r2 += ar2;
	};
	
	if(r3!=NULL) //r3 is eq (2.5) in 2100.08693 EXACTLY
	{
		double ar3 = 0.0;
		#pragma omp parallel for
		for(uint i=0; i<nev; i++)
			for(uint j=0; j<nev; j++)
				Atr[i*nev + j] = exp(1.0i*t*(E[i] - E[j]))*A[i*nev + j];
		#pragma omp parallel for
		for(uint i=0; i<nev; i++)
			for(uint j=0; j<nev; j++)
				Br[i*nev + j] = B[i*nev + j] + 0.0i;
					
		A_eq_B_mult_C(tmp, Atr, Br, nev);
		#pragma omp parallel for reduction(+ : ar3)
		for(uint i=0; i<nev; i++)
			for(uint j=0; j<nev; j++)
			{
				t_complex Mij = tmp[i*nev + j]*exp(-0.25*beta*(E[i] + E[j]));
				ar3 += norm(Mij);
			};
		*r3 += ar3;
	};
	
	if(r1!=NULL) *r1 /= (-Z); 
	if(r2!=NULL) *r2 /= (-Z);
	if(r3!=NULL) *r3 /= (-Z);
			
	delete [] Atr; delete [] tmp; delete [] Br;
}

void       HamiltonianSystem::thermal_vevs_regularized(double beta, double& aA, double& aB, double& Z)
{
	aA = 0.0; aB = 0.0; Z = 0.0;
	
	double rA = 0.0;
	#pragma omp parallel for reduction(+ : rA)
	for(uint ie1=0; ie1<nev; ie1++)
		for(uint ie2=0; ie2<nev; ie2++)
			rA += A[ie1*nev + ie2]*A[ie2*nev + ie1]*exp(-0.5*beta*(E[ie1] + E[ie2]));
	aA += rA;
	
	double rB = 0.0;		
	#pragma omp parallel for reduction(+ : rB)
	for(uint ie1=0; ie1<nev; ie1++)
		for(uint ie2=0; ie2<nev; ie2++)
			rB += B[ie1*nev + ie2]*B[ie2*nev + ie1]*exp(-0.5*beta*(E[ie1] + E[ie2]));
	aB += rB;
			
	double aZ = 0.0;
	#pragma omp parallel for reduction(+ : aZ)
	for(uint ie=0; ie<nev; ie++)
		aZ += exp(-beta*E[ie]);
	Z += aZ;
	
	aA /= (+Z); aB /= (-Z); 
}

void HamiltonianSystem::thermal_vevs(double beta, double& aA, double& aB, double& Z)
{
	aA = 0.0; aB = 0.0; Z = 0.0;
	
	double rA = 0.0;
	#pragma omp parallel for reduction(+ : rA)
	for(uint ie1=0; ie1<nev; ie1++)
		for(uint ie2=0; ie2<nev; ie2++)
			rA += A[ie1*nev + ie2]*A[ie2*nev + ie1]*exp(-beta*E[ie1]);
	aA += rA;
	
	double rB = 0.0;		
	#pragma omp parallel for reduction(+ : rB)
	for(uint ie1=0; ie1<nev; ie1++)
		for(uint ie2=0; ie2<nev; ie2++)
			rB += B[ie1*nev + ie2]*B[ie2*nev + ie1]*exp(-beta*E[ie1]);
	aB += rB;
			
	double aZ = 0.0;
	#pragma omp parallel for reduction(+ : aZ)
	for(uint ie=0; ie<nev; ie++)
		aZ += exp(-beta*E[ie]);
	Z += aZ;
	
	aA /= (+Z); aB /= (-Z); 
}

void   HamiltonianSystem::resize_OTOC_data(uint new_nev)
{
	if(new_nev>nev)
	{
		cerr << ansi::red << "new_nev = " << new_nev << " is larger than nev = " << nev << ansi::reset << flush << endl;
		return; 
	};
	
	double* nA = new double[new_nev*new_nev];
	double* nB = new double[new_nev*new_nev];
	
	for(uint i=0; i<new_nev; i++)
		for(uint j=0; j<new_nev; j++)
		{
			nA[i*new_nev + j] = A[i*nev + j];
			nB[i*new_nev + j] = B[i*nev + j];
		};

	delete [] A; A = nA;
	delete [] B; B = nB;
	nev = new_nev;
}

void SupersymmetricSystem::SymmetryTransformation(uint ig, double* in, double* out)
{
	//Construct npsi after acting with a group element number ig
	//First we apply the coordinate transformation
	std::fill(out, out + 2*N, 0.0);
	//#pragma omp parallel for
	for(uint parity_sector=0; parity_sector<2; parity_sector ++)
		for(uint idx=0; idx<NB; idx++)
		{
			uint k1, k2, nk1, nk2, new_parity_sector = 0, nidx; double phase = 1.0;
				
			idx2cart(idx, parity_sector, k1, k2);
			phase = D4d_k1k2_fwd(ig, k1,  k2, nk1, nk2);
			nidx = cart2idx_full(nk1, nk2, new_parity_sector);
			out[(    new_parity_sector)*N      + nidx] += phase*D4d_E1[4*ig + 2*0 + 0]*in[parity_sector*N + idx];
		
			idx2cart(idx, 1-parity_sector, k1, k2);
			phase = D4d_k1k2_fwd(ig, k1,  k2, nk1, nk2);
			nidx = cart2idx_full(nk1, nk2, new_parity_sector);
			out[(    new_parity_sector)*N      + nidx] += phase*D4d_E1[4*ig + 2*0 + 1]*in[parity_sector*N + NB + idx];
		
			idx2cart(idx,   parity_sector, k1, k2);
			phase = D4d_k1k2_fwd(ig, k1,  k2, nk1, nk2);
			nidx = cart2idx_full(nk1, nk2, new_parity_sector);
			out[(1 - new_parity_sector)*N + NB + nidx] += phase*D4d_E1[4*ig + 2*1 + 0]*in[parity_sector*N +      idx];
		
			idx2cart(idx, 1-parity_sector, k1, k2);
			phase = D4d_k1k2_fwd(ig, k1,  k2, nk1, nk2);
			nidx = cart2idx_full(nk1, nk2, new_parity_sector);
			out[(1 - new_parity_sector)*N + NB + nidx] += phase*D4d_E1[4*ig + 2*1 + 1]*in[parity_sector*N + NB + idx];
		};
}

void HamiltonianSystem::ClassifyIrreps(string datadir, uint parity_sector, bool print_irrep_info)
{
	uint NF = 2*N;
	double* npsi   = new double[NG*NF]; double* tmp    = new double[NF];
	double* nschur = new double[NIRREPS];
	
	char irreps_fname[512];
   	sprintf(irreps_fname, "%s/eigensystems/irreps_%s_M%u_L%1.4lf_p%u.dat", datadir.c_str(), model_name, M, L, parity_sector);
   	cout << "\t\t" << ansi::white << "File to save irrep info: " << ansi::magenta << irreps_fname << ansi::reset << endl << flush; 
   	
	FILE* f = fopen(irreps_fname, "wb");
	if(f==NULL){cerr << ansi::red << "Cannot open the file " << irreps_fname << ansi::reset << endl << flush;};
	
	double max_tnorm = 0.0;
	
	for(uint ie=0; ie<nev; ie++)
	{
		for(uint ig=0; ig<NG; ig++) //loop over all group elements
			SymmetryTransformation(ig, psi + N*ie, npsi + ig*NF, parity_sector);
		
		for(uint irrep=0; irrep<NIRREPS; irrep++)
		{
			uint    D  = irrep_dims[irrep];
			t_irrep T  = irreps[irrep];
			double* wpsi = new double[D*D*NF];
			
			std::fill(wpsi, wpsi + D*D*NF, 0.0);
			
			for(uint ig=0; ig<NG; ig++)
				for(uint ij=0; ij<D*D; ij++)
					A_pluseq_bB(wpsi + NF*ij, T[D*D*ig + ij], npsi + ig*NF, NF);
					
			nschur[irrep] = norm(wpsi, D*D*NF);
			delete [] wpsi;
		};
		
		double maxnschur = nschur[0]; uint maxirrep = 0; double tnorm = nschur[0]*nschur[0];
		for(uint irrep=1; irrep<NIRREPS; irrep++)
		{
			tnorm += nschur[irrep]*nschur[irrep];
			if(nschur[irrep] > maxnschur){ maxnschur = nschur[irrep]; maxirrep = irrep; };
		};
		tnorm = sqrt(tnorm - maxnschur*maxnschur);
		
		max_tnorm = std::max(max_tnorm, tnorm);
				
		cint res = (tnorm < 1.0E-6 ? (cint)(+1)*(cint)(maxirrep+1) :  cint(-1)*(cint)(maxirrep+1) );
					
		fwrite(&res, sizeof(cint), 1, f); fflush(f);
		
		if(print_irrep_info)
		{
			char resstr[64];
			sprintf(resstr, "Unknown");
			if(std::abs(res)==1) sprintf(resstr, "A1%s", (res<0? "???" : ""));
			if(std::abs(res)==2) sprintf(resstr, "A2%s", (res<0? "???" : ""));
			if(std::abs(res)==3) sprintf(resstr, "B1%s", (res<0? "???" : ""));
			if(std::abs(res)==4) sprintf(resstr, "B2%s", (res<0? "???" : ""));
			if(std::abs(res)==5) sprintf(resstr, "E0%s", (res<0? "???" : ""));
			if(std::abs(res)==6) sprintf(resstr, "E1%s", (res<0? "???" : ""));
			if(std::abs(res)==7) sprintf(resstr, "E2%s", (res<0? "???" : ""));
			
			cout << "\t" << ansi::cyan << "E[" << ansi::green << std::setw(3) << ie << ansi::cyan << "] = " << ansi::yellow << std::fixed << std::setprecision(3) << std::setw(8) << E[ie] << ansi::cyan << " : " << ansi::magenta << resstr << ansi::reset << endl << flush;
		};
		
	};
	
	fclose(f);
	
	cout << ansi::green << "\t Max. tnorm = " << ansi::magenta << max_tnorm << ansi::reset << endl << endl << flush;

	delete [] npsi; delete [] tmp; delete [] nschur;
}

