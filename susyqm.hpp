#ifndef _SUSYQM_HPP_
#define _SUSYQM_HPP_

#include<iostream>
#include<iomanip>
#include<chrono>
#include<fstream>
#include<algorithm>
#include<vector>
#include<climits>

#include "arpackpp/arrssym.h"
#include "linalg.hpp"
#include "ansi_io.hpp"
#include "timing.hpp"

#include "C4v.hpp"
#include "D4d.hpp"

typedef const double* t_irrep; 

class HamiltonianSystem
{
	public:
		HamiltonianSystem(int M, double L);
		~HamiltonianSystem(){};
		int          M     = 0;
		double       L     = 1.0;
		int          nev   = 0;
		uint         NB    = 1;												//Size of the bosonic Hilbert space - for x2y2
		uint64       NB2   = 1;												//NB*NB
		uint         N     = 1;												//Size of the Hilbert space
		uint64       N2	   = 1;										    	//N*N
		uint         NG    = 1;												//Number of symmetry transformations
		uint         NIRREPS = 0;											//Number of irreps of the symmetry group
		char         model_name[16] = "";
		t_irrep*     irreps     = NULL;
		uint*        irrep_dims = NULL;
		const uint*  group_mult_table = NULL;
		/* Enumeration of states */
		uint64       cart2idx(uint k1, uint k2);
		uint64       cart2idx_full(uint k1, uint k2, uint &parity_sector);
		void         idx2cart(uint64 idx, uint parity_sector, uint &k1, uint &k2);
		/* Some checks for symmetry transformation */
		void                CheckSymmetryGroupImplementation();
		/* Hamiltonian definition */
		virtual void        H(double* in, uint parity_sector) = 0;						//The matrix in contains N*N Hamiltonian matrix
		virtual void        H(double* in, double* out, uint parity_sector) = 0;			//Hamiltonian, in and out are real since the Hamiltonian is real
		void                H(double* in, double* out);								//Acts on the combined Hilbert space of all parity sectors, for testing pur
		virtual double*     H(uint parity_sector);										//Returns N*N matrix
		/* Eigenstate algorithms */
		        void FindAllEigenstates(uint parity_sector);
		         int FindLowestEigenstates(uint nev, uint parity_sector, double precision=1.0E-10, uint ncv=0);
		        void CheckEigensystem(uint parity_sector, double* eigensys_err=NULL, double* orthogonality_err=NULL);
		        void CheckHermiticity(uint parity_sector, double* err=NULL);
		        void ReorthogonalizeEigenstates(double precision);
		        void ClassifyIrreps(string datadir, uint parity_sector, bool print_irrep_info = false);
		virtual void SymmetryTransformation(uint ig, double* in, double* out) = 0;
	            void SymmetryTransformation(uint ig, double* in, double* out, uint parity_sector);
		/* Placeholders for the eigenspectrum */						
		double*     E    = NULL;									//Eigenenergies
		double*     psi  = NULL;									//The real-value matrix of eigenvectors, H = O.E.O^T, E is treated as diagonal matrix
		//Saving/reading eigensystems
		uint         read_eigensystem(string datadir, uint parity_sector, bool  read_evecs=true);
		void        write_eigensystem(string datadir, uint parity_sector, bool write_evecs=true);
		//Basic matrix elements in the 2D harmonic oscillator basis with length parameter L ( <0|x^2|0>=L^2/2 )
		double       x(uint k, uint l);
		double      ip(uint k, uint l);
		double     xsq(uint k, uint l);
		double     psq(uint k, uint l);
		//Operators to be used in OTOCs
		virtual void  x2(double* in, double* out, uint parity_sector) = 0;
		virtual void ip2(double* in, double* out, uint parity_sector) = 0;
		//OTOC structures
		void     init_OTOCS(uint parity_sector); //iA, iB are codes of operators: 0 = x1, 1 = x2, 2 = p1, 3 = p2
		void    write_OTOC_data(string datadir, uint parity_sector);
		int      read_OTOC_data(string datadir, uint parity_sector, uint fnev=0);
		void   resize_OTOC_data(uint new_nev);
		double* A = NULL; //Matrix elements of operator A
		double* B = NULL; //Matrix elements of operator B
		void       OTOC_regularized(double beta, double t, double* r1=NULL, double* r2=NULL, double* r3=NULL);
		void       thermal_vevs_regularized(double beta, double& aA, double& aB, double& Z);
		void       thermal_vevs(            double beta, double& aA, double& aB, double& Z);
};

class BosonicSystem : public HamiltonianSystem
{
	public:
		BosonicSystem(int M, double L, bool free_system=false);
		~BosonicSystem(){};
		bool free_system = false;
		/* Defining the Hamiltonian */
		void        H(double* in, uint parity_sector);						//The matrix in contains N*N Hamiltonian matrix
		void        H(double* in, double* out, uint parity_sector);			//Hamiltonian, in and out are real since the Hamiltonian is real
		/* Defining the operators to be used in OTOCs - and also the SUSY Hamiltonian */
		void pluseq_x1(  double* in, double* out, uint parity_sector_in, double c = 1.0);
		void pluseq_x2(  double* in, double* out, uint parity_sector, double c = 1.0);
		void pluseq_ip2( double* in, double* out, uint parity_sector, double c = 1.0);
		virtual void x2( double* in, double* out, uint parity_sector);
		virtual void ip2(double* in, double* out, uint parity_sector);
		void SymmetryTransformation(uint ig, double* in, double* out);
};

class SupersymmetricSystem : public BosonicSystem
{
	public:
		SupersymmetricSystem(int M, double L);
		~SupersymmetricSystem(){};
		// Defining the Hamiltonian 
		void        H(double* in, uint parity_sector);							//The matrix in contains N*N Hamiltonian matrix
		void        H(double* in, double* out, uint parity_sector);			//Hamiltonian, in and out are real since the Hamiltonian is real
		// Defining the operators to be used in OTOCs
		void  x2(double* in, double* out, uint parity_sector);
		void ip2(double* in, double* out, uint parity_sector);
		void SymmetryTransformation(uint ig, double* in, double* out);
};

#endif
