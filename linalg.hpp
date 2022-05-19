#ifndef _LINALG_HPP_
#define _LINALG_HPP_

#include<iostream>
#include<algorithm>
#include<complex>
#include<cfloat>
#include<cstdint>
#include<cblas.h>
#include<lapacke.h>
#include <ctime>
#include <random>

using namespace std; 

typedef complex<double>       t_complex;
typedef unsigned int               uint;
typedef uint64_t                 uint64;
typedef signed char                cint;

#define SQR(_x) ((_x)*(_x))

#ifndef M_PI
constexpr auto M_PI = 3.14159265358979323846;
#endif

double      scalar_prod(double*    psi1, double*    psi2, uint n);
t_complex   scalar_prod(t_complex* psi1, t_complex* psi2, uint n);

double      norm(double*    psi, uint n);
double      norm(t_complex* psi, uint n);

double      norm_diff(double* psi1, double* psi2, uint n);
double      norm_diff(t_complex* psi1, t_complex* psi2, uint n);
double      unitarity_norm(t_complex* U, uint n);

void		A_pluseq_a_psi_dirprod_chi(double* A, double a, double* psi, double* chi, uint n); // A += a*psi*chi (outer product of psi and chi)

void        rescale(double* A, double a, uint n);
void		rescale(t_complex* A, t_complex a, uint n);

void        double2complex(t_complex* out, double* in, uint n);
t_complex*  double2complex(double* in, uint n);

double      unitarity_err(t_complex* U, uint n);

void  sort_eigensystem(           double* evals, double* evecs, uint n, uint nev);
void check_eigensystem(double* A, double* evals, double* evecs, uint n, uint nev, double* evec_err=NULL, double* ortho_err=NULL);

template <typename T>
double      orthonormality_norm(T* basis, uint nv, uint n)
{
	double err = 0.0;
	for(uint i=0; i<nv; i++)
		for(uint j=0; j<nv; j++)
		{
			T sp = scalar_prod(&(basis[i*n]), &(basis[j*n]), n);
			sp = (i==j? sp - 1.0 : sp);
			err += std::abs(sp)*std::abs(sp);
		};
	return sqrt(err);
}

template <typename T>
void       psi_eq_A_mult_chi(T* psi, T* A, T* chi, uint n) //Matrix-vector multiplication
{
	#ifdef _LINALG_OMP
	#pragma omp parallel for
	#endif
	for(uint i=0; i<n; i++)
	{
		psi[i] = 0.0;
		for(uint j=0; j<n; j++)
			psi[i] += A[i*n + j]*chi[j];
	};
};

template <typename T>
T*         psi_eq_A_mult_chi(T* A,   T* chi, uint n)
{
	T* r = new T[n];
	psi_eq_A_mult_chi(r, A, chi, n);
	return r;
};

template<typename T>
void    aA_plus_bB(T* out, T a, T* A, T b, T* B, uint n)
{
	#ifdef _LINALG_OMP
	#pragma omp parallel for
	#endif
	for(uint i=0; i<n; i++)
		out[i] = a*A[i] + b*B[i];
};

void  A_pluseq_bB(double* A, double b, double *B, uint n);
void  A_pluseq_bB(t_complex* A, t_complex b, t_complex *B, uint n);

void  A_eq_B_mult_C(double* A, double* B, double* C, uint n); //Matrix multiplication C = A*B
void  A_eq_B_mult_C(t_complex* A, t_complex* B, t_complex* C, uint n); //Matrix multiplication C = A*B

template <typename T>
T* A_eq_B_mult_C(T* B, T* C, uint n) //Matrix multiplication, returns A
{
	T* r = new T[n*n];
	A_eq_B_mult_C(r, B, C, n);
	return r;
};

template <typename T>
T tr(T* A, uint n)
{
	T r = 0.0;
	for(uint i=0; i<n; i++)
		r += A[i*n + i];
	return r;
}

template <typename T>
void identity_matrix(T* A, uint n)
{
	for(uint i=0; i<n; i++)
		for(uint j=0; j<n; j++)
			A[i*n + j] = (i==j? 1.0 : 0.0);
}

template <typename T>
T*  identity_matrix(uint n)
{
	T* r = new T[n*n];
	identity_matrix(r, n);
	return r;
}

void hermitian_conjugate(t_complex* out, t_complex* in, uint n);
void hermitian_conjugate(t_complex* A, uint n);
void           transpose(t_complex* A, uint n);
void           transpose(double*    A, uint n);

void  commutator(t_complex* C, t_complex* A, t_complex* B, uint n); //C = [A, B] = A*B - B*A

double*  commutator(double* A, double* B, uint n); //C = [A, B] = A*B - B*A

extern std::ranlux48 rng_engine;
extern std::normal_distribution<double> normdist;

void         rand_vec(double*    out, uint n);
void         rand_vec(t_complex* out, uint n);
double* rand_vec(uint n);

#endif
