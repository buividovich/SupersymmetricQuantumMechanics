#ifndef _C4V_HPP_
#define _C4V_HPP_

typedef unsigned int uint;

//This is the symmetry group of the bosonic system

const unsigned int C4v_Mult_Table[64] = 
{
  0, 1, 2, 3, 4, 5, 6, 7 ,
  1, 2, 3, 0, 7, 4, 5, 6 ,
  2, 3, 0, 1, 6, 7, 4, 5 ,
  3, 0, 1, 2, 5, 6, 7, 4 ,
  4, 5, 6, 7, 0, 1, 2, 3 ,
  5, 6, 7, 4, 3, 0, 1, 2 ,
  6, 7, 4, 5, 2, 3, 0, 1 ,
  7, 4, 5, 6, 1, 2, 3, 0 
};

//Group representation in the basis of harmonic oscillator functions
inline double C4v_k1k2_fwd(uint ig, uint k1, uint k2, uint &nk1, uint &nk2)
{
	double phase = 1.0;
	
	if(ig==0){nk1 = k1; nk2 = k2; phase = 1.0;                         };
	if(ig==1){nk1 = k2; nk2 = k1; phase = ((k1     )%2==0? 1.0 : -1.0);};
	if(ig==2){nk1 = k1; nk2 = k2; phase = ((k1 + k2)%2==0? 1.0 : -1.0);};
	if(ig==3){nk1 = k2; nk2 = k1; phase = ((k2     )%2==0? 1.0 : -1.0);};
	
	if(ig==4){nk1 = k1; nk2 = k2; phase = ((k2     )%2==0? 1.0 : -1.0);};
	if(ig==5){nk1 = k2; nk2 = k1; phase = 1.0;                         };
	if(ig==6){nk1 = k1; nk2 = k2; phase = ((k1     )%2==0? 1.0 : -1.0);};
	if(ig==7){nk1 = k2; nk2 = k1; phase = ((k1 + k2)%2==0? 1.0 : -1.0);};
	
	return phase;
}

inline double C4v_k1k2_bwd(uint ig, uint k1, uint k2, uint &nk1, uint &nk2)
{
	double phase = 1.0;
	
	if(ig==0 ){nk1 = k1; nk2 = k2; phase = 1.0;                         };
	if(ig==1 ){nk1 = k2; nk2 = k1; phase = ((k2     )%2==0? 1.0 : -1.0);};
	if(ig==2 ){nk1 = k1; nk2 = k2; phase = ((k1 + k2)%2==0? 1.0 : -1.0);};
	if(ig==3 ){nk1 = k2; nk2 = k1; phase = ((k1     )%2==0? 1.0 : -1.0);};
	
	if(ig==4 ){nk1 = k1; nk2 = k2; phase = ((k2     )%2==0? 1.0 : -1.0);};
	if(ig==5 ){nk1 = k2; nk2 = k1; phase = 1.0;                         };
	if(ig==6 ){nk1 = k1; nk2 = k2; phase = ((k1     )%2==0? 1.0 : -1.0);};
	if(ig==7 ){nk1 = k2; nk2 = k1; phase = ((k1 + k2)%2==0? 1.0 : -1.0);};
	
	return phase;
}

//Here are the Abelian irreps of C4v

const double C4v_A1[8] = {+1, +1,  +1, +1, +1, +1, +1, +1};
const double C4v_A2[8] = {+1, +1,  +1, +1, -1, -1, -1, -1};
const double C4v_B1[8] = {+1, -1,  +1, -1, +1, -1, +1, -1};
const double C4v_B2[8] = {+1, -1,  +1, -1, -1, +1, -1, +1};

//Here is the only non-Abelian, 2d irrep of the C4v group

const double C4v_E0[32] = 
{
     1,   0,
     0,   1
 ,
     0,   1,
    -1,   0
 ,
    -1,   0,
     0,  -1
 ,
     0,  -1,
     1,   0
 ,
     1,   0,
     0,  -1
 ,
     0,   1,
     1,   0
 ,
    -1,   0,
     0,   1
 ,
     0,  -1,
    -1,   0
};

#endif

 //Warm-up - let's test the Schur's orthogonality relations
 //Loop over all combinations of fundamental indices
 /*for(uint aij=0; aij<16; aij++)
 {
 	uint ij = aij;
 	uint j2 = ij%2; ij /= 2;
 	uint i2 = ij%2; ij /= 2;
 	uint j1 = ij%2; ij /= 2;
 	uint i1 = ij%2; ij /= 2;
 	
 	double S = 0.0;
 	for(uint ig=0; ig<16; ig++)
 		S += D4d_E1[ig][i1][j1]*D4d_E2[ig][i2][j2];
 		
 	if(std::abs(S)>1.0E-10)
 		cout << ansi::yellow << i1 << j1 << i2 << j2 << ": " << ansi::cyan << S << ansi::reset << endl << flush;
 };*/
