#ifndef _D4D_HPP_
#define _D4D_HPP_

typedef unsigned int uint;

#define IS2 0.70710678118654752440084436210485

//This is the symmetry group of the supersymmetric system

const unsigned int D4d_Mult_Table[256] = 
{
   0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
   1,  2,  3,  4,  5,  6,  7,  0, 10, 15,  9,  8, 14, 11, 13, 12,
   2,  3,  4,  5,  6,  7,  0,  1,  9, 12, 15, 10, 13,  8, 11, 14,
   3,  4,  5,  6,  7,  0,  1,  2, 15, 14, 12,  9, 11, 10,  8, 13,
   4,  5,  6,  7,  0,  1,  2,  3, 12, 13, 14, 15,  8,  9, 10, 11,
   5,  6,  7,  0,  1,  2,  3,  4, 14, 11, 13, 12, 10, 15,  9,  8,
   6,  7,  0,  1,  2,  3,  4,  5, 13,  8, 11, 14,  9, 12, 15, 10,
   7,  0,  1,  2,  3,  4,  5,  6, 11, 10,  8, 13, 15, 14, 12,  9,
   8, 11, 13, 14, 12, 15,  9, 10,  0,  6,  7,  1,  4,  2,  3,  5,
   9, 10,  8, 11, 13, 14, 12, 15,  2,  0,  1,  3,  6,  4,  5,  7,
  10,  8, 11, 13, 14, 12, 15,  9,  1,  7,  0,  2,  5,  3,  4,  6,
  11, 13, 14, 12, 15,  9, 10,  8,  7,  5,  6,  0,  3,  1,  2,  4,
  12, 15,  9, 10,  8, 11, 13, 14,  4,  2,  3,  5,  0,  6,  7,  1,
  13, 14, 12, 15,  9, 10,  8, 11,  6,  4,  5,  7,  2,  0,  1,  3,
  14, 12, 15,  9, 10,  8, 11, 13,  5,  3,  4,  6,  1,  7,  0,  2,
  15,  9, 10,  8, 11, 13, 14, 12,  3,  1,  2,  4,  7,  5,  6,  0
};

//Group representation in the basis of harmonic oscillator functions
inline double D4d_k1k2_fwd(uint ig, uint k1, uint k2, uint &nk1, uint &nk2)
{
	double phase = 1.0;
	
	if(ig==0  || ig==4 ){nk1 = k1; nk2 = k2; phase = 1.0;                         };
	if(ig==1  || ig==5 ){nk1 = k2; nk2 = k1; phase = ((k1     )%2==0? 1.0 : -1.0);};
	if(ig==2  || ig==6 ){nk1 = k1; nk2 = k2; phase = ((k1 + k2)%2==0? 1.0 : -1.0);};
	if(ig==3  || ig==7 ){nk1 = k2; nk2 = k1; phase = ((k2     )%2==0? 1.0 : -1.0);};
	
	if(ig==8  || ig==12){nk1 = k1; nk2 = k2; phase = ((k1     )%2==0? 1.0 : -1.0);};
	if(ig==9  || ig==13){nk1 = k1; nk2 = k2; phase = ((k2     )%2==0? 1.0 : -1.0);};
	if(ig==10 || ig==14){nk1 = k2; nk2 = k1; phase = 1.0;						  }; 
	if(ig==11 || ig==15){nk1 = k2; nk2 = k1; phase = ((k1 + k2)%2==0? 1.0 : -1.0);};
	
	return phase;
}

inline double D4d_k1k2_bwd(uint ig, uint k1, uint k2, uint &nk1, uint &nk2)
{
	double phase = 1.0;
	
	if(ig==0  || ig==4 ){nk1 = k1; nk2 = k2; phase = 1.0;                         };
	if(ig==1  || ig==5 ){nk1 = k2; nk2 = k1; phase = ((k2     )%2==0? 1.0 : -1.0);};
	if(ig==2  || ig==6 ){nk1 = k1; nk2 = k2; phase = ((k1 + k2)%2==0? 1.0 : -1.0);};
	if(ig==3  || ig==7 ){nk1 = k2; nk2 = k1; phase = ((k1     )%2==0? 1.0 : -1.0);};
	
	if(ig==8  || ig==12){nk1 = k1; nk2 = k2; phase = ((k1     )%2==0? 1.0 : -1.0);};
	if(ig==9  || ig==13){nk1 = k1; nk2 = k2; phase = ((k2     )%2==0? 1.0 : -1.0);};
	if(ig==10 || ig==14){nk1 = k2; nk2 = k1; phase = 1.0;   					  };
	if(ig==11 || ig==15){nk1 = k2; nk2 = k1; phase = ((k1 + k2)%2==0? 1.0 : -1.0);};
	
	return phase;
}

//The four Abelian irreps of D4d
const double D4d_A1[16] = {+1, +1, +1, +1, +1, +1, +1, +1, +1, +1, +1, +1, +1, +1, +1, +1};
const double D4d_A2[16] = {+1, +1, +1, +1, +1, +1, +1, +1, -1, -1, -1, -1, -1, -1, -1, -1};
const double D4d_B1[16] = {+1, -1, +1, -1, +1, -1, +1, -1, +1, +1, -1, -1, +1, +1, -1, -1};
const double D4d_B2[16] = {+1, -1, +1, -1, +1, -1, +1, -1, -1, -1, +1, +1, -1, -1, +1, +1};

//Here are the three non-Abelian, 2d irreps of the D4d group
const double D4d_E0[64] = 
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
    -1,   0,
     0,   1
 ,
     1,   0,
     0,  -1
 ,
     0,   1,
     1,   0
 ,
     0,  -1,
    -1,   0
 ,
    -1,   0,
     0,   1
 ,
     1,   0,
     0,  -1
 ,
     0,   1,
     1,   0
 ,
     0,  -1,
    -1,   0
};

const double D4d_E1[64] = 
{
     1,   0,
     0,   1
 ,
   IS2,-IS2,
   IS2, IS2
 ,
     0,  -1,
     1,   0
 ,
  -IS2,-IS2,
   IS2,-IS2
 ,
    -1,   0,
     0,  -1
 ,
  -IS2, IS2,
  -IS2,-IS2
 ,
     0,   1,
    -1,   0
 ,
   IS2, IS2,
  -IS2, IS2
 ,
     1,   0,
     0,  -1
 ,
     0,   1,
     1,   0
 ,
   IS2, IS2,
   IS2,-IS2
 ,
   IS2,-IS2,
  -IS2,-IS2
 ,
    -1,   0,
     0,   1
 ,
     0,  -1,
    -1,   0
 ,
  -IS2,-IS2,
  -IS2, IS2
 ,
  -IS2, IS2,
   IS2, IS2
};

const double D4d_E2[64] = 
{
     1,   0,
     0,   1
 ,
  -IS2, IS2,
  -IS2,-IS2
 ,
     0,  -1,
     1,   0
 ,
   IS2, IS2,
  -IS2, IS2
 ,
    -1,   0,
     0,  -1
 ,
   IS2,-IS2,
   IS2, IS2
 ,
     0,   1,
    -1,   0
 ,
  -IS2,-IS2,
   IS2,-IS2
 ,
     1,   0,
     0,  -1
 ,
     0,   1,
     1,   0
 ,
  -IS2,-IS2,
  -IS2, IS2
 ,
  -IS2, IS2,
   IS2, IS2
 ,
    -1,   0,
     0,   1
 ,
     0,  -1,
    -1,   0
 ,
   IS2, IS2,
   IS2,-IS2
 ,
   IS2,-IS2,
  -IS2,-IS2
};

#endif
