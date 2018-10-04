

#ifndef GDFGM_DATA_H
#define GDFGM_DATA_H

// inputs
real_t gdfgm_Am[NI*NX*NX];
real_t gdfgm_Bm[NI*NX*NU];
real_t gdfgm_beq[(NI+1)*NX];

real_t gdfgm_Hinv[NP]; // H has a diagonal structure 
real_t gdfgm_f[NP];

real_t gdfgm_lb[NP];
real_t gdfgm_ub[NP];


// intermediate results
real_t gdfgm_Ld[(NI+1)*NX*NX];
real_t gdfgm_Ll[NI*NX*NX];
real_t gdfgm_Lb[2*NX*ND]; 

// a-priori information
real_t gdfgm_optSol[NP];
real_t gdfgm_lam[ND]; // also result

// results
real_t gdfgm_sol[NP];
real_t gdfgm_fval;
real_t gdfgm_mu[2*NP];

int gdfgm_nIter;



#endif