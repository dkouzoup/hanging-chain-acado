

#ifndef GDFGM_CORE_H
#define GDFGM_CORE_H

void runFGMiterations(real_t Hinv[NP], real_t f[NP], 
		      real_t Am[NI*NX*NX], real_t Bm[NI*NX*NU], real_t beq[(NI+1)*NX], 
		      real_t ub[NP], real_t lb[NP], 
		      real_t Ld[(NI+1)*NX*NX], real_t Ll[NI*NX*NX], 
		      real_t z[NP], real_t zunc[NP], real_t lambda[ND], real_t yInit[ND],
		      int *nIter, gdfgm_options opt);

#endif