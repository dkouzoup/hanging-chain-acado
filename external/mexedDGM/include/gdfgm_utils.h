

#ifndef GDFGM_UTILS_H
#define GDFGM_UTILS_H

#include "gdfgm_dimensions.h"

/** Definition of the floating point data type. */
#ifndef REAL_T
#define REAL_T
#ifdef __USE_SINGLE_PRECISION__
  typedef float real_t;
#else
  typedef double real_t;
#endif	/* __USE_SINGLE_PRECISION__ */
#endif
  
#ifdef USE_EXTERNAL_LIBRARIES
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
extern void dpotrf_( const char* uplo, const int* n, real_t* a, const int* lda, 
		     int *info );
#endif
#endif

/* Function declarations */

void buildMatrixL(real_t *Hinv, real_t *Am, real_t *Bm,
        real_t *Ld, real_t *Ll );

void mainDiag(real_t *res, real_t *A, real_t *B, real_t *Q, real_t *R);

void offDiag(int sign, int dim, int ld, real_t *res, real_t *M, real_t *D);

void blockCholeskyFactorization(real_t *Ld, real_t *Ll);

void bandedFactorization(real_t *Ld, real_t *Ll, real_t *Lb);

//void fullFactorization(real_t *Ld, real_t *Ll, real_t *L);

void AeqTimesVector(real_t *Am, real_t *Bm, real_t *z, real_t *res);

void AeqTransposeTimesVector(real_t *Am, real_t *Bm, real_t *y, real_t *res);

void LmTimesVector(real_t *Ld, real_t *Ll, real_t *z, real_t *res);

void blockForwardSubstitution(real_t *Ld, real_t *Ll, real_t *x);

void blockBackwardSubstitution(real_t *Ld, real_t *Ll, real_t *x);

#endif
