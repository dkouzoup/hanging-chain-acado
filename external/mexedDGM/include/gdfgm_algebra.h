
#ifndef GDFGM_ALGEBRA_H
#define GDFGM_ALGEBRA_H

#ifdef USE_EXTERNAL_LIBRARIES
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
extern void dpotrf_( const char* uplo, const int* n, real_t* a, const int* lda, 
		     int *info );
#endif
#endif

/** Definition of the floating point data type. */
#ifndef REAL_T
#define REAL_T
#ifdef __USE_SINGLE_PRECISION__
  typedef float real_t;
#else
  typedef double real_t;
#endif	/* __USE_SINGLE_PRECISION__ */
#endif

      
/* Function declarations */

void initVector(real_t *vec, int dim);

void copyVector(const real_t *vsrc, real_t *vdest, const int dim);

void transposeMatrix(real_t *M, int nRows, int nCols, int major);

void matrixVectorProductColMajor(const real_t *M, const real_t *v, real_t *res, const int nRows, const int nCols);


//myblas - replacing cblas (for easier changes in gdfgm algorithm)

void myblas_dcopy(const int N, const real_t *X, const int incX, real_t *Y, const int incY);

void myblas_dscal(const int N, const real_t alpha, real_t *X, const int incX);

void myblas_daxpy(const int N, const real_t alpha, const real_t *X, const int incX, real_t *Y, const int incY);

real_t myblas_dnrm2(const int N, const real_t *X, const int incX);

//void myblas_dpotrf_(const char *uplow, const int *_n, real_t *a, const int *_lda, int *info);
void myblas_dpotrf_(char *uplow, int *_n, real_t *a, int *_lda, int *info);

void myblas_dgemm(const int Order, const int TransA, const int TransB, const int M, 
		  const int N, const int K, const real_t alpha, const real_t *A, 
		  const int lda, const real_t *B, const int ldb, const real_t beta,
		  real_t *C , const int ldc);

void myblas_dtrsm(const int Order, const int Side, const int Uplo, const int TransA,
                 const int Diag, const int M, const int N, const real_t alpha, 
		 const real_t *A, const int lda, real_t *B, const int ldb);

#ifndef USE_EXTERNAL_LIBRARIES
enum CBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO      {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG      {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE      {CblasLeft=141, CblasRight=142};
#endif

#ifdef USE_EXTERNAL_LIBRARIES
void myblas_dgemv (const int Order, const int TransA, 
        const int M, const int N, const real_t alpha, const real_t *A, 
        const int lda, const real_t *X, const int incX, const real_t beta, 
        real_t *Y, const int incY);
#else
void myblas_dgemv (const int Order, const int TransA, 
        const int M, const int N, const real_t alpha, const real_t *A, 
        const int lda, const real_t *X, const int incX, const real_t beta, 
        real_t *Y, const int incY) ;
#endif

#endif