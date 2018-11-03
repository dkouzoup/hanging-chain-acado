
// Linear algebra for the generalized dual fast gradient method

#include "include/gdfgm_algebra.h"

#include <math.h>
#include <stdio.h>

// initialize vector to 0.0
void initVector(real_t *vec, int dim)
{
    
#ifndef USE_EXTERNAL_LIBRARIES
    int i;
    
    for(i=0;i<dim;i++)
    {
        vec[i] = 0.0;
    }
#else
    cblas_dscal(dim,0.0,vec,1);
#endif
   
}

void copyVector(const real_t *vsrc, real_t *vdest, const int dim)
{

#ifndef USE_EXTERNAL_LIBRARIES
    int i;

    for(i=0;i<dim;i++)
    {
        vdest[i] = vsrc[i];
    }
#else
    cblas_dcopy(dim, vsrc, 1, vdest, 1); // currently not assigning val, but 0 (TODO)
#endif    

}


// transpose matrix (for compatibility between row and column major)
void transposeMatrix(real_t *M, int nRows, int nCols, int major)
{
    // major = 1 if input is saved in row major, 2 for column major
    // result overwrittes M
    
    real_t Mtemp[nRows*nCols];
    
    int r,c,i;
    
    if (major == 1)
    {
        for(r=0;r<nRows;r++)
        {
            for(c=0;c<nCols;c++)
            {
                Mtemp[c*nRows+r] = M[r*nCols+c];
            }
        }
    }
    else if (major == 2)
    {
        for(c=0;c<nCols;c++)
        {
            for(r=0;r<nRows;r++)
            {
                Mtemp[r*nCols+c] = M[c*nRows+r];
            }
        }
    }
    
    for(i=1;i<nRows*nCols;i++)
    {
        M[i]=Mtemp[i];
    }
    
}

/* column major matrix-vector product
void matrixVectorProductColMajor(const real_t *M, const real_t *v, real_t *res, const int nRows, const int nCols)
{
#ifndef USE_EXTERNAL_LIBRARIES
    
    int i,j;
    
    for(i=0;i<nRows;i++)
    {
        res[i] = M[i]*v[0];
    }
    for(j=1;j<nCols;j++)
    {
        for(i=0;i<nRows;i++)
        {
            res[i] += M[j*nRows+i]*v[j];
        }
        
    }
    
#else
    cblas_dgemv(CblasColMajor, CblasNoTrans, nRows, nCols, 1.0, M, nRows, v, 1, 0.0, res, 1);
#endif

}*/


/* * * * * * * * * CBLAS replacement functions * * * * * * * * * * * * * */

void myblas_dcopy(const int N, const real_t *X, const int incX, real_t *Y, const int incY)
{

#ifndef USE_EXTERNAL_LIBRARIES
    
    int i;
    
    for(i=0;i<N;i++)
    {
        Y[i*incY] = X[i*incX]; 
    }
    
#else
    cblas_dcopy(N, X, incX, Y, incY);
#endif 

}

// X = alpha*X  ( X must be of dimension N*incX - function not used)
void myblas_dscal(const int N, const real_t alpha, real_t *X, const int incX)
{

#ifndef USE_EXTERNAL_LIBRARIES
    
    int i;

    // vector X must be initialized first!
    
    for(i=0;i<N;i++)
    {
        X[i*incX] = alpha*X[i*incX]; 
    }
    
#else
    cblas_dscal(N, alpha, X, incX);
#endif 

}


// Y = alpha*X + Y
void myblas_daxpy(const int N, const real_t alpha, const real_t *X, const int incX, real_t *Y, const int incY)
{

#ifndef USE_EXTERNAL_LIBRARIES

    int i;

    for(i=0;i<N;i++)
    {
        Y[i*incY] += alpha*X[i*incX];
    }

#else
    cblas_daxpy(N, alpha, X, incX, Y, incY);
#endif

}


// Y = alpha*(A*X) + beta*Y (WILL NOT BE USED - too complicated for no reason - INCOMPLETE)
#ifdef USE_EXTERNAL_LIBRARIES 
void myblas_dgemv (const int Order, const int TransA, 
        const int M, const int N, const real_t alpha, const real_t *A, 
        const int lda, const real_t *X, const int incX, const real_t beta, 
        real_t *Y, const int incY)
#else
void myblas_dgemv (const int Order, const int TransA, 
        const int M, const int N, const real_t alpha, const real_t *A, 
        const int lda, const real_t *X, const int incX, const real_t beta, 
        real_t *Y, const int incY)  
#endif        
{
#ifndef USE_EXTERNAL_LIBRARIES
  
     int i,j;
  
     if (Order == CblasColMajor) 
     {
	  if (TransA == CblasTrans) //inefficient, better perform mult directly..
	  {
	       for (i = 0; i<N; i++)
		    Y[i] *= beta/alpha;

	       for (i = 0; i < N; i++)
		    for (j = 0; j < M; j++)
			 Y[i] += A[i*M+j]*X[j];
	       
	       for (i = 0; i < N; i++)
		    Y[i] *= alpha;
	  } 
	  else
	  {
	       // in first column scale existing Y with (alpha)/(beta)
	       for(i=0;i<M;i++)
	       {
		    Y[i] *= beta/alpha;
		    Y[i] += A[i]*X[0];
	       }
	       // intermediate 
	       for(j=1;j<N;j++)
	       {
		    //Y[i] *= beta/alpha;
		    for(i=0;i<M;i++)
		    {
			 Y[i] += A[j*M+i]*X[j];
		    }
		    
	       }
	       // in last column multiply with alpha
	       for(i=0;i<M;i++)
	       {
		    Y[i] *= alpha;
	       }
	  }
	  
     }
     else if (Order == CblasRowMajor)
     {
	  if (TransA == CblasTrans)
	  {
	       //transposeMatrix(A, M, N, 1);
	  }
	  
	  // TODO: complete multiplication
	  
     }
#else
     cblas_dgemv(Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
#endif    
}

// Cholsesky factorization (Code taken from qpOASES_e 3.0)
void myblas_dpotrf_(const char *uplow, const int *_n, real_t *a,
		    const int *_lda, int *info)
{
#if 0 //def USE_EXTERNAL_LIBRARIES
  // TODO: this seg faults on linux (not on mac), disabled
  dpotrf_(uplow, _n, a, _lda, info);
  
#else 
  real_t sum;
  int i, j, k;
  int n = (int)(*_n);
  int lda = (int)(*_lda);
  real_t tmp;
 
  for( i=0; i<n; ++i )
    {
      /* j == i */
      sum = a[i + lda*i];
      
      for( k=(i-1); k>=0; --k )
	sum -= a[i+lda*k] * a[i+lda*k];
      
      if ( sum > 0.0 )
	a[i+lda*i] = sqrt( sum );
      else
	{
	  a[0] = sum; /* tunnel negative diagonal element to caller */
	  if (info != 0)
	    *info = (int)i+1;
	  return;
	}
      
      for( j=(i+1); j<n; ++j )
	{
	  sum = a[i*lda + j];
	  
	  for( k=(i-1); k>=0; --k )
	    sum -= a[i+lda*k] * a[j+lda*k];
	  
	  a[j+lda*i] = sum / a[i+lda*i];
	}
    }
  if (info != 0)
    *info = 0;

#endif
}

// 2-norm of a vector
real_t myblas_dnrm2(const int N, const real_t *X, const int incX)
{
#ifdef USE_EXTERNAL_LIBRARIES
     
     return cblas_dnrm2(N, X, incX);
     
#else
     
     int i;
     real_t sum = 0.0;

     for (i = 0; i < N; i++)
	  sum += X[i]*X[i];
  
     return sqrt (sum);
#endif
     
}

// Matrix matrix multiplication
void myblas_dgemm(const int Order, const int TransA, const int TransB, const int M, 
		  const int N, const int K, const real_t alpha, const real_t *A, 
		  const int lda, const real_t *B, const int ldb, const real_t beta,
		  real_t *C , const int ldc)
{
#ifdef USE_EXTERNAL_LIBRARIES    

     cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta,
		 C, ldc);
#else
     int i, j, k;
     // NOTE: this is tailored to our application only

     for (i = 0; i < M; i++)
	  for (j = 0; j < N; j++)
	       for (k = 0; k < K; k++)
		    C[i + N*j] -= A[i + k*M] * B[k*K + j]; 
		    
#endif
}

// Matrix substitution 
void myblas_dtrsm(const int Order, const int Side, const int Uplo, const int TransA,
		  const int Diag, const int M, const int N, const real_t alpha, 
		  const real_t *A, const int lda, real_t *B, const int ldb)
{
#ifdef USE_EXTERNAL_LIBRARIES    
     
     cblas_dtrsm(Order, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);

#else
     // NOTE: this is tailored to our application only
     
     int row, col, n;
  
     for (n = 0; n < N; n++)
     {
	  for (col = 0; col < M; col++)
	  {
	       for (row = 0; row < col; row++)
	       {
		    B[col*N + n] -= A[col + row*M] * B[row*N + n];    
	       }
	       B[col*N + n] = B[col*N + n] / A[col + col*M];
	  }
	  
     }
    
#endif
}
