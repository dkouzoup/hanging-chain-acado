
#include <stdio.h>
#include <math.h>

#include "include/gdfgm_utils.h"
#include "include/gdfgm_dimensions.h"
#include "include/gdfgm_algebra.h"

void buildMatrixL(real_t *Hinv, real_t *Am, real_t *Bm, 
        real_t *Ld, real_t *Ll)
{
    int i,j,n; 
    
    for(i=0;i<NX;i++)
    {
        // Fill in first NX diag elements with inv(Q_0)
        Ld[i*(NX+1)] = Hinv[i];
    }
    
    // Fill in second block and corresponding offdiagonal block    
    mainDiag(&Ld[NX*NX], &Am[0], &Bm[0], &Hinv[0], &Hinv[NX]);
    
    for(i=0;i<NX;i++)
    {
        Ld[NX*NX+i*NX+i] += Hinv[(NX+NU)+i];
    }
    
    // Add A_n*Qinv block to the left 
    // *rest of offdiagonal blocks have a minus sign
    offDiag(1, NX, NX, &Ll[0], &Am[0], &Hinv[0]);
    
    // Fill in from third block and on
    for(n=1;n<NI;n++)
    {
        mainDiag(&Ld[(n+1)*NX*NX], &Am[n*NX*NX], &Bm[n*NX*NU], 
                &Hinv[n*(NX+NU)], &Hinv[n*(NX+NU)+NX]);
        
        for(i=0;i<NX;i++)
        {
            Ld[(n+1)*NX*NX+i*NX+i] += Hinv[(n+1)*(NX+NU)+i];
        }
        
        // add -A_n*Qinv block to the left
        offDiag(-1, NX, NX, &Ll[n*NX*NX], &Am[n*NX*NX], &Hinv[n*(NX+NU)]);
    }
     
}


void mainDiag(real_t *res, real_t *A, real_t *B, real_t *Q, real_t *R)
{
    int i,j,k;
    
    for(j=0;j<NX;j++) // columns
        
    {
        for(i=j;i<NX;i++) // rows
            
        {
            for(k=0;k<NX;k++)
            {
                res[j*NX+i] += A[k*NX+i]*A[k*NX+j]*Q[k];
            }
            for(k=0;k<NU;k++)
            {
                res[j*NX+i] += B[k*NX+i]*B[k*NX+j]*R[k];
            }
        }
    }
}


void offDiag(int sign, int dim, int ld, real_t *res, real_t *M, real_t *D)
{
    int i,j;
    
    if(sign == 1)
    {
        for(j=0;j<dim;j++)
        {
            for(i=0;i<dim;i++)
            {
                res[j*ld+i] = M[j*ld+i]*D[j];
            }
        }
    }
    else if(sign == -1)
    {
        for(j=0;j<dim;j++)
        {
            for(i=0;i<dim;i++)
            {
                res[j*ld+i] = -M[j*ld+i]*D[j];
            }
        }
    }
    
}

void AeqTimesVector(real_t *Am, real_t *Bm, real_t *z, real_t *res)
{
    int i,n;
        
    // First  NX elements stay the same due to identity matrix
    for(i=0;i<NX;i++)
    {
        res[i] = z[i];
    }
	   
    for (i = NX; i < ND; i++)
      res[i] = 0.0;
    
    for(n=0;n<NI;n++)
    {
        //A*x{n}
        myblas_dgemv(CblasColMajor, CblasNoTrans, NX, NX, 1.0, 
                &Am[n*NX*NX], NX, &z[n*(NX+NU)], 
                1, 1.0, &res[(n+1)*NX], 1);
        
        //+B*u_{n}
        myblas_dgemv(CblasColMajor, CblasNoTrans, NX, NU, 1.0,
                &Bm[n*NX*NU], NX, &z[n*(NX+NU)+NX],
                1, 1.0, &res[(n+1)*NX], 1);
        
        //-x_{n+1}
        myblas_daxpy(NX, -1.0, &z[(n+1)*(NX+NU)] , 1, &res[(n+1)*NX], 1);
    }
  
}

void AeqTransposeTimesVector(real_t *Am, real_t *Bm, 
			     real_t *y, real_t *res)
{
    int i,n;
    
    // Initialize result
    for (i = 0; i < NP; i++)
	 res[i] = 0.0;
    
    // Fill in first NX+NU elements
    myblas_dgemv(CblasColMajor, CblasTrans, NX, NX, 1.0,
            &Am[0], NX, &y[NX], 1, 1.0, &res[0], 1);
    
    myblas_daxpy(NX, 1.0, &y[0], 1, &res[0], 1);
    
    myblas_dgemv(CblasColMajor, CblasTrans, NX, NU, 1.0,
            &Bm[0], NX, &y[NX], 1, 1.0, &res[NX], 1);
    
    // Fill in middle elements
    for(n=1;n<NI;n++)
    {
        myblas_dgemv(CblasColMajor, CblasTrans, NX, NX, 1.0, &Am[n*NX*NX], 
                NX, &y[(n+1)*NX], 1, 1.0, &res[n*(NX+NU)], 1);
        
        myblas_daxpy(NX, -1.0, &y[n*NX], 1, &res[n*(NX+NU)], 1);
       
        myblas_dgemv(CblasColMajor, CblasTrans, NX, NU, 1.0, &Bm[n*NX*NU], 
                NX, &y[(n+1)*NX], 1, 1.0, &res[n*(NX+NU)+NX], 1);
    }
    
    // Fill in last NX elements
    myblas_daxpy(NX, -1.0, &y[n*NX], 1, &res[n*(NX+NU)], 1);
    
}


void blockCholeskyFactorization(real_t *Ld, real_t *Ll)
{
    char low = 'L';
    int nblk = NX;
    int info;
        
    int i,n;
    
    // first block is just square root of first diagonal block
    for(i=0;i<NX;i++)
    {
        Ld[i*NX+i] = sqrt(Ld[i*NX+i]);
        
    }
    
    for(n=1;n<=NI;n++)
    {
        // Matrix substitution
        myblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasTrans, 
                CblasNonUnit, NX, NX, 1.0, &Ld[(n-1)*NX*NX], NX, 
                &Ll[(n-1)*NX*NX], NX);
        
        // Substitution of diagonal block
        // *this operation adds also false data to upper triangular part 
        //  that should NOT be referenced      
        myblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, NX, NX, NX,
                -1.0, &Ll[(n-1)*NX*NX], NX, &Ll[(n-1)*NX*NX], NX,
                1.0, &Ld[n*NX*NX], NX); 
        
        // Cholesky decomposition
        myblas_dpotrf_(&low,&nblk,&Ld[n*NX*NX],&nblk,&info);
    }
     
}

void blockForwardSubstitution(real_t *Ld, real_t *Ll, real_t *x)
{
     int row, col, n;

     // Substitute first diagonal block   
     for (row = 0; row < NX; row++)
     {
	  for (col = 0; col <= row-1; col++)
	  {
	       x[row] -= Ld[row + col*NX] * x[col];    
	  }
	  x[row] = x[row] / Ld[row + row*NX];
     }
     
     // Substitute remaining diagonal and off-diagonal blocks
     for (n = 1; n <= NI; n++)
     {
	  for (row = 0; row < NX; row++)
	  {
	       for (col = 0; col < NX; col++)
	       {
		    x[n*NX + row] -= Ll[(n-1)*NX*NX + row + col*NX] * x[(n-1)*NX + col];
	       }
	       for (col = 0; col <= row-1; col++)
	       {
		    x[n*NX + row] -= Ld[n*NX*NX + row + col*NX] * x[n*NX + col];
	       }
	       x[n*NX + row] = x[n*NX + row] / Ld[n*NX*NX + row + row*NX];
	  }
     }
}

void blockBackwardSubstitution(real_t *Ld, real_t *Ll, real_t *x)
{
     int row, col, n;

     // Substitute last diagonal block   
     for (row = NX-1; row >= 0; row--)
     {
	  for (col = NX-1; col > row; col--)
	  {
	       x[NI*NX + row] -= Ld[NI*NX*NX + col + row*NX] * x[NI*NX + col];    
	  }
	  x[NI*NX + row] = x[NI*NX + row] / Ld[NI*NX*NX + row + row*NX];
     }
     
     // Substitute remaining diagonal and off-diagonal blocks
     for (n = NI-1; n >= 0; n--)
     {
	  for (row = NX-1; row >= 0; row--)
	  {
	       for (col = NX-1; col >= 0; col--)
	       {
		    x[n*NX + row] -= Ll[(n)*NX*NX + col + row*NX] * x[(n+1)*NX + col];
	       }
	       for (col = NX-1; col > row; col--)
	       {
		    x[n*NX + row] -= Ld[n*NX*NX + col + row*NX] * x[n*NX + col];
	       }
	       x[n*NX + row] = x[n*NX + row] / Ld[n*NX*NX + row + row*NX];
	  }
     }
}



void bandedFactorization(real_t *Ld, real_t *Ll, real_t *Lb)
{
    int i,n,k,m;
    
    
    for(i=0;i<NX;i++) // fill first NX diagonals
    {
        m = i; // pointer to current element of result
        
        for(n=0;n<=NI;n++)
        {
            for(k=0;k<NX-i;k++) // fill in elements from diag blocks
            {
                Lb[m] = Ld[n*NX*NX+i+k*NX+k];
                m = m + 2*NX;
            }
            for(k=0;k<i;k++) // fill in elements from off-diag blocks
            {
                Lb[m] = Ll[(n+1)*NX*NX-i*NX+k*NX+k];
                m = m + 2*NX;
            }
            
        }
        
    }
    
    for(i=0;i<NX;i++) // fill in next NX diagonals
    {
        m = i+NX;
        
        for(n=0;n<NI;n++)
        {
            for(k=0;k<NX-i;k++)
            {
                Lb[m] = Ll[n*NX*NX+i+k*NX+k];
                m = m + 2*NX;
            }
            
            m = m + i*2*NX;
            
        }
        
    }
    
}


/*void fullFactorization(real_t *Ld, real_t *Ll, real_t *L)
{
    int n,i,j;
    
    for(j=0;j<NX;j++)
    {
        for(i=0;i<NX;i++)
        {
            L[j*ND+i] = Ld[j*NX+i];
        }
    }
    
    for(n=0;n<NI;n++)
    {
        for(j=0;j<NX;j++)
        {
            for(i=0;i<NX;i++)
            {
                L[(n+1)*NX*(ND+1)+j*ND+i] = Ld[(n+1)*NX*NX + j*NX+i];
            }
        }
        for(j=0;j<NX;j++)
        {
            for(i=0;i<NX;i++)
            {
                L[n*NX*(ND+1)+NX+j*ND+i] = Ll[n*NX*NX + j*NX+i];
            }
        }
        
        
    }
    
}
*/
/*
void LmTimesVector(real_t *Ld, real_t *Ll, real_t *z, real_t *res)
{
    int i,n;

    // Fill in first NX elements
    cblas_dsymv(CblasColMajor, CblasLower, NX, 1.0, &Ld[0], 
            NX, &z[0], 1, 0.0, &res[0], 1);
    
    myblas_dgemv(CblasColMajor, CblasTrans, NX, NX, 1.0, &Ll[0], 
            NX, &z[NX], 1, 1.0, &res[0], 1); 
    
    // Fill in the rest of the vector
    for(n=0;n<NI-1;n++)
    {           
        myblas_dgemv(CblasColMajor, CblasNoTrans, NX, NX, 1.0, 
                &Ll[n*NX*NX], NX, &z[n*NX], 1, 0.0, &res[(n+1)*NX], 1);
        
        cblas_dsymv(CblasColMajor, CblasLower, NX, 1.0, &Ld[(n+1)*NX*NX],
                NX, &z[(n+1)*NX], 1, 1.0, &res[(n+1)*NX], 1);
        
        myblas_dgemv(CblasColMajor, CblasTrans, NX, NX, 1.0, 
                &Ll[(n+1)*NX*NX], NX, &z[(n+2)*NX], 1, 1.0, 
                &res[(n+1)*NX], 1);
    }
    
    myblas_dgemv(CblasColMajor, CblasNoTrans, NX, NX, 1.0,
            &Ll[n*NX*NX], NX, &z[n*NX], 1, 0.0, &res[(n+1)*NX], 1);
    
    cblas_dsymv(CblasColMajor, CblasLower, NX, 1.0, &Ld[(n+1)*NX*NX],
            NX, &z[(n+1)*NX], 1, 1.0, &res[(n+1)*NX], 1);
    
    

}
*/

/* 
void fullToBandedFactorization(real_t *Ld, real_t *Ll, real_t *Lb)
{
    int i,n,k,m;
    
    
    for(i=0;i<NX;i++) // fill first NX diagonals
    {
        m = i*ND; // pointer to current element of result
        
        for(n=0;n<=NI;n++)
        {
            for(k=0;k<NX-i;k++) // fill in elements from diag blocks
            {
                Lb[m++] = Ld[n*NX*NX+i+k*NX+k];
            }
            for(k=0;k<i;k++) // fill in elements from off-diag blocks
            {
                Lb[m++] = Ll[(n+1)*NX*NX-i*NX+k*NX+k];
            }
            
        }
        
    }
    
    for(i=0;i<NX;i++) // fill in next NX diagonals
    {
        m = (i+NX)*ND;
        
        for(n=0;n<NI;n++)
        {
            for(k=0;k<NX-i;k++)
            {
                Lb[m++] = Ll[n*NX*NX+i+k*NX+k];
            }
            
            m = m + i;
            
        }
        
    }
    
}
 */
