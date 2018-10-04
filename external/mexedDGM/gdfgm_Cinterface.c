

/* Functions for C implementation */

#include <stdio.h>

#include "include/gdfgm_dimensions.h"
#include "include/gdfgm_types.h"
#include "include/gdfgm_Cinterface.h"
#include "include/gdfgm_algebra.h"
#include "include/gdfgm_data.h"

#ifndef USE_IN_MATLAB

// Print matrix with specified format and way of storage 
/*
void printMatrix(	const char* name,
                    real_t* mat,
                    unsigned nRows,
                    unsigned nCols,
                    unsigned major) // 1 for row 2 for column major
{
    unsigned r, c;
    printf("%s: \n", name);
    
    if (major == 1)
    {
        for (r = 0; r < nRows; ++r)
        {
            for(c = 0; c < nCols; ++c)
                printf("\t%f", mat[r * nCols + c]);
            printf("\n");
        }
    }
    else if (major == 2)
    {
        for (r = 0; r < nRows; ++r)
        {
            for(c = 0; c < nCols; ++c)
                printf("\t%f", mat[c* nRows + r]);
            printf("\n");
        }
        
    }
    
}
*/

void gfgm_init(gdfgm_data *dat)
{
    
    // initialize data to zero
    initVector(gdfgm_Am, NI*NX*NX);
    initVector(gdfgm_Bm, NI*NX*NU);
    initVector(gdfgm_beq,(NI+1)*NX);
    initVector(gdfgm_Hinv,NP);
    initVector(gdfgm_f,NP);
    initVector(gdfgm_ub,NP);
    initVector(gdfgm_lb,NP);

    initVector(gdfgm_Ld, (NI+1)*NX*NX);
    initVector(gdfgm_Ll, NI*NX*NX);
    initVector(gdfgm_Lb, 2*NX*ND);
    
    initVector(gdfgm_optSol, NP);   // TODO: for stopping criterion
    initVector(gdfgm_lam, ND);      // TODO: for warm-start
    
    initVector(gdfgm_sol, NP);
    initVector(gdfgm_mu, 2*NP);
    
    gdfgm_fval = 0.0;
    
    // setup structure of pointers
    dat->Am   = gdfgm_Am;
    dat->Bm   = gdfgm_Bm;
    dat->beq  = gdfgm_beq;
    dat->Hinv = gdfgm_Hinv;
    dat->f    = gdfgm_f;
    dat->ub   = gdfgm_ub;
    dat->lb   = gdfgm_lb;
    
    dat->Ld   = gdfgm_Ld;
    dat->Ll   = gdfgm_Ll;
    dat->Lb   = gdfgm_Lb;
    
    dat->optSol = gdfgm_optSol;
    dat->lam    = gdfgm_lam;
    
    dat->sol  = gdfgm_sol;
    dat->mu   = gdfgm_mu;
    dat->fval = &gdfgm_fval;

}


void gfgm_update(gdfgm_data *dat, real_t *H, real_t *G, real_t *C, real_t *c, real_t *lb, real_t *ub, real_t *x0relative)
{
    int i,j,k;
    real_t Atemp[NX*NX];
    real_t Btemp[NX*NU];
    real_t ABtemp[NX*(NX+NU)];
    
    // set intermediate previous results to zero
    initVector(gdfgm_Ld, (NI+1)*NX*NX);
    initVector(gdfgm_Ll, NI*NX*NX);
    initVector(gdfgm_Lb, 2*NX*ND);
    initVector(gdfgm_sol, NP);
    
    for(k=0;k<NI;k++)
    {
        // extract [A_k B_k]
        for(j=0;j<NX*(NX+NU);j++)
        {
            ABtemp[j] = C[k*(NX*(NX+NU))+j]; 
        }
        
        // extract A_k
        for(i=0;i<NX;i++)
        {
            for(j=0;j<NX;j++)
            {
                Atemp[i*NX+j] = ABtemp[i*(NX+NU)+j];
            }
        }
        
        // extract B_k
        for(i=0;i<NX;i++)
        {
            for(j=0;j<NU;j++)
            {
                Btemp[i*NU+j] = ABtemp[i*(NX+NU)+j+NX];
            }
        }        
        
        // transpose A_k and B_k
        transposeMatrix(Atemp, NX, NX, 1); 
        transposeMatrix(Btemp, NX, NU, 1); 
        
        //printMatrix("ABtemp",ABtemp,NX,NX+NU,1);
        //printMatrix("Atemp",Atemp,NX,NX,2);
        //printMatrix("Btemp",Btemp,NX,NU,2);
        
        // store A_k to Am
        for(j=0;j<NX*NX;j++)
        {
            (*dat).Am[k*NX*NX+j] = Atemp[j]; 
        }
        
        // store B_k to Bm
        for(j=0;j<NX*NU;j++)
        {
            (*dat).Bm[k*NX*NU+j] = Btemp[j];
        }

        // extract H_k, f_k, ub_k and lb_k
        for(j=0;j<NX+NU;j++)
        {
            (*dat).Hinv[k*(NX+NU)+j] = 1.0/H[k*(NX+NU)*(NX+NU)+ j*(NX+NU) + j];
            (*dat).f[k*(NX+NU)+j]    = G[k*(NX+NU) + j]; 
            (*dat).ub[k*(NX+NU)+j]   = ub[k*(NX+NU)+ j];
            (*dat).lb[k*(NX+NU)+j]   = lb[k*(NX+NU)+ j];
        }
        
        // build beq
        for(j=0;j<NX;j++)
        {
            (*dat).beq[(k+1)*NX + j] = -c[k*NX + j];
            
        }

           
    }
    
    //extract H_N and f_N
    for(j=0;j<NX;j++)
    {
        (*dat).Hinv[NI*(NX+NU)+j] = 1.0/H[NI*(NX+NU)*(NX+NU)+ j*(NX) + j];
        (*dat).f[NI*(NX+NU)+j]    = G[NI*(NX+NU) + j];
        (*dat).ub[NI*(NX+NU)+j]   = ub[NI*(NX+NU) + j];
        (*dat).lb[NI*(NX+NU)+j]   = lb[NI*(NX+NU) + j];

    }
    
    // initial value embedding
    for(j=0;j<NX;j++)
    {
        (*dat).beq[j] = x0relative[j];
    }
    
}

#endif
