
#include <stdio.h>
#include <math.h>

#include "include/gdfgm_dimensions.h"
#include "include/gdfgm_utils.h"
#include "include/gdfgm_algebra.h"
#include "include/gdfgm_types.h"
#include "include/timing.h"
#include "include/gdfgm_logData.h"
#include "include/gdfgm_core.h"


int gdfgm_solve (real_t *sol, real_t *fval, real_t *Hinv, real_t *f, 
        real_t *Am, real_t *Bm, real_t *beq, real_t *ub, real_t *lb, 
        real_t *Ld, real_t *Ll, real_t *Lb, real_t *lam, real_t *mu, 
        real_t *optSol, gdfgm_options opt)
{
    
    // Declarations   
    int i,j;
    int nIter;
    
    real_t lambda[ND];                         // dual vectors
    real_t z[NP], zunc[NP];                    // primal vector (unconstrained and clipped)
   
    #if LOG_DATA > 2
    // initialization of accumulative timers
    gdfgm_AeqTransTimesVector_time = 0;
    gdfgm_AeqTimesVector_time = 0;
    gdfgm_solves_time = 0;
    gdfgm_termCondition_time = 0;
    #endif

    #if LOG_DATA > 0
    tic(&gdfgm_initialization_tmr);      
    #endif

    // R = Aeq*Hinv*Aeq'
    #if LOG_DATA > 1
    tic(&gdfgm_buildL_tmr);
    #endif
    buildMatrixL(Hinv, Am, Bm, Ld, Ll);
    #if LOG_DATA > 1
    gdfgm_buildL_time = toc(&gdfgm_buildL_tmr);
    #endif

    // Banded Cholesky factorization 
    #if LOG_DATA > 1  
    tic(&gdfgm_cholesky_tmr);
    #endif
    blockCholeskyFactorization(Ld, Ll);
    #if LOG_DATA > 1
    gdfgm_cholesky_time = toc(&gdfgm_cholesky_tmr);
    #endif
    
    // Build banded 
    #if LOG_DATA > 1
    tic(&gdfgm_convertToBanded_tmr);
    #endif
    #ifdef USE_BAND_SOLVER
    bandedFactorization(Ld, Ll, Lb);
    #endif
    #if LOG_DATA > 1
    gdfgm_convertToBanded_time = toc(&gdfgm_convertToBanded_tmr);
    #endif

    #if LOG_DATA > 0
    gdfgm_initialization_time = toc(&gdfgm_initialization_tmr);
    #endif

    // Start iterations
    #if LOG_DATA > 0
    tic(&gdfgm_iterations_tmr);
    #endif
    
    runFGMiterations(Hinv,f, Am, Bm, beq, ub, lb, 
		     Ld, Ll, z, zunc, lambda, lam, &nIter, opt);
    
    #if LOG_DATA > 0
    gdfgm_iterations_time = toc(&gdfgm_iterations_tmr);
    #endif

    // Copy primal vector to solution
    myblas_dcopy(NP,z,1,sol,1);

    // Copy dual vector for warm-starting
    myblas_dcopy(ND,lambda,1,lam,1);
            
    // Calculate multipliers for inequality constraints
    if (opt.calcMu == YES)
    {
        myblas_dcopy(NP,zunc,1,z,1);
        
        myblas_daxpy(NP, -1.0, ub, 1, z, 1);
        for(i=0;i<NP;i++)
        {
            z[i]  = max(z[i],0);
            mu[i] = z[i]/Hinv[i];
        }
        
        myblas_dcopy(NP,zunc,1,z,1);
        
        myblas_daxpy(NP, -1.0, lb, 1, z, 1);
        for(i=0;i<NP;i++)
        {
            z[i]     = max(-z[i],0);
            mu[i+NP] = z[i]/Hinv[i];
        }
    }


    // Set objective value
    fval[0] = -1; // TODO
        
    return nIter;
    
}
