

#include "include/gdfgm_dimensions.h"
#include "include/gdfgm_utils.h"
#include "include/gdfgm_algebra.h"
#include "include/gdfgm_types.h"
#include "include/gdfgm_core.h"
#include "include/gdfgm_logData.h"

#include <math.h>

#if LOG_DATA > 2
#include "include/timing.h"
#endif

void runFGMiterations(real_t Hinv[NP], real_t f[NP], 
		      real_t Am[NI*NX*NX], real_t Bm[NI*NX*NU], real_t beq[(NI+1)*NX], 
		      real_t ub[NP], real_t lb[NP], 
		      real_t Ld[(NI+1)*NX*NX], real_t Ll[NI*NX*NX], 
		      real_t z[NP], real_t zunc[NP], real_t lambda[ND], real_t yInit[ND],
		      int *nIter, gdfgm_options opt)
{
    int iter;
    int terminate;
    int i;
    
    //real_t zunc[NP];                         // primal vector (unconstrained)
    real_t lambda_old[ND], y[ND];              // dual vectors
    
#ifdef STOP_GRADIENT_MAP
    real_t errorNorm;                          // norm of gradient map for termination condition
    real_t err[ND];                            // error vector
#endif
    
#ifdef STOP_OPTIMAL_SOLUTION    
    real_t dist[NP];                           // distance from opt sol
    real_t optSolNorm;
#endif

    real_t beta;                               /// coef. for dual FGM
    
#ifdef USE_PRECALCULATED_BETAS
#include "include/betas.h"
#else
    real_t alpha = 1.0;
    real_t alpha_new;                          
#endif
    
    // Initialize First-Order Method
    iter = 0;
    terminate = 0;
    
#ifdef STOP_OPTIMAL_SOLUTION
    optSolNorm = myblas_dnrm2(NP,optSol,1);
#endif
    
    // Warm-start dual multipliers
    myblas_dcopy(ND,yInit,1,y,1);
    
    // Initialize vectors to zero
    initVector(lambda_old, ND);
    initVector(lambda, ND);
    initVector(zunc, NP);
    
    while (iter < opt.maxIt && terminate == 0)
    {
	// z = Aeq'*y
#if LOG_DATA > 2
	tic(&gdfgm_AeqTransTimesVector_tmr);
#endif
	AeqTransposeTimesVector(Am, Bm, y, zunc);
#if LOG_DATA > 2
	gdfgm_AeqTransTimesVector_time += toc(&gdfgm_AeqTransTimesVector_tmr);
#endif
	
	// z += f
	myblas_daxpy(NP, 1.0, f, 1, zunc, 1);
	
	// z = -Hinv*z, clipped to feasible set
	for(i=0;i<NP;i++)
        {
	    zunc[i] = -Hinv[i]*zunc[i];
	    z[i]    = max(min(zunc[i],ub[i]),lb[i]);
        }
	
	// lambda = Aeq*z
#if LOG_DATA > 2
	tic(&gdfgm_AeqTimesVector_tmr);
#endif
	AeqTimesVector(Am, Bm, z, lambda);
#if LOG_DATA > 2
	gdfgm_AeqTimesVector_time += toc(&gdfgm_AeqTimesVector_tmr);
#endif
	
	// lambda -= beq
	myblas_daxpy(ND, -1.0, beq, 1, lambda, 1);
	
#ifdef STOP_GRADIENT_MAP
	myblas_dcopy(ND,lambda,1,err,1); // store gradient map
#endif
	
	// Backward and forward solves
	// lambda *= Ldinv
#if LOG_DATA > 2    
	tic(&gdfgm_solves_tmr);
#endif
	
#if defined(USE_BAND_SOLVER) && defined(USE_EXTERNAL_LIBRARIES)
	cblas_dtbsv(CblasColMajor,CblasLower,CblasNoTrans,CblasNonUnit,
		    ND,2*NX-1,&Lb[0],2*NX,lambda,1);
	cblas_dtbsv(CblasColMajor,CblasLower, CblasTrans,CblasNonUnit,
		    ND,2*NX-1,&Lb[0],2*NX,lambda,1);
#else
#ifdef USE_BAND_SOLVER
#error Band-solver required BLAS but USE_EXTERNAL_LIBRARIES not defined.
#endif
	blockForwardSubstitution(Ld, Ll, lambda);	
	blockBackwardSubstitution(Ld, Ll, lambda);
#endif
	
#if LOG_DATA > 2    
        gdfgm_solves_time += toc(&gdfgm_solves_tmr);
#endif
	
        // lambda += y
        myblas_daxpy(ND, 1.0, y, 1, lambda, 1);
        
        // Check termination condition
        // norm(Ld*(y - lambda),2) < tol
#if LOG_DATA > 2
        tic(&gdfgm_termCondition_tmr);
#endif
	
#ifdef STOP_GRADIENT_MAP      
	errorNorm = myblas_dnrm2(ND,err,1);
        
	if(errorNorm < opt.tol)
	    terminate = 1;
#endif
	
#ifdef STOP_OPTIMAL_SOLUTION    
	
	myblas_dcopy(NP,z,1,dist,1);
	myblas_daxpy(NP,-1.0,optSol,1,dist,1);
        
	errorNorm = myblas_dnrm2(NP,dist,1);
	errorNorm = errorNorm/optSolNorm;
	if(errorNorm < opt.tol)
	    terminate = 1; 
#endif
	
#if LOG_DATA > 2    
        gdfgm_termCondition_time += toc(&gdfgm_termCondition_tmr);
#endif
	
        // Update alpha and beta sequences
#ifndef USE_PRECALCULATED_BETAS        
        alpha_new = (alpha/2.0)*(sqrt(alpha*alpha+4.0)-alpha);
        beta      = (alpha*(1.0-alpha))/(alpha*alpha + alpha_new);
        alpha     = alpha_new;
#else
        beta = betas[iter];
#endif
	
        // y = lambda + beta*(lambda-lambda_old)
        
        // y = 0
        initVector(y, ND); 
        // y = -beta*lambda_old
        myblas_daxpy(ND,-beta,lambda_old,1,y,1);
        // y = (beta+1)*lambda - beta*lambda_old
        myblas_daxpy(ND,beta+1,lambda,1,y,1); 
        
        // lambda_old <-- lambda
        myblas_dcopy(ND,lambda,1,lambda_old,1);

	iter++;
    }
    
    *nIter = iter;
}
