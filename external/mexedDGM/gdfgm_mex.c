
#include <stdio.h>

#include "include/timing.h"
#include "include/gdfgm_dimensions.h"
#include "include/gdfgm_algebra.h"
#include "include/gdfgm_types.h"

#include "mex.h"

/* Input Arguments */

#define	H_IN   prhs[0]
#define	F_IN   prhs[1]
#define A_IN   prhs[2]
#define B_IN   prhs[3]
#define BEQ_IN prhs[4]
#define UB_IN  prhs[5]
#define LB_IN  prhs[6]
#define OPT_IN prhs[7] // options structure
#define LAM_IN prhs[8] // warm-start for lambdas
#define SOL_IN prhs[9] // optimal solution

/* Output Arguments */

#define SOL_OUT  plhs[0]
#define FVAL_OUT plhs[1]
#define CPU_OUT  plhs[2]
#define ITER_OUT plhs[3]
#define LAM_OUT  plhs[4]
#define MU_OUT   plhs[5]

void readOptionsStruct(const mxArray* optionsMatlabPtr, gdfgm_options *optionsPtr); 
void getOptionValue( const mxArray* optionsMatlabPtr, const char* optionString, double **optionValue);

int gdfgm_solve (double *sol, double *fval, double *H, double *f,
        double *Am, double *Bm, double *beq, double *ub, double *lb,
        double *Ld, double *Ll, double *Lb, double *lam, double *mu,
        double *optSol, gdfgm_options opt);


/* Main function */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )     
{ 

    // Declarations
    int i, j, niter;

    double *H, *f, *Am, *Bm, *beq, *ub, *lb, *lam,*it ;       // inputs
    double *sol, *tol, *fval, *cputime, *iter, *lambda, *mu;  // outputs
    
    double *optSol = 0; // pointer to known optimal solution
    double Hinv[NP], Ld[(NI+1)*NX*NX], Ll[NI*NX*NX], Lb[2*NX*ND]; // data
    int termCond;

    const mxArray *optionsMatlabPtr;
    
    gdfgm_options  opt;
    gdfgm_options* optPtr = &opt;
            
    timer tmr1, tmr2;
    real_t time1,time2;
    
    
    tic(&tmr1);
    
    // Check number of inputs
    if(nrhs<8 || nrhs >10) {
        mexErrMsgIdAndTxt("MyToolbox:DFGM:nrhs",
                "Wrong number of inputs.");
    }
      
    // Initialization of pointers
    H   = mxGetPr(H_IN);
    f   = mxGetPr(F_IN);
    Am  = mxGetPr(A_IN);
    Bm  = mxGetPr(B_IN);
    beq = mxGetPr(BEQ_IN);
    ub  = mxGetPr(UB_IN);
    lb  = mxGetPr(LB_IN);
    lam = mxGetPr(LAM_IN);
            
    
    //TODO check if structure first
    // setup options struct
    optionsMatlabPtr = OPT_IN;
    readOptionsStruct(optionsMatlabPtr, optPtr);
    

    #ifdef STOP_OPTIMAL_SOLUTION
        if(nrhs == 10)
        {
            optSol = mxGetPr(SOL_IN);
        }
        else
        {
            mexErrMsgTxt( "ERROR: Optimal solution missing." );
        }
    #endif
    
    
    // Check number of inputs again, after setting options
    if (opt.warmStart == YES)
    {
        if (nrhs < 9)
        {
            mexErrMsgTxt( "ERROR: Multipliers for warmStart are missing." );
        }
        // TODO: check also size of vector
    }
    

    // Memory allocation for outputs
    SOL_OUT  = mxCreateDoubleMatrix(NP, 1, mxREAL);
    FVAL_OUT = mxCreateDoubleMatrix(1,  1, mxREAL);
    CPU_OUT  = mxCreateDoubleMatrix(1,  1, mxREAL);
    ITER_OUT = mxCreateDoubleMatrix(1,  1, mxREAL);
    LAM_OUT  = mxCreateDoubleMatrix(ND, 1, mxREAL);
    MU_OUT   = mxCreateDoubleMatrix(2*NP, 1, mxREAL);
        
    
    // Pointers to the newly allocated memory
    sol     = mxGetPr(SOL_OUT);
    fval    = mxGetPr(FVAL_OUT);
    cputime = mxGetPr(CPU_OUT);
    iter    = mxGetPr(ITER_OUT);
    lambda  = mxGetPr(LAM_OUT);
    mu      = mxGetPr(MU_OUT);
    
    // Initialization of intermediate results
    // *can be done during model simulation and thus excluded from timing    
    initVector(Ld, (NI+1)*NX*NX);
    initVector(Ll, NI*NX*NX);
    initVector(Lb, 2*NX*ND);
    //initVector(L,ND*ND);

    
    // Calculation of Hinv 
    // *can be done ahead of time and thus excluded from timing
    for(i=0;i<NP;i++)
    {
        Hinv[i] = 1/H[i*NP+i];
    }
    
    // Copy dual vector to output lambda
    if (opt.warmStart == YES)
    {
        myblas_dcopy(ND,lam,1,lambda,1);
    }
    else
    {
        initVector(lambda,ND);
    }
    
    //printf("Maximum iterations %d\n",opt.maxIt);
    //printf("Termination condition %d\n",opt.termCond);

    
    tic(&tmr2);
    niter = gdfgm_solve(sol, fval, Hinv, f, Am, Bm, beq, ub, lb, Ld, Ll, Lb, 
            lambda, mu, optSol, opt);
    time2 = toc(&tmr2); 

    *cputime = time2;
    *iter    = niter;
    
    time1 = toc(&tmr1);
    

    //printf("\nTime elapsed for mex function: %f\n",time1);
    //printf("\nTime elapsed for inner function: %f\n",time2);

    return;
    
}

// Transfer data from MATLAB structure to C structure
void readOptionsStruct(const mxArray* optionsMatlabPtr, gdfgm_options *optionsPtr)
{
    mxArray *optionFieldPtr;
    double  *optionValuePtr;
    
    // Read maximum number of iterations
    getOptionValue(optionsMatlabPtr, "maximumIterations", &optionValuePtr);
    optionsPtr->maxIt = (int)*optionValuePtr;
    // Read tolerance
    getOptionValue(optionsMatlabPtr, "tolerance", &optionValuePtr);
    optionsPtr->tol = (double)*optionValuePtr;
    // Read mode for termination condition
    //getOptionValue(optionsMatlabPtr, "terminationCondition", &optionValuePtr);
    //optionsPtr->termCond = (stop_t)*optionValuePtr;
    // Read flag for warm-start
    getOptionValue(optionsMatlabPtr, "warmStart", &optionValuePtr);
    optionsPtr->warmStart = (boolean_t)*optionValuePtr;
    // Read flag for computation of mu (multipliers of ineq. constraints)
    getOptionValue(optionsMatlabPtr, "calculateAllMultipliers", &optionValuePtr);
    optionsPtr->calcMu = (boolean_t)*optionValuePtr;
    // Read flag for external libraries
    getOptionValue(optionsMatlabPtr, "useExternalLibraries", &optionValuePtr);
    optionsPtr->useExt = (boolean_t)*optionValuePtr;
    
}

// Read in a MATLAB structure field (function adapted from qpDUNES)
void getOptionValue( const mxArray* optionsMatlabPtr, const char* optionString, double **optionValue)
{
	mxArray* optionFieldPtr = mxGetField( optionsMatlabPtr,0,optionString );

    if( !mxIsEmpty(optionFieldPtr) )
	{
		if ( ( mxGetM( optionFieldPtr ) != 1 ) || ( mxGetN( optionFieldPtr ) != 1 ) ) {
			mexPrintf( "Error reading options string '%s'.", optionString );
			mexErrMsgTxt( "ERROR: Option value has to be a numerical constant." );
		}

		*optionValue = mxGetPr( optionFieldPtr );
	}
}

