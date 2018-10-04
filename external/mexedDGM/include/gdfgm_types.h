

#ifndef GDFGM_TYPES_H
#define GDFGM_TYPES_H

#ifndef USE_IN_MATLAB
    #define STOP_GRADIENT_MAP = 1 //default option at the moment
#endif


#ifndef REAL_T
#define REAL_T
#ifdef __USE_SINGLE_PRECISION__
  typedef float real_t;
#else
  typedef double real_t;
#endif	/* __USE_SINGLE_PRECISION__ */
#endif

    
#define YES 1
#define NO 0

/* options structure for generalized dual fast gradient method */
typedef struct
{
	int     maxIt;      // maximum number of iterations
	real_t  tol;        // tolerance in termination condition
        int     warmStart;  // warmstart method
	int     calcMu;     // calculate multipliers of inequality constraints
	int     useExt;     // use external packages (currently BLAS, LAPACK)

} gdfgm_options;


/* QP data structure for generalized dual fast gradient method */
typedef struct
{
    // objective
    real_t *Hinv;
    real_t *f;
    
    // equality constraints
    real_t *Am;
    real_t *Bm;
    real_t *beq;
    
    // bounds
    real_t *ub;
    real_t *lb;
    
    // intermediate results
    real_t *Ld;
    real_t *Ll;
    real_t *Lb;
    
    // a-priori information
    real_t *optSol;
    real_t *lam;
    
    // results
    real_t *sol;
    real_t *fval;
    real_t *mu;
    
    //info
    real_t cpuTime;
    int nIter;
    
} gdfgm_data;


#endif
