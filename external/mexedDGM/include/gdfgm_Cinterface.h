/** Definition of the floating point data type. */

#ifndef CINTERFACE_H
#define CINTERFACE_H

#ifndef REAL_T
#define REAL_T
#ifdef __USE_SINGLE_PRECISION__
  typedef float real_t;
#else
  typedef double real_t;
#endif	/* __USE_SINGLE_PRECISION__ */
#endif


#ifndef USE_IN_MATLAB

void printMatrix(const char* name, real_t* mat, unsigned nRows,
        unsigned nCols, unsigned major);

// the generalized fast gradient method solver
int gdfgm_solve (real_t *sol, real_t *fval, real_t *H, real_t *f,
        real_t *Am, real_t *Bm, real_t *beq, real_t *ub, real_t *lb,
        real_t *Ld, real_t *Ll, real_t *Lb, real_t *lam, real_t *mu,
        real_t *optSol, gdfgm_options opt);

// solver initialization (in ACADO)
void gfgm_init(gdfgm_data *dat);

// solver data update (in ACADO)
void gfgm_update(gdfgm_data *dat, real_t *H, real_t *G, real_t *C, 
        real_t *c, real_t *lb, real_t *ub, real_t *x0);


#endif


#endif
