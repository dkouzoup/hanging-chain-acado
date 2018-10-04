

#ifndef GDFGM_LOGDATA_H
#define GDFGM_LOGDATA_H

#ifndef USE_IN_MATLAB
#ifndef NUM_STEPS
#include "../../acado_common.h"
#endif
#endif

#if LOG_DATA > 0 // LOG ONLY INIT AND ITER TIMES

timer gdfgm_total_tmr; // timing of solver outside the function
timer gdfgm_initialization_tmr; // time of the whole initialization
timer gdfgm_iterations_tmr; // time for all iterations

// scalar values for updates within solver
double gdfgm_total_time;
double gdfgm_initialization_time;
double gdfgm_iterations_time;

// matrices for logging results of each RTI iteration
double Gdfgm_total_time[NUM_STEPS];
double Gdfgm_initialization_time[NUM_STEPS];
double Gdfgm_iterations_time[NUM_STEPS];

int Gdfgm_iterations_number[NUM_STEPS];

#endif

#if LOG_DATA > 1 // LOG IN ADDITION ALL INITIALIZATION PARTS

timer gdfgm_cholesky_tmr; // block Cholesky factorization of matrix L (data stored back in L)

timer gdfgm_buildL_tmr; // build matrix L (i.e. Ld and Ll)

timer gdfgm_convertToBanded_tmr; // convert Cholesky factor of L to banded form

double gdfgm_cholesky_time;
double gdfgm_buildL_time;
double gdfgm_convertToBanded_time;

double Gdfgm_cholesky_time[NUM_STEPS];
double Gdfgm_buildL_time[NUM_STEPS];
double Gdfgm_convertToBanded_time[NUM_STEPS];

#endif

#if LOG_DATA > 2 // LOG IN ADDITION ALL ITERATION PARTS

timer gdfgm_AeqTransTimesVector_tmr; // Aeq'*v

timer gdfgm_AeqTimesVector_tmr; // Aeq*v

timer gdfgm_solves_tmr; // backward and forward solves

timer gdfgm_termCondition_tmr; // evaluation of termination condition

double gdfgm_AeqTransTimesVector_time;
double gdfgm_AeqTimesVector_time;
double gdfgm_solves_time;
double gdfgm_termCondition_time;

double Gdfgm_AeqTransTimesVector_time[NUM_STEPS];
double Gdfgm_AeqTimesVector_time[NUM_STEPS];
double Gdfgm_solves_time[NUM_STEPS];
double Gdfgm_termCondition_time[NUM_STEPS];

#endif

#ifdef LOG_LOOP
#ifndef USE_IN_MATLAB
timer integration_tmr; //system simulation for closed-loop (don't care)
timer feedback_tmr; // calculation of feedback control
timer preparation_tmr; // preparation of next QP with all data BUT the initial condition
timer oneloop_tmr; // time for one whole loop

double Integration_time[NUM_STEPS];
double Feedback_time[NUM_STEPS];
double Preparation_time[NUM_STEPS];
double Oneloop_time[NUM_STEPS];
#endif
#endif

#endif