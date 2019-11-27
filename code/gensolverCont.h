/* Produced by CVXGEN, 2017-08-04 19:30:31 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.h. */
/* Description: Header file with relevant definitions. */
#ifndef GENSOLVERCONT_H
#define GENSOLVERCONT_H
#include <array>
#include <ctime>
#include <cstdlib>
/* Uncomment the next line to remove all library dependencies. */
/*#define ZERO_LIBRARY_MODE */
#ifdef MATLAB_MEX_FILE
/* Matlab functions. MATLAB_MEX_FILE will be defined by the mex compiler. */
/* If you are not using the mex compiler, this functionality will not intrude, */
/* as it will be completely disabled at compile-time. */
//#include "mex.h"
#else
#ifndef ZERO_LIBRARY_MODE
#include <iostream>
#include <stdio.h>
#endif
#endif
/* Space must be allocated somewhere (testsolver.c, csolve.c or your own */
/* program) for the global variables vars, params, work and settings. */
/* At the bottom of this file, they are externed. */
#ifndef ZERO_LIBRARY_MODE
#include <cmath>
#include <math.h>
#define pm(A, m, n) printmatrix(#A, A, m, n, 1)
#endif
using namespace std;
class GensolverCont {
private:
typedef struct Params_t {
  int caseNumber;
  double c2[1];
  double c1[1];
  double c0[1];
  double betaSC[1];
  double PgAPPSC[1];
  double gammaSC[1];
  double BSC[1];
  double lambda_2SC[1];
  double rho[1];
  double Pg_N_init[1];
  double Pg_N_avg[1];
  double ug_N[1];
  double Vg_N_avg[1];
  double Thetag_N_avg[1];
  double vg_N[1];
  double PgMin[1];
  double PgMax[1];
} Params;
typedef struct Vars_t {
  double *Pg; /* 1 rows. */
  double *Thetag; /* 1 rows. */
} Vars;
typedef struct Workspace_t {
  double h[2];
  double s_inv[2];
  double s_inv_z[2];
  double *b;
  double q[2];
  double rhs[6];
  double x[6];
  double *s;
  double *z;
  double *y;
  double lhs_aff[6];
  double lhs_cc[6];
  double buffer[6];
  double buffer2[6];
  double KKT[10];
  double L[4];
  double d[6];
  double v[6];
  double d_inv[6];
  double gap;
  double optval;
  double ineq_resid_squared;
  double eq_resid_squared;
  double block_33[1];
  /* Pre-op symbols. */
  double frac_957443481600;
  double frac_121674190848;
  double quad_583541456896[1];
  double quad_155015135232[1];
  double quad_497818112000[1];
  int converged;
} Workspace;
typedef struct Settings_t {
  double resid_tol;
  double eps;
  int max_iters;
  int refine_steps;
  int better_start;
  /* Better start obviates the need for s_init and z_init. */
  double s_init;
  double z_init;
  int verbose;
  /* Show extra details of the iterative refinement steps. */
  int verbose_refinement;
  int debug;
  /* For regularization. Minimum value of abs(D_ii) in the kkt D factor. */
  double kkt_reg;
} Settings;
Vars vars;
Params params;
Workspace work;
Settings settings;
long global_seed;
clock_t tic_timestart;
double *Piterate;
double *Thiterate;
public:
GensolverCont( int, double, double, double, double, double, double, double ); // Gensolver object constructor
~GensolverCont(); // Gensolver object destructor
double getPMax(); // Gets the maximum power generating limit
double getPMin(); // Gets the minimum power generating limit
double getQuadCoeff(); // Gets the Quadratic cost coefficient
double getLinCoeff(); // Gets the Linear cost coefficient
double getConstCoeff(); // Gets the Constant coefficient
/* Function definitions in ldl.c: */
void ldl_solve(double *target, double *var);
void ldl_factor(void);
double check_factorization(void);
void matrix_multiply(double *result, double *source);
double check_residual(double *target, double *multiplicand);
void fill_KKT(void);

/* Function definitions in matrix_support.c: */
void multbymA(double *lhs, double *rhs);
void multbymAT(double *lhs, double *rhs);
void multbymG(double *lhs, double *rhs);
void multbymGT(double *lhs, double *rhs);
void multbyP(double *lhs, double *rhs);
void fillq(void);
void fillh(void);
void fillb(void);
void pre_ops(void);

/* Function definitions in solver.c: */
double eval_gap(void);
void set_defaults(void);
void setup_pointers(void);
void setup_indexing(void);
void set_start(void);
double eval_objv(void);
void fillrhs_aff(void);
void fillrhs_cc(void);
void refine(double *target, double *var);
double calc_ineq_resid_squared(void);
double calc_eq_resid_squared(void);
void better_start(void);
void fillrhs_start(void);
long solve(void);

/* Function definitions in testsolver.c: */
void mainsolve(int, double, double, double, double, double, double, double, double, double, double);
void load_default_data(int, double, double, double, double, double, double, double, double, double, double);
double getPSol(void);
double *getThetaPtr(void);
double getObj(void);
double getBeta();
double getGamma();

/* Function definitions in util.c: */
void tic(void);
float toc(void);
float tocq(void);
void printmatrix(char *name, double *A, int m, int n, int sparse);
double unif(double lower, double upper);
float ran1(long*idum, int reset);
float randn_internal(long *idum, int reset);
double randn(void);
void reset_rand(void);
};

#endif
