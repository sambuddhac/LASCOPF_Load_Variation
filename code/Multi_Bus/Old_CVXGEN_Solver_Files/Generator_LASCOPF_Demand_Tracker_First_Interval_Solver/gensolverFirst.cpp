//solver.cpp
/* Produced by CVXGEN, 2016-05-03 18:52:49 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.c. */
/* Description: Main solver file. */
#include "gensolverFirst.h"
#include <array>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <cmath>
using namespace std;
double GensolverFirst::eval_gap(void) {
  int i;
  double gap;
  gap = 0;
  for (i = 0; i < 6; i++)
    gap += work.z[i]*work.s[i];
  return gap;
}
void GensolverFirst::set_defaults(void) {
  settings.resid_tol = 1e-6;
  settings.eps = 1e-4;
  settings.max_iters = 25;
  settings.refine_steps = 1;
  settings.s_init = 1;
  settings.z_init = 1;
  settings.debug = 0;
  settings.verbose = 0;
  settings.verbose_refinement = 0;
  settings.better_start = 1;
  settings.kkt_reg = 1e-7;
}
void GensolverFirst::setup_pointers(void) {
  work.y = work.x + 13;
  work.s = work.x + 13;
  work.z = work.x + 19;
  vars.Pg = work.x + 0;
  vars.PgNext = work.x + 1;
  vars.Thetag = work.x + 2;
}
void GensolverFirst::setup_indexing(void) {
  setup_pointers();
}
void GensolverFirst::set_start(void) {
  int i;
  for (i = 0; i < 13; i++)
    work.x[i] = 0;
  for (i = 0; i < 0; i++)
    work.y[i] = 0;
  for (i = 0; i < 6; i++)
    work.s[i] = (work.h[i] > 0) ? work.h[i] : settings.s_init;
  for (i = 0; i < 6; i++)
    work.z[i] = settings.z_init;
}
double GensolverFirst::eval_objv(void) {
  int i;
  double objv;
  /* Borrow space in work.rhs. */
  multbyP(work.rhs, work.x);
  objv = 0;
  for (i = 0; i < 13; i++)
    objv += work.x[i]*work.rhs[i];
  objv *= 0.5;
  for (i = 0; i < 13; i++)
    objv += work.q[i]*work.x[i];
  objv += params.c0[0]+work.quad_608032620544[0]+work.quad_643747205120[0]+work.frac_121674190848*(work.quad_324888432640[0]+work.quad_331507691520[0]);
  return objv;
}
void GensolverFirst::fillrhs_aff(void) {
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work.rhs;
  r2 = work.rhs + 13;
  r3 = work.rhs + 19;
  r4 = work.rhs + 25;
  /* r1 = -A^Ty - G^Tz - Px - q. */
  multbymAT(r1, work.y);
  multbymGT(work.buffer, work.z);
  for (i = 0; i < 13; i++)
    r1[i] += work.buffer[i];
  multbyP(work.buffer, work.x);
  for (i = 0; i < 13; i++)
    r1[i] -= work.buffer[i] + work.q[i];
  /* r2 = -z. */
  for (i = 0; i < 6; i++)
    r2[i] = -work.z[i];
  /* r3 = -Gx - s + h. */
  multbymG(r3, work.x);
  for (i = 0; i < 6; i++)
    r3[i] += -work.s[i] + work.h[i];
  /* r4 = -Ax + b. */
  multbymA(r4, work.x);
  for (i = 0; i < 0; i++)
    r4[i] += work.b[i];
}
void GensolverFirst::fillrhs_cc(void) {
  int i;
  double *r2;
  double *ds_aff, *dz_aff;
  double mu;
  double alpha;
  double sigma;
  double smu;
  double minval;
  r2 = work.rhs + 13;
  ds_aff = work.lhs_aff + 13;
  dz_aff = work.lhs_aff + 19;
  mu = 0;
  for (i = 0; i < 6; i++)
    mu += work.s[i]*work.z[i];
  /* Don't finish calculating mu quite yet. */
  /* Find min(min(ds./s), min(dz./z)). */
  minval = 0;
  for (i = 0; i < 6; i++)
    if (ds_aff[i] < minval*work.s[i])
      minval = ds_aff[i]/work.s[i];
  for (i = 0; i < 6; i++)
    if (dz_aff[i] < minval*work.z[i])
      minval = dz_aff[i]/work.z[i];
  /* Find alpha. */
  if (-1 < minval)
      alpha = 1;
  else
      alpha = -1/minval;
  sigma = 0;
  for (i = 0; i < 6; i++)
    sigma += (work.s[i] + alpha*ds_aff[i])*
      (work.z[i] + alpha*dz_aff[i]);
  sigma /= mu;
  sigma = sigma*sigma*sigma;
  /* Finish calculating mu now. */
  mu *= 0.16666666666666666;
  smu = sigma*mu;
  /* Fill-in the rhs. */
  for (i = 0; i < 13; i++)
    work.rhs[i] = 0;
  for (i = 19; i < 25; i++)
    work.rhs[i] = 0;
  for (i = 0; i < 6; i++)
    r2[i] = work.s_inv[i]*(smu - ds_aff[i]*dz_aff[i]);
}
void GensolverFirst::refine(double *target, double *var) {
  int i, j;
  double *residual = work.buffer;
  double norm2;
  double *new_var = work.buffer2;
  for (j = 0; j < settings.refine_steps; j++) {
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 25; i++) {
      residual[i] = residual[i] - target[i];
      norm2 += residual[i]*residual[i];
    }
#ifndef ZERO_LIBRARY_MODE
    if (settings.verbose_refinement) {
      if (j == 0)
        printf("Initial residual before refinement has norm squared %.6g.\n", norm2);
      else
        printf("After refinement we get squared norm %.6g.\n", norm2);
    }
#endif
    /* Solve to find new_var = KKT \ (target - A*var). */
    ldl_solve(residual, new_var);
    /* Update var += new_var, or var += KKT \ (target - A*var). */
    for (i = 0; i < 25; i++) {
      var[i] -= new_var[i];
    }
  }
#ifndef ZERO_LIBRARY_MODE
  if (settings.verbose_refinement) {
    /* Check the residual once more, but only if we're reporting it, since */
    /* it's expensive. */
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 25; i++) {
      residual[i] = residual[i] - target[i];
      norm2 += residual[i]*residual[i];
    }
    if (j == 0)
      printf("Initial residual before refinement has norm squared %.6g.\n", norm2);
    else
      printf("After refinement we get squared norm %.6g.\n", norm2);
  }
#endif
}
double GensolverFirst::calc_ineq_resid_squared(void) {
  /* Calculates the norm ||-Gx - s + h||. */
  double norm2_squared;
  int i;
  /* Find -Gx. */
  multbymG(work.buffer, work.x);
  /* Add -s + h. */
  for (i = 0; i < 6; i++)
    work.buffer[i] += -work.s[i] + work.h[i];
  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 6; i++)
    norm2_squared += work.buffer[i]*work.buffer[i];
  return norm2_squared;
}
double GensolverFirst::calc_eq_resid_squared(void) {
  /* Calculates the norm ||-Ax + b||. */
  double norm2_squared;
  int i;
  /* Find -Ax. */
  multbymA(work.buffer, work.x);
  /* Add +b. */
  for (i = 0; i < 0; i++)
    work.buffer[i] += work.b[i];
  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 0; i++)
    norm2_squared += work.buffer[i]*work.buffer[i];
  return norm2_squared;
}
void GensolverFirst::better_start(void) {
  /* Calculates a better starting point, using a similar approach to CVXOPT. */
  /* Not yet speed optimized. */
  int i;
  double *x, *s, *z, *y;
  double alpha;
  work.block_33[0] = -1;
  /* Make sure sinvz is 1 to make hijacked KKT system ok. */
  for (i = 0; i < 6; i++)
    work.s_inv_z[i] = 1;
  fill_KKT();
  ldl_factor();
  fillrhs_start();
  /* Borrow work.lhs_aff for the solution. */
  ldl_solve(work.rhs, work.lhs_aff);
  /* Don't do any refinement for now. Precision doesn't matter too much. */
  x = work.lhs_aff;
  s = work.lhs_aff + 13;
  z = work.lhs_aff + 19;
  y = work.lhs_aff + 25;
  /* Just set x and y as is. */
  for (i = 0; i < 13; i++)
    work.x[i] = x[i];
  for (i = 0; i < 0; i++)
    work.y[i] = y[i];
  /* Now complete the initialization. Start with s. */
  /* Must have alpha > max(z). */
  alpha = -1e99;
  for (i = 0; i < 6; i++)
    if (alpha < z[i])
      alpha = z[i];
  if (alpha < 0) {
    for (i = 0; i < 6; i++)
      work.s[i] = -z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 6; i++)
      work.s[i] = -z[i] + alpha;
  }
  /* Now initialize z. */
  /* Now must have alpha > max(-z). */
  alpha = -1e99;
  for (i = 0; i < 6; i++)
    if (alpha < -z[i])
      alpha = -z[i];
  if (alpha < 0) {
    for (i = 0; i < 6; i++)
      work.z[i] = z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 6; i++)
      work.z[i] = z[i] + alpha;
  }
}
void GensolverFirst::fillrhs_start(void) {
  /* Fill rhs with (-q, 0, h, b). */
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work.rhs;
  r2 = work.rhs + 13;
  r3 = work.rhs + 19;
  r4 = work.rhs + 25;
  for (i = 0; i < 13; i++)
    r1[i] = -work.q[i];
  for (i = 0; i < 6; i++)
    r2[i] = 0;
  for (i = 0; i < 6; i++)
    r3[i] = work.h[i];
  for (i = 0; i < 0; i++)
    r4[i] = work.b[i];
}
long GensolverFirst::solve(void) {
  int i;
  int iter;
  double *dx, *ds, *dy, *dz;
  double minval;
  double alpha;
  work.converged = 0;
  setup_pointers();
  pre_ops();
#ifndef ZERO_LIBRARY_MODE
  if (settings.verbose)
    printf("iter     objv        gap       |Ax-b|    |Gx+s-h|    step\n");
#endif
  fillq();
  fillh();
  fillb();
  if (settings.better_start)
    better_start();
  else
    set_start();
  for (iter = 0; iter < settings.max_iters; iter++) {
    for (i = 0; i < 6; i++) {
      work.s_inv[i] = 1.0 / work.s[i];
      work.s_inv_z[i] = work.s_inv[i]*work.z[i];
    }
    work.block_33[0] = 0;
    fill_KKT();
    ldl_factor();
    /* Affine scaling directions. */
    fillrhs_aff();
    ldl_solve(work.rhs, work.lhs_aff);
    refine(work.rhs, work.lhs_aff);
    /* Centering plus corrector directions. */
    fillrhs_cc();
    ldl_solve(work.rhs, work.lhs_cc);
    refine(work.rhs, work.lhs_cc);
    /* Add the two together and store in aff. */
    for (i = 0; i < 25; i++)
      work.lhs_aff[i] += work.lhs_cc[i];
    /* Rename aff to reflect its new meaning. */
    dx = work.lhs_aff;
    ds = work.lhs_aff + 13;
    dz = work.lhs_aff + 19;
    dy = work.lhs_aff + 25;
    /* Find min(min(ds./s), min(dz./z)). */
    minval = 0;
    for (i = 0; i < 6; i++)
      if (ds[i] < minval*work.s[i])
        minval = ds[i]/work.s[i];
    for (i = 0; i < 6; i++)
      if (dz[i] < minval*work.z[i])
        minval = dz[i]/work.z[i];
    /* Find alpha. */
    if (-0.99 < minval)
      alpha = 1;
    else
      alpha = -0.99/minval;
    /* Update the primal and dual variables. */
    for (i = 0; i < 13; i++)
      work.x[i] += alpha*dx[i];
    for (i = 0; i < 6; i++)
      work.s[i] += alpha*ds[i];
    for (i = 0; i < 6; i++)
      work.z[i] += alpha*dz[i];
    for (i = 0; i < 0; i++)
      work.y[i] += alpha*dy[i];
    work.gap = eval_gap();
    work.eq_resid_squared = calc_eq_resid_squared();
    work.ineq_resid_squared = calc_ineq_resid_squared();
#ifndef ZERO_LIBRARY_MODE
    if (settings.verbose) {
      work.optval = eval_objv();
      printf("%3d   %10.3e  %9.2e  %9.2e  %9.2e  % 6.4f\n",
          iter+1, work.optval, work.gap, sqrt(work.eq_resid_squared),
          sqrt(work.ineq_resid_squared), alpha);
    }
#endif
    /* Test termination conditions. Requires optimality, and satisfied */
    /* constraints. */
    if (   (work.gap < settings.eps)
        && (work.eq_resid_squared <= settings.resid_tol*settings.resid_tol)
        && (work.ineq_resid_squared <= settings.resid_tol*settings.resid_tol)
       ) {
      work.converged = 1;
      work.optval = eval_objv();
      return iter+1;
    }
  }
  return iter;
}

//util.cpp
/* Produced by CVXGEN, 2016-05-03 18:52:49 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: util.c. */
/* Description: Common utility file for all cvxgen code. */

void GensolverFirst::tic(void) {
  tic_timestart = clock();
}
float GensolverFirst::toc(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  printf("time: %8.2f.\n", (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC);
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}
float GensolverFirst::tocq(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}
void GensolverFirst::printmatrix(char *name, double *A, int m, int n, int sparse) {
  int i, j;
  printf("%s = [...\n", name);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      if ((sparse == 1) && (A[i+j*m] == 0))
        printf("         0");
      else
        printf("  % 9.4f", A[i+j*m]);
    printf(",\n");
  }
  printf("];\n");
}
double GensolverFirst::unif(double lower, double upper) {
  return lower + ((upper - lower)*rand())/RAND_MAX;
}
/* Next function is from numerical recipes in C. */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
float GensolverFirst::ran1(long*idum, int reset) {
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (reset) {
    iy = 0;
  }
  if (*idum<=0||!iy) {
    if (-(*idum)<1)*idum=1;
    else *idum=-(*idum);
    for (j=NTAB+7; j>=0; j--) {
      k = (*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum<0)*idum+=IM;
      if (j<NTAB)iv[j]=*idum;
    }
    iy = iv[0];
  }
  k = (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if (*idum<0)*idum += IM;
  j = iy/NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy)> RNMX) return RNMX;
  else return temp;
}
/* Next function is from numerical recipes in C. */
float GensolverFirst::randn_internal(long *idum, int reset) {
  static int iset=0;
  static float gset;
  float fac, rsq, v1, v2;
  if (reset) {
    iset = 0;
  }
  if (iset==0) {
    do {
      v1 = 2.0*ran1(idum, reset)-1.0;
      v2 = 2.0*ran1(idum, reset)-1.0;
      rsq = v1*v1+v2*v2;
    } while(rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}
double GensolverFirst::randn(void) {
  return randn_internal(&global_seed, 0);
}
void GensolverFirst::reset_rand(void) {
  srand(15);
  global_seed = 1;
  randn_internal(&global_seed, 1);
}

//matrix_support.cpp
/* Produced by CVXGEN, 2016-05-03 18:52:49 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */
void GensolverFirst::multbymA(double *lhs, double *rhs) {
}
void GensolverFirst::multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
  lhs[3] = 0;
  lhs[4] = 0;
  lhs[5] = 0;
  lhs[6] = 0;
  lhs[7] = 0;
  lhs[8] = 0;
  lhs[9] = 0;
  lhs[10] = 0;
  lhs[11] = 0;
  lhs[12] = 0;
}
void GensolverFirst::multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1);
  lhs[1] = -rhs[0]*(1);
  lhs[2] = -rhs[0]*(1)-rhs[1]*(-1);
  lhs[3] = -rhs[0]*(-1)-rhs[1]*(1);
  lhs[4] = -rhs[0]*(-1);
  lhs[5] = -rhs[0]*(1);
}
void GensolverFirst::multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1)-rhs[1]*(1)-rhs[2]*(1)-rhs[3]*(-1)-rhs[4]*(-1)-rhs[5]*(1);
  lhs[1] = -rhs[2]*(-1)-rhs[3]*(1);
  lhs[2] = 0;
  lhs[3] = 0;
  lhs[4] = 0;
  lhs[5] = 0;
  lhs[6] = 0;
  lhs[7] = 0;
  lhs[8] = 0;
  lhs[9] = 0;
  lhs[10] = 0;
  lhs[11] = 0;
  lhs[12] = 0;
}
void GensolverFirst::multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2*(params.c2[0]+work.frac_479733190656+work.frac_121674190848*work.quad_877793640448[0]));
  lhs[1] = rhs[1]*(2*work.frac_479733190656);
  lhs[2] = rhs[2]*(2*work.frac_121674190848);
  lhs[3] = rhs[3]*(2*work.frac_121674190848);
  lhs[4] = rhs[4]*(2*work.frac_121674190848);
  lhs[5] = rhs[5]*(2*work.frac_121674190848);
  lhs[6] = rhs[6]*(2*work.frac_121674190848);
  lhs[7] = rhs[7]*(2*work.frac_121674190848);
  lhs[8] = rhs[8]*(2*work.frac_121674190848);
  lhs[9] = rhs[9]*(2*work.frac_121674190848);
  lhs[10] = rhs[10]*(2*work.frac_121674190848);
  lhs[11] = rhs[11]*(2*work.frac_121674190848);
  lhs[12] = rhs[12]*(2*work.frac_121674190848);
}
void GensolverFirst::fillq(void) {
  work.q[0] = params.c1[0]-2*params.PgNu[0]*work.frac_479733190656+params.gamma[0]*params.B[0]+params.lambda_1[0]+2*work.frac_121674190848*((-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*params.ones[0]+(-params.Pg_N_init[1]+params.Pg_N_avg[1]+params.ug_N[1])*params.ones[1]+(-params.Pg_N_init[2]+params.Pg_N_avg[2]+params.ug_N[2])*params.ones[2]+(-params.Pg_N_init[3]+params.Pg_N_avg[3]+params.ug_N[3])*params.ones[3]+(-params.Pg_N_init[4]+params.Pg_N_avg[4]+params.ug_N[4])*params.ones[4]+(-params.Pg_N_init[5]+params.Pg_N_avg[5]+params.ug_N[5])*params.ones[5]+(-params.Pg_N_init[6]+params.Pg_N_avg[6]+params.ug_N[6])*params.ones[6]+(-params.Pg_N_init[7]+params.Pg_N_avg[7]+params.ug_N[7])*params.ones[7]+(-params.Pg_N_init[8]+params.Pg_N_avg[8]+params.ug_N[8])*params.ones[8]+(-params.Pg_N_init[9]+params.Pg_N_avg[9]+params.ug_N[9])*params.ones[9]+(-params.Pg_N_init[10]+params.Pg_N_avg[10]+params.ug_N[10])*params.ones[10]);
  work.q[1] = -2*params.PgNextNu[0]*work.frac_479733190656+params.gamma[0]*params.D[0]+params.lambda_2[0];
  work.q[2] = 2*work.frac_121674190848*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0]);
  work.q[3] = 2*work.frac_121674190848*(-params.Vg_N_avg[1]-params.Thetag_N_avg[1]+params.vg_N[1]);
  work.q[4] = 2*work.frac_121674190848*(-params.Vg_N_avg[2]-params.Thetag_N_avg[2]+params.vg_N[2]);
  work.q[5] = 2*work.frac_121674190848*(-params.Vg_N_avg[3]-params.Thetag_N_avg[3]+params.vg_N[3]);
  work.q[6] = 2*work.frac_121674190848*(-params.Vg_N_avg[4]-params.Thetag_N_avg[4]+params.vg_N[4]);
  work.q[7] = 2*work.frac_121674190848*(-params.Vg_N_avg[5]-params.Thetag_N_avg[5]+params.vg_N[5]);
  work.q[8] = 2*work.frac_121674190848*(-params.Vg_N_avg[6]-params.Thetag_N_avg[6]+params.vg_N[6]);
  work.q[9] = 2*work.frac_121674190848*(-params.Vg_N_avg[7]-params.Thetag_N_avg[7]+params.vg_N[7]);
  work.q[10] = 2*work.frac_121674190848*(-params.Vg_N_avg[8]-params.Thetag_N_avg[8]+params.vg_N[8]);
  work.q[11] = 2*work.frac_121674190848*(-params.Vg_N_avg[9]-params.Thetag_N_avg[9]+params.vg_N[9]);
  work.q[12] = 2*work.frac_121674190848*(-params.Vg_N_avg[10]-params.Thetag_N_avg[10]+params.vg_N[10]);
}
void GensolverFirst::fillh(void) {
  work.h[0] = -params.PgMin[0];
  work.h[1] = params.PgMax[0];
  work.h[2] = -params.RgMin[0];
  work.h[3] = params.RgMax[0];
  work.h[4] = -(params.RgMin[0]+params.PgPrev[0]);
  work.h[5] = -(-params.PgPrev[0]-params.RgMax[0]);
}
void GensolverFirst::fillb(void) {
}
void GensolverFirst::pre_ops(void) {
  work.frac_479733190656 = params.beta[0];
  work.frac_479733190656 /= 2;
  work.frac_121674190848 = params.rho[0];
  work.frac_121674190848 /= 2;
  work.quad_877793640448[0] = params.ones[0]*params.ones[0]+params.ones[1]*params.ones[1]+params.ones[2]*params.ones[2]+params.ones[3]*params.ones[3]+params.ones[4]*params.ones[4]+params.ones[5]*params.ones[5]+params.ones[6]*params.ones[6]+params.ones[7]*params.ones[7]+params.ones[8]*params.ones[8]+params.ones[9]*params.ones[9]+params.ones[10]*params.ones[10];
  work.quad_608032620544[0] = params.PgNu[0]*work.frac_479733190656*params.PgNu[0];
  work.quad_643747205120[0] = params.PgNextNu[0]*work.frac_479733190656*params.PgNextNu[0];
  work.quad_324888432640[0] = ((-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*(-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])+(-params.Pg_N_init[1]+params.Pg_N_avg[1]+params.ug_N[1])*(-params.Pg_N_init[1]+params.Pg_N_avg[1]+params.ug_N[1])+(-params.Pg_N_init[2]+params.Pg_N_avg[2]+params.ug_N[2])*(-params.Pg_N_init[2]+params.Pg_N_avg[2]+params.ug_N[2])+(-params.Pg_N_init[3]+params.Pg_N_avg[3]+params.ug_N[3])*(-params.Pg_N_init[3]+params.Pg_N_avg[3]+params.ug_N[3])+(-params.Pg_N_init[4]+params.Pg_N_avg[4]+params.ug_N[4])*(-params.Pg_N_init[4]+params.Pg_N_avg[4]+params.ug_N[4])+(-params.Pg_N_init[5]+params.Pg_N_avg[5]+params.ug_N[5])*(-params.Pg_N_init[5]+params.Pg_N_avg[5]+params.ug_N[5])+(-params.Pg_N_init[6]+params.Pg_N_avg[6]+params.ug_N[6])*(-params.Pg_N_init[6]+params.Pg_N_avg[6]+params.ug_N[6])+(-params.Pg_N_init[7]+params.Pg_N_avg[7]+params.ug_N[7])*(-params.Pg_N_init[7]+params.Pg_N_avg[7]+params.ug_N[7])+(-params.Pg_N_init[8]+params.Pg_N_avg[8]+params.ug_N[8])*(-params.Pg_N_init[8]+params.Pg_N_avg[8]+params.ug_N[8])+(-params.Pg_N_init[9]+params.Pg_N_avg[9]+params.ug_N[9])*(-params.Pg_N_init[9]+params.Pg_N_avg[9]+params.ug_N[9])+(-params.Pg_N_init[10]+params.Pg_N_avg[10]+params.ug_N[10])*(-params.Pg_N_init[10]+params.Pg_N_avg[10]+params.ug_N[10]));
  work.quad_331507691520[0] = ((-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])+(-params.Vg_N_avg[1]-params.Thetag_N_avg[1]+params.vg_N[1])*(-params.Vg_N_avg[1]-params.Thetag_N_avg[1]+params.vg_N[1])+(-params.Vg_N_avg[2]-params.Thetag_N_avg[2]+params.vg_N[2])*(-params.Vg_N_avg[2]-params.Thetag_N_avg[2]+params.vg_N[2])+(-params.Vg_N_avg[3]-params.Thetag_N_avg[3]+params.vg_N[3])*(-params.Vg_N_avg[3]-params.Thetag_N_avg[3]+params.vg_N[3])+(-params.Vg_N_avg[4]-params.Thetag_N_avg[4]+params.vg_N[4])*(-params.Vg_N_avg[4]-params.Thetag_N_avg[4]+params.vg_N[4])+(-params.Vg_N_avg[5]-params.Thetag_N_avg[5]+params.vg_N[5])*(-params.Vg_N_avg[5]-params.Thetag_N_avg[5]+params.vg_N[5])+(-params.Vg_N_avg[6]-params.Thetag_N_avg[6]+params.vg_N[6])*(-params.Vg_N_avg[6]-params.Thetag_N_avg[6]+params.vg_N[6])+(-params.Vg_N_avg[7]-params.Thetag_N_avg[7]+params.vg_N[7])*(-params.Vg_N_avg[7]-params.Thetag_N_avg[7]+params.vg_N[7])+(-params.Vg_N_avg[8]-params.Thetag_N_avg[8]+params.vg_N[8])*(-params.Vg_N_avg[8]-params.Thetag_N_avg[8]+params.vg_N[8])+(-params.Vg_N_avg[9]-params.Thetag_N_avg[9]+params.vg_N[9])*(-params.Vg_N_avg[9]-params.Thetag_N_avg[9]+params.vg_N[9])+(-params.Vg_N_avg[10]-params.Thetag_N_avg[10]+params.vg_N[10])*(-params.Vg_N_avg[10]-params.Thetag_N_avg[10]+params.vg_N[10]));
}

//ldl.cpp
/* Produced by CVXGEN, 2016-05-03 18:52:49 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: ldl.c. */
/* Description: Basic test harness for solver.c. */

/* Be sure to place ldl_solve first, so storage schemes are defined by it. */
void GensolverFirst::ldl_solve(double *target, double *var) {
  int i;
  /* Find var = (L*diag(work.d)*L') \ target, then unpermute. */
  /* Answer goes into var. */
  /* Forward substitution. */
  /* Include permutation as we retrieve from target. Use v so we can unpermute */
  /* later. */
  work.v[0] = target[2];
  work.v[1] = target[3];
  work.v[2] = target[4];
  work.v[3] = target[5];
  work.v[4] = target[6];
  work.v[5] = target[7];
  work.v[6] = target[8];
  work.v[7] = target[9];
  work.v[8] = target[10];
  work.v[9] = target[11];
  work.v[10] = target[12];
  work.v[11] = target[13];
  work.v[12] = target[14];
  work.v[13] = target[15];
  work.v[14] = target[16];
  work.v[15] = target[17];
  work.v[16] = target[18];
  work.v[17] = target[19]-work.L[0]*work.v[11];
  work.v[18] = target[20]-work.L[1]*work.v[12];
  work.v[19] = target[23]-work.L[2]*work.v[15];
  work.v[20] = target[24]-work.L[3]*work.v[16];
  work.v[21] = target[0]-work.L[4]*work.v[17]-work.L[5]*work.v[18]-work.L[6]*work.v[19]-work.L[7]*work.v[20];
  work.v[22] = target[1];
  work.v[23] = target[21]-work.L[8]*work.v[13]-work.L[9]*work.v[21]-work.L[10]*work.v[22];
  work.v[24] = target[22]-work.L[11]*work.v[14]-work.L[12]*work.v[21]-work.L[13]*work.v[22]-work.L[14]*work.v[23];
  /* Diagonal scaling. Assume correctness of work.d_inv. */
  for (i = 0; i < 25; i++)
    work.v[i] *= work.d_inv[i];
  /* Back substitution */
  work.v[23] -= work.L[14]*work.v[24];
  work.v[22] -= work.L[10]*work.v[23]+work.L[13]*work.v[24];
  work.v[21] -= work.L[9]*work.v[23]+work.L[12]*work.v[24];
  work.v[20] -= work.L[7]*work.v[21];
  work.v[19] -= work.L[6]*work.v[21];
  work.v[18] -= work.L[5]*work.v[21];
  work.v[17] -= work.L[4]*work.v[21];
  work.v[16] -= work.L[3]*work.v[20];
  work.v[15] -= work.L[2]*work.v[19];
  work.v[14] -= work.L[11]*work.v[24];
  work.v[13] -= work.L[8]*work.v[23];
  work.v[12] -= work.L[1]*work.v[18];
  work.v[11] -= work.L[0]*work.v[17];
  /* Unpermute the result, from v to var. */
  var[0] = work.v[21];
  var[1] = work.v[22];
  var[2] = work.v[0];
  var[3] = work.v[1];
  var[4] = work.v[2];
  var[5] = work.v[3];
  var[6] = work.v[4];
  var[7] = work.v[5];
  var[8] = work.v[6];
  var[9] = work.v[7];
  var[10] = work.v[8];
  var[11] = work.v[9];
  var[12] = work.v[10];
  var[13] = work.v[11];
  var[14] = work.v[12];
  var[15] = work.v[13];
  var[16] = work.v[14];
  var[17] = work.v[15];
  var[18] = work.v[16];
  var[19] = work.v[17];
  var[20] = work.v[18];
  var[21] = work.v[23];
  var[22] = work.v[24];
  var[23] = work.v[19];
  var[24] = work.v[20];
#ifndef ZERO_LIBRARY_MODE
  if (settings.debug) {
    printf("Squared norm for solution is %.8g.\n", check_residual(target, var));
  }
#endif
}
void GensolverFirst::ldl_factor(void) {
  work.d[0] = work.KKT[0];
  if (work.d[0] < 0)
    work.d[0] = settings.kkt_reg;
  else
    work.d[0] += settings.kkt_reg;
  work.d_inv[0] = 1/work.d[0];
  work.v[1] = work.KKT[1];
  work.d[1] = work.v[1];
  if (work.d[1] < 0)
    work.d[1] = settings.kkt_reg;
  else
    work.d[1] += settings.kkt_reg;
  work.d_inv[1] = 1/work.d[1];
  work.v[2] = work.KKT[2];
  work.d[2] = work.v[2];
  if (work.d[2] < 0)
    work.d[2] = settings.kkt_reg;
  else
    work.d[2] += settings.kkt_reg;
  work.d_inv[2] = 1/work.d[2];
  work.v[3] = work.KKT[3];
  work.d[3] = work.v[3];
  if (work.d[3] < 0)
    work.d[3] = settings.kkt_reg;
  else
    work.d[3] += settings.kkt_reg;
  work.d_inv[3] = 1/work.d[3];
  work.v[4] = work.KKT[4];
  work.d[4] = work.v[4];
  if (work.d[4] < 0)
    work.d[4] = settings.kkt_reg;
  else
    work.d[4] += settings.kkt_reg;
  work.d_inv[4] = 1/work.d[4];
  work.v[5] = work.KKT[5];
  work.d[5] = work.v[5];
  if (work.d[5] < 0)
    work.d[5] = settings.kkt_reg;
  else
    work.d[5] += settings.kkt_reg;
  work.d_inv[5] = 1/work.d[5];
  work.v[6] = work.KKT[6];
  work.d[6] = work.v[6];
  if (work.d[6] < 0)
    work.d[6] = settings.kkt_reg;
  else
    work.d[6] += settings.kkt_reg;
  work.d_inv[6] = 1/work.d[6];
  work.v[7] = work.KKT[7];
  work.d[7] = work.v[7];
  if (work.d[7] < 0)
    work.d[7] = settings.kkt_reg;
  else
    work.d[7] += settings.kkt_reg;
  work.d_inv[7] = 1/work.d[7];
  work.v[8] = work.KKT[8];
  work.d[8] = work.v[8];
  if (work.d[8] < 0)
    work.d[8] = settings.kkt_reg;
  else
    work.d[8] += settings.kkt_reg;
  work.d_inv[8] = 1/work.d[8];
  work.v[9] = work.KKT[9];
  work.d[9] = work.v[9];
  if (work.d[9] < 0)
    work.d[9] = settings.kkt_reg;
  else
    work.d[9] += settings.kkt_reg;
  work.d_inv[9] = 1/work.d[9];
  work.v[10] = work.KKT[10];
  work.d[10] = work.v[10];
  if (work.d[10] < 0)
    work.d[10] = settings.kkt_reg;
  else
    work.d[10] += settings.kkt_reg;
  work.d_inv[10] = 1/work.d[10];
  work.v[11] = work.KKT[11];
  work.d[11] = work.v[11];
  if (work.d[11] < 0)
    work.d[11] = settings.kkt_reg;
  else
    work.d[11] += settings.kkt_reg;
  work.d_inv[11] = 1/work.d[11];
  work.L[0] = (work.KKT[12])*work.d_inv[11];
  work.v[12] = work.KKT[13];
  work.d[12] = work.v[12];
  if (work.d[12] < 0)
    work.d[12] = settings.kkt_reg;
  else
    work.d[12] += settings.kkt_reg;
  work.d_inv[12] = 1/work.d[12];
  work.L[1] = (work.KKT[14])*work.d_inv[12];
  work.v[13] = work.KKT[15];
  work.d[13] = work.v[13];
  if (work.d[13] < 0)
    work.d[13] = settings.kkt_reg;
  else
    work.d[13] += settings.kkt_reg;
  work.d_inv[13] = 1/work.d[13];
  work.L[8] = (work.KKT[16])*work.d_inv[13];
  work.v[14] = work.KKT[17];
  work.d[14] = work.v[14];
  if (work.d[14] < 0)
    work.d[14] = settings.kkt_reg;
  else
    work.d[14] += settings.kkt_reg;
  work.d_inv[14] = 1/work.d[14];
  work.L[11] = (work.KKT[18])*work.d_inv[14];
  work.v[15] = work.KKT[19];
  work.d[15] = work.v[15];
  if (work.d[15] < 0)
    work.d[15] = settings.kkt_reg;
  else
    work.d[15] += settings.kkt_reg;
  work.d_inv[15] = 1/work.d[15];
  work.L[2] = (work.KKT[20])*work.d_inv[15];
  work.v[16] = work.KKT[21];
  work.d[16] = work.v[16];
  if (work.d[16] < 0)
    work.d[16] = settings.kkt_reg;
  else
    work.d[16] += settings.kkt_reg;
  work.d_inv[16] = 1/work.d[16];
  work.L[3] = (work.KKT[22])*work.d_inv[16];
  work.v[11] = work.L[0]*work.d[11];
  work.v[17] = work.KKT[23]-work.L[0]*work.v[11];
  work.d[17] = work.v[17];
  if (work.d[17] > 0)
    work.d[17] = -settings.kkt_reg;
  else
    work.d[17] -= settings.kkt_reg;
  work.d_inv[17] = 1/work.d[17];
  work.L[4] = (work.KKT[24])*work.d_inv[17];
  work.v[12] = work.L[1]*work.d[12];
  work.v[18] = work.KKT[25]-work.L[1]*work.v[12];
  work.d[18] = work.v[18];
  if (work.d[18] > 0)
    work.d[18] = -settings.kkt_reg;
  else
    work.d[18] -= settings.kkt_reg;
  work.d_inv[18] = 1/work.d[18];
  work.L[5] = (work.KKT[26])*work.d_inv[18];
  work.v[15] = work.L[2]*work.d[15];
  work.v[19] = work.KKT[27]-work.L[2]*work.v[15];
  work.d[19] = work.v[19];
  if (work.d[19] > 0)
    work.d[19] = -settings.kkt_reg;
  else
    work.d[19] -= settings.kkt_reg;
  work.d_inv[19] = 1/work.d[19];
  work.L[6] = (work.KKT[28])*work.d_inv[19];
  work.v[16] = work.L[3]*work.d[16];
  work.v[20] = work.KKT[29]-work.L[3]*work.v[16];
  work.d[20] = work.v[20];
  if (work.d[20] > 0)
    work.d[20] = -settings.kkt_reg;
  else
    work.d[20] -= settings.kkt_reg;
  work.d_inv[20] = 1/work.d[20];
  work.L[7] = (work.KKT[30])*work.d_inv[20];
  work.v[17] = work.L[4]*work.d[17];
  work.v[18] = work.L[5]*work.d[18];
  work.v[19] = work.L[6]*work.d[19];
  work.v[20] = work.L[7]*work.d[20];
  work.v[21] = work.KKT[31]-work.L[4]*work.v[17]-work.L[5]*work.v[18]-work.L[6]*work.v[19]-work.L[7]*work.v[20];
  work.d[21] = work.v[21];
  if (work.d[21] < 0)
    work.d[21] = settings.kkt_reg;
  else
    work.d[21] += settings.kkt_reg;
  work.d_inv[21] = 1/work.d[21];
  work.L[9] = (work.KKT[32])*work.d_inv[21];
  work.L[12] = (work.KKT[33])*work.d_inv[21];
  work.v[22] = work.KKT[34];
  work.d[22] = work.v[22];
  if (work.d[22] < 0)
    work.d[22] = settings.kkt_reg;
  else
    work.d[22] += settings.kkt_reg;
  work.d_inv[22] = 1/work.d[22];
  work.L[10] = (work.KKT[35])*work.d_inv[22];
  work.L[13] = (work.KKT[36])*work.d_inv[22];
  work.v[13] = work.L[8]*work.d[13];
  work.v[21] = work.L[9]*work.d[21];
  work.v[22] = work.L[10]*work.d[22];
  work.v[23] = work.KKT[37]-work.L[8]*work.v[13]-work.L[9]*work.v[21]-work.L[10]*work.v[22];
  work.d[23] = work.v[23];
  if (work.d[23] > 0)
    work.d[23] = -settings.kkt_reg;
  else
    work.d[23] -= settings.kkt_reg;
  work.d_inv[23] = 1/work.d[23];
  work.L[14] = (-work.L[12]*work.v[21]-work.L[13]*work.v[22])*work.d_inv[23];
  work.v[14] = work.L[11]*work.d[14];
  work.v[21] = work.L[12]*work.d[21];
  work.v[22] = work.L[13]*work.d[22];
  work.v[23] = work.L[14]*work.d[23];
  work.v[24] = work.KKT[38]-work.L[11]*work.v[14]-work.L[12]*work.v[21]-work.L[13]*work.v[22]-work.L[14]*work.v[23];
  work.d[24] = work.v[24];
  if (work.d[24] > 0)
    work.d[24] = -settings.kkt_reg;
  else
    work.d[24] -= settings.kkt_reg;
  work.d_inv[24] = 1/work.d[24];
#ifndef ZERO_LIBRARY_MODE
  if (settings.debug) {
    printf("Squared Frobenius for factorization is %.8g.\n", check_factorization());
  }
#endif
}
double GensolverFirst::check_factorization(void) {
  /* Returns the squared Frobenius norm of A - L*D*L'. */
  double temp, residual;
  /* Only check the lower triangle. */
  residual = 0;
  temp = work.KKT[31]-1*work.d[21]*1-work.L[4]*work.d[17]*work.L[4]-work.L[5]*work.d[18]*work.L[5]-work.L[6]*work.d[19]*work.L[6]-work.L[7]*work.d[20]*work.L[7];
  residual += temp*temp;
  temp = work.KKT[34]-1*work.d[22]*1;
  residual += temp*temp;
  temp = work.KKT[0]-1*work.d[0]*1;
  residual += temp*temp;
  temp = work.KKT[1]-1*work.d[1]*1;
  residual += temp*temp;
  temp = work.KKT[2]-1*work.d[2]*1;
  residual += temp*temp;
  temp = work.KKT[3]-1*work.d[3]*1;
  residual += temp*temp;
  temp = work.KKT[4]-1*work.d[4]*1;
  residual += temp*temp;
  temp = work.KKT[5]-1*work.d[5]*1;
  residual += temp*temp;
  temp = work.KKT[6]-1*work.d[6]*1;
  residual += temp*temp;
  temp = work.KKT[7]-1*work.d[7]*1;
  residual += temp*temp;
  temp = work.KKT[8]-1*work.d[8]*1;
  residual += temp*temp;
  temp = work.KKT[9]-1*work.d[9]*1;
  residual += temp*temp;
  temp = work.KKT[10]-1*work.d[10]*1;
  residual += temp*temp;
  temp = work.KKT[11]-1*work.d[11]*1;
  residual += temp*temp;
  temp = work.KKT[13]-1*work.d[12]*1;
  residual += temp*temp;
  temp = work.KKT[15]-1*work.d[13]*1;
  residual += temp*temp;
  temp = work.KKT[17]-1*work.d[14]*1;
  residual += temp*temp;
  temp = work.KKT[19]-1*work.d[15]*1;
  residual += temp*temp;
  temp = work.KKT[21]-1*work.d[16]*1;
  residual += temp*temp;
  temp = work.KKT[12]-work.L[0]*work.d[11]*1;
  residual += temp*temp;
  temp = work.KKT[14]-work.L[1]*work.d[12]*1;
  residual += temp*temp;
  temp = work.KKT[16]-work.L[8]*work.d[13]*1;
  residual += temp*temp;
  temp = work.KKT[18]-work.L[11]*work.d[14]*1;
  residual += temp*temp;
  temp = work.KKT[20]-work.L[2]*work.d[15]*1;
  residual += temp*temp;
  temp = work.KKT[22]-work.L[3]*work.d[16]*1;
  residual += temp*temp;
  temp = work.KKT[23]-work.L[0]*work.d[11]*work.L[0]-1*work.d[17]*1;
  residual += temp*temp;
  temp = work.KKT[25]-work.L[1]*work.d[12]*work.L[1]-1*work.d[18]*1;
  residual += temp*temp;
  temp = work.KKT[37]-work.L[8]*work.d[13]*work.L[8]-1*work.d[23]*1-work.L[9]*work.d[21]*work.L[9]-work.L[10]*work.d[22]*work.L[10];
  residual += temp*temp;
  temp = work.KKT[38]-work.L[11]*work.d[14]*work.L[11]-1*work.d[24]*1-work.L[12]*work.d[21]*work.L[12]-work.L[13]*work.d[22]*work.L[13]-work.L[14]*work.d[23]*work.L[14];
  residual += temp*temp;
  temp = work.KKT[27]-work.L[2]*work.d[15]*work.L[2]-1*work.d[19]*1;
  residual += temp*temp;
  temp = work.KKT[29]-work.L[3]*work.d[16]*work.L[3]-1*work.d[20]*1;
  residual += temp*temp;
  temp = work.KKT[24]-1*work.d[17]*work.L[4];
  residual += temp*temp;
  temp = work.KKT[26]-1*work.d[18]*work.L[5];
  residual += temp*temp;
  temp = work.KKT[32]-work.L[9]*work.d[21]*1;
  residual += temp*temp;
  temp = work.KKT[35]-work.L[10]*work.d[22]*1;
  residual += temp*temp;
  temp = work.KKT[33]-work.L[12]*work.d[21]*1;
  residual += temp*temp;
  temp = work.KKT[36]-work.L[13]*work.d[22]*1;
  residual += temp*temp;
  temp = work.KKT[28]-1*work.d[19]*work.L[6];
  residual += temp*temp;
  temp = work.KKT[30]-1*work.d[20]*work.L[7];
  residual += temp*temp;
  return residual;
}
void GensolverFirst::matrix_multiply(double *result, double *source) {
  /* Finds result = A*source. */
  result[0] = work.KKT[31]*source[0]+work.KKT[24]*source[19]+work.KKT[26]*source[20]+work.KKT[32]*source[21]+work.KKT[33]*source[22]+work.KKT[28]*source[23]+work.KKT[30]*source[24];
  result[1] = work.KKT[34]*source[1]+work.KKT[35]*source[21]+work.KKT[36]*source[22];
  result[2] = work.KKT[0]*source[2];
  result[3] = work.KKT[1]*source[3];
  result[4] = work.KKT[2]*source[4];
  result[5] = work.KKT[3]*source[5];
  result[6] = work.KKT[4]*source[6];
  result[7] = work.KKT[5]*source[7];
  result[8] = work.KKT[6]*source[8];
  result[9] = work.KKT[7]*source[9];
  result[10] = work.KKT[8]*source[10];
  result[11] = work.KKT[9]*source[11];
  result[12] = work.KKT[10]*source[12];
  result[13] = work.KKT[11]*source[13]+work.KKT[12]*source[19];
  result[14] = work.KKT[13]*source[14]+work.KKT[14]*source[20];
  result[15] = work.KKT[15]*source[15]+work.KKT[16]*source[21];
  result[16] = work.KKT[17]*source[16]+work.KKT[18]*source[22];
  result[17] = work.KKT[19]*source[17]+work.KKT[20]*source[23];
  result[18] = work.KKT[21]*source[18]+work.KKT[22]*source[24];
  result[19] = work.KKT[12]*source[13]+work.KKT[23]*source[19]+work.KKT[24]*source[0];
  result[20] = work.KKT[14]*source[14]+work.KKT[25]*source[20]+work.KKT[26]*source[0];
  result[21] = work.KKT[16]*source[15]+work.KKT[37]*source[21]+work.KKT[32]*source[0]+work.KKT[35]*source[1];
  result[22] = work.KKT[18]*source[16]+work.KKT[38]*source[22]+work.KKT[33]*source[0]+work.KKT[36]*source[1];
  result[23] = work.KKT[20]*source[17]+work.KKT[27]*source[23]+work.KKT[28]*source[0];
  result[24] = work.KKT[22]*source[18]+work.KKT[29]*source[24]+work.KKT[30]*source[0];
}
double GensolverFirst::check_residual(double *target, double *multiplicand) {
  /* Returns the squared 2-norm of lhs - A*rhs. */
  /* Reuses v to find the residual. */
  int i;
  double residual;
  residual = 0;
  matrix_multiply(work.v, multiplicand);
  for (i = 0; i < 13; i++) {
    residual += (target[i] - work.v[i])*(target[i] - work.v[i]);
  }
  return residual;
}
void GensolverFirst::fill_KKT(void) {
  work.KKT[31] = 2*(params.c2[0]+work.frac_479733190656+work.frac_121674190848*work.quad_877793640448[0]);
  work.KKT[34] = 2*work.frac_479733190656;
  work.KKT[0] = 2*work.frac_121674190848;
  work.KKT[1] = 2*work.frac_121674190848;
  work.KKT[2] = 2*work.frac_121674190848;
  work.KKT[3] = 2*work.frac_121674190848;
  work.KKT[4] = 2*work.frac_121674190848;
  work.KKT[5] = 2*work.frac_121674190848;
  work.KKT[6] = 2*work.frac_121674190848;
  work.KKT[7] = 2*work.frac_121674190848;
  work.KKT[8] = 2*work.frac_121674190848;
  work.KKT[9] = 2*work.frac_121674190848;
  work.KKT[10] = 2*work.frac_121674190848;
  work.KKT[11] = work.s_inv_z[0];
  work.KKT[13] = work.s_inv_z[1];
  work.KKT[15] = work.s_inv_z[2];
  work.KKT[17] = work.s_inv_z[3];
  work.KKT[19] = work.s_inv_z[4];
  work.KKT[21] = work.s_inv_z[5];
  work.KKT[12] = 1;
  work.KKT[14] = 1;
  work.KKT[16] = 1;
  work.KKT[18] = 1;
  work.KKT[20] = 1;
  work.KKT[22] = 1;
  work.KKT[23] = work.block_33[0];
  work.KKT[25] = work.block_33[0];
  work.KKT[37] = work.block_33[0];
  work.KKT[38] = work.block_33[0];
  work.KKT[27] = work.block_33[0];
  work.KKT[29] = work.block_33[0];
  work.KKT[24] = -1;
  work.KKT[26] = 1;
  work.KKT[32] = 1;
  work.KKT[35] = -1;
  work.KKT[33] = -1;
  work.KKT[36] = 1;
  work.KKT[28] = -1;
  work.KKT[30] = 1;
}

//testsolver.cpp
/* Produced by CVXGEN, 2016-05-03 18:52:49 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */

//Vars vars;
//Params params;
//Workspace work;
//Settings settings;
GensolverFirst::GensolverFirst( int numberOfCases, double cost2, double cost1, double cost0, double PgenMax, double PgenMin, double RgenMax, double RgenMin, double Beta, double Gamma, double PgenPrevious ) {
  params.caseNumber = numberOfCases;
  params.contNum[0] = 0.0;
  params.c2[0] = cost2;
  params.c1[0] = cost1;
  params.c0[0] = cost0;
  params.beta[0] = Beta;
  params.PgNu[0] = 0.0;
  params.PgNextNu[0] = 0.0;
  params.gamma[0] = Gamma;
  params.B[0] = 0.0;
  params.D[0] = 0.0;
  params.lambda_1[0] = 0.0;
  params.lambda_2[0] = 0.0;
  params.rho[0] = 1.0;
  params.PgMin[0] = PgenMin;
  params.PgMax[0] = PgenMax;
  params.RgMin[0] = RgenMin;
  params.RgMax[0] = RgenMax;
  params.PgPrev[0] = PgenPrevious;
for ( int i = 0; i <= 10; ++i ) {
  params.ones[i] = 0.0;
  params.Pg_N_init[i] = 0.0;
  params.Pg_N_avg[i] = 0.0;
  params.ug_N[i] = 0.0;
  params.Vg_N_avg[i] = 0.0;
  params.Thetag_N_avg[i] = 0.0;
  params.vg_N[i] = 0.0;
}
  global_seed = 1;
  Piterate = NULL;
  PgNextiterate = NULL;
  Thiterate = NULL;
  /*for ( int i = 0; i <= 3; ++i ) {
   Thiterate[i] = 0.0;
  }*/
  set_defaults();
  setup_indexing();
}
GensolverFirst::~GensolverFirst() {
}
#define NUMTESTS 0
void GensolverFirst::mainsolve(int APPItCount, double gsRho, double Pgprev[], double Pgavrg[], double Pprice[], double Apriceavg[], double Aavg[], double Aprice[], double PgenAPP, double PgenNextAPP, double BAPP, double DAPP, double LambAPP1, double LambAPP2) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif
  set_defaults();
  setup_indexing();
  load_default_data(APPItCount, gsRho, Pgprev, Pgavrg, Pprice, Apriceavg, Aavg, Aprice, PgenAPP, PgenNextAPP, BAPP, DAPP, LambAPP1, LambAPP2);
  /* Solve problem instance for the record. */
  settings.verbose = 0;
  num_iters = solve();
#ifndef ZERO_LIBRARY_MODE
#if (NUMTESTS > 0)
  /* Now solve multiple problem instances for timing purposes. */
  settings.verbose = 0;
  tic();
  for (i = 0; i < NUMTESTS; i++) {
    solve();
  }
  time = tocq();
  printf("Timed %d solves over %.3f seconds.\n", NUMTESTS, time);
  time_per = time / NUMTESTS;
  if (time_per > 1) {
    printf("Actual time taken per solve: %.3g s.\n", time_per);
  } else if (time_per > 1e-3) {
    printf("Actual time taken per solve: %.3g ms.\n", 1e3*time_per);
  } else {
    printf("Actual time taken per solve: %.3g us.\n", 1e6*time_per);
  }
#endif
#endif
}

void GensolverFirst::load_default_data(int APPItCount, double gsRho, double Pgenprev[], double Pgenavg[], double Powerprice[], double Angpriceavg[], double Angavg[], double Angprice[], double PgenAPP, double PgenNextAPP, double BAPP, double DAPP, double LambAPP1, double LambAPP2) {
 if (APPItCount==1) {
  params.rho[0] = gsRho;
  params.PgNu[0] = 0.0;
  params.PgNextNu[0] = 0.0;
  params.B[0] = 0.0;
  params.D[0] = 0.0;
  params.lambda_1[0] = 0.0;
  params.lambda_2[0] = 0.0;
 }
 else {
  params.rho[0] = gsRho;
  params.PgNu[0] = PgenAPP;
  params.PgNextNu[0] = PgenNextAPP;
  params.B[0] = BAPP;
  params.D[0] = DAPP;
  params.lambda_1[0] = LambAPP1;
  params.lambda_2[0] = LambAPP2;
 }
for ( int i = 0; i <= params.caseNumber; ++i ) {
  params.ones[i] = 1.0;
  params.Pg_N_init[i] = Pgenprev[i];
  params.Pg_N_avg[i] = Pgenavg[i];
  params.ug_N[i] = Powerprice[i];
  params.Vg_N_avg[i] = Angpriceavg[i];
  params.Thetag_N_avg[i] = Angavg[i];
  params.vg_N[i] = Angprice[i];
}
  set_defaults();
  setup_indexing();
}

double GensolverFirst::getPSol(void) {
 Piterate = (vars.Pg);
 return *Piterate;
}

double GensolverFirst::getPNextSol(void) {
 Piterate = (vars.PgNext);
 return *PgNextiterate;
}

double GensolverFirst::getObj(void) {
 return (params.c2[0])*(*Piterate)*(*Piterate)+(params.c1[0])*(*Piterate)+(params.c0[0]);
}

double *GensolverFirst::getThetaPtr(void){
Thiterate = (vars.Thetag);
return Thiterate;
}

double GensolverFirst::getRMax() // start of getPMax function
{
	return params.RgMax[0];
} // end of getMax

double GensolverFirst::getRMin() // start of getPMin function
{
	return params.RgMin[0];
} // end of getPMin

double GensolverFirst::getPMax() // start of getPMax function
{
	return params.PgMax[0];
} // end of getMax

double GensolverFirst::getPMin() // start of getPMin function
{
	return params.PgMin[0];
} // end of getPMin

double GensolverFirst::getQuadCoeff() // start of getQuadCoeff function
{
	return params.c2[0]; 
} // end of getQuadCoeff

double GensolverFirst::getLinCoeff() // start of getLinCoeff function
{
	return params.c1[0]; 
} // end of getLinCoeff

double GensolverFirst::getConstCoeff() // start of getConstCoeff function
{
	return params.c0[0]; 
} // end of getConstCoeff
