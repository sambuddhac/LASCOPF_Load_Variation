//solver.cpp
/* Produced by CVXGEN, 2019-02-05 18:01:16 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.c. */
/* Description: Main solver file. */
#include "gensolverFirstBase.h"
#include <array>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <cmath>
using namespace std;
double GensolverFirstBase::eval_gap(void) {
  int i;
  double gap;
  gap = 0;
  for (i = 0; i < 6; i++)
    gap += work.z[i]*work.s[i];
  return gap;
}
void GensolverFirstBase::set_defaults(void) {
  settings.resid_tol = 1e-6;
  settings.eps = 1e-4;
  settings.max_iters = 25;
  settings.refine_steps = 1;
  settings.s_init = 1;
  settings.z_init = 1;
  settings.debug = 0;
  settings.verbose = 1;
  settings.verbose_refinement = 0;
  settings.better_start = 1;
  settings.kkt_reg = 1e-7;
}
void GensolverFirstBase::setup_pointers(void) {
  work.y = work.x + 3;
  work.s = work.x + 3;
  work.z = work.x + 9;
  vars.Pg = work.x + 0;
  vars.PgNext = work.x + 1;
  vars.Thetag = work.x + 2;
}
void GensolverFirstBase::setup_indexing(void) {
  setup_pointers();
}
void GensolverFirstBase::set_start(void) {
  int i;
  for (i = 0; i < 3; i++)
    work.x[i] = 0;
  for (i = 0; i < 0; i++)
    work.y[i] = 0;
  for (i = 0; i < 6; i++)
    work.s[i] = (work.h[i] > 0) ? work.h[i] : settings.s_init;
  for (i = 0; i < 6; i++)
    work.z[i] = settings.z_init;
}
double GensolverFirstBase::eval_objv(void) {
  int i;
  double objv;
  /* Borrow space in work.rhs. */
  multbyP(work.rhs, work.x);
  objv = 0;
  for (i = 0; i < 3; i++)
    objv += work.x[i]*work.rhs[i];
  objv *= 0.5;
  for (i = 0; i < 3; i++)
    objv += work.q[i]*work.x[i];
  objv += params.c0[0]+work.quad_608032620544[0]+work.quad_643747205120[0]+work.quad_311816339456[0]+work.quad_155015135232[0]+work.quad_497818112000[0];
  return objv;
}
void GensolverFirstBase::fillrhs_aff(void) {
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work.rhs;
  r2 = work.rhs + 3;
  r3 = work.rhs + 9;
  r4 = work.rhs + 15;
  /* r1 = -A^Ty - G^Tz - Px - q. */
  multbymAT(r1, work.y);
  multbymGT(work.buffer, work.z);
  for (i = 0; i < 3; i++)
    r1[i] += work.buffer[i];
  multbyP(work.buffer, work.x);
  for (i = 0; i < 3; i++)
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
void GensolverFirstBase::fillrhs_cc(void) {
  int i;
  double *r2;
  double *ds_aff, *dz_aff;
  double mu;
  double alpha;
  double sigma;
  double smu;
  double minval;
  r2 = work.rhs + 3;
  ds_aff = work.lhs_aff + 3;
  dz_aff = work.lhs_aff + 9;
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
  for (i = 0; i < 3; i++)
    work.rhs[i] = 0;
  for (i = 9; i < 15; i++)
    work.rhs[i] = 0;
  for (i = 0; i < 6; i++)
    r2[i] = work.s_inv[i]*(smu - ds_aff[i]*dz_aff[i]);
}
void GensolverFirstBase::refine(double *target, double *var) {
  int i, j;
  double *residual = work.buffer;
  double norm2;
  double *new_var = work.buffer2;
  for (j = 0; j < settings.refine_steps; j++) {
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 15; i++) {
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
    for (i = 0; i < 15; i++) {
      var[i] -= new_var[i];
    }
  }
#ifndef ZERO_LIBRARY_MODE
  if (settings.verbose_refinement) {
    /* Check the residual once more, but only if we're reporting it, since */
    /* it's expensive. */
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 15; i++) {
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
double GensolverFirstBase::calc_ineq_resid_squared(void) {
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
double GensolverFirstBase::calc_eq_resid_squared(void) {
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
void GensolverFirstBase::better_start(void) {
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
  s = work.lhs_aff + 3;
  z = work.lhs_aff + 9;
  y = work.lhs_aff + 15;
  /* Just set x and y as is. */
  for (i = 0; i < 3; i++)
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
void GensolverFirstBase::fillrhs_start(void) {
  /* Fill rhs with (-q, 0, h, b). */
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work.rhs;
  r2 = work.rhs + 3;
  r3 = work.rhs + 9;
  r4 = work.rhs + 15;
  for (i = 0; i < 3; i++)
    r1[i] = -work.q[i];
  for (i = 0; i < 6; i++)
    r2[i] = 0;
  for (i = 0; i < 6; i++)
    r3[i] = work.h[i];
  for (i = 0; i < 0; i++)
    r4[i] = work.b[i];
}
long GensolverFirstBase::solve(void) {
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
    for (i = 0; i < 15; i++)
      work.lhs_aff[i] += work.lhs_cc[i];
    /* Rename aff to reflect its new meaning. */
    dx = work.lhs_aff;
    ds = work.lhs_aff + 3;
    dz = work.lhs_aff + 9;
    dy = work.lhs_aff + 15;
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
    for (i = 0; i < 3; i++)
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
/* Produced by CVXGEN, 2019-02-05 18:01:17 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: util.c. */
/* Description: Common utility file for all cvxgen code. */
void GensolverFirstBase::tic(void) {
  tic_timestart = clock();
}
float GensolverFirstBase::toc(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  printf("time: %8.2f.\n", (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC);
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}
float GensolverFirstBase::tocq(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}
void GensolverFirstBase::printmatrix(char *name, double *A, int m, int n, int sparse) {
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
double GensolverFirstBase::unif(double lower, double upper) {
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
float GensolverFirstBase::ran1(long*idum, int reset) {
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
float GensolverFirstBase::randn_internal(long *idum, int reset) {
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
double GensolverFirstBase::randn(void) {
  return randn_internal(&global_seed, 0);
}
void GensolverFirstBase::reset_rand(void) {
  srand(15);
  global_seed = 1;
  randn_internal(&global_seed, 1);
}

//matrix_support.cpp
/* Produced by CVXGEN, 2019-02-05 18:01:13 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */
void GensolverFirstBase::multbymA(double *lhs, double *rhs) {
}
void GensolverFirstBase::multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
}
void GensolverFirstBase::multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1);
  lhs[1] = -rhs[0]*(1);
  lhs[2] = -rhs[0]*(1)-rhs[1]*(-1);
  lhs[3] = -rhs[0]*(-1)-rhs[1]*(1);
  lhs[4] = -rhs[0]*(-1);
  lhs[5] = -rhs[0]*(1);
}
void GensolverFirstBase::multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1)-rhs[1]*(1)-rhs[2]*(1)-rhs[3]*(-1)-rhs[4]*(-1)-rhs[5]*(1);
  lhs[1] = -rhs[2]*(-1)-rhs[3]*(1);
  lhs[2] = 0;
}
void GensolverFirstBase::multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2*(params.c2[0]+work.frac_479733190656+work.frac_500585390080+work.frac_121674190848));
  lhs[1] = rhs[1]*(2*work.frac_479733190656);
  lhs[2] = rhs[2]*(2*work.frac_121674190848);
}
void GensolverFirstBase::fillq(void) {
  work.q[0] = params.c1[0]-2*params.PgNu[0]*work.frac_479733190656-2*params.PgNuInner[0]*work.frac_500585390080+params.gammaSC[0]*(params.BSC[0]+params.BSC[1]+params.BSC[2]+params.BSC[3]+params.BSC[4]+params.BSC[5]+params.BSC[6]+params.BSC[7]+params.BSC[8]+params.BSC[9]+params.BSC[10]+params.BSC[11]+params.BSC[12]+params.BSC[13]+params.BSC[14]+params.BSC[15]+params.BSC[16]+params.BSC[17]+params.BSC[18]+params.BSC[19]+params.BSC[20]+params.BSC[21]+params.BSC[22]+params.BSC[23]+params.BSC[24]+params.BSC[25]+params.BSC[26]+params.BSC[27]+params.BSC[28]+params.BSC[29]+params.BSC[30]+params.BSC[31]+params.BSC[32]+params.BSC[33]+params.BSC[34]+params.BSC[35]+params.BSC[36]+params.BSC[37]+params.BSC[38]+params.BSC[39]+params.BSC[40]+params.BSC[41]+params.BSC[42]+params.BSC[43]+params.BSC[44]+params.BSC[45]+params.BSC[46]+params.BSC[47]+params.BSC[48]+params.BSC[49]+params.BSC[50]+params.BSC[51]+params.BSC[52]+params.BSC[53]+params.BSC[54]+params.BSC[55]+params.BSC[56]+params.BSC[57]+params.BSC[58]+params.BSC[59]+params.BSC[60]+params.BSC[61]+params.BSC[62]+params.BSC[63]+params.BSC[64]+params.BSC[65]+params.BSC[66]+params.BSC[67]+params.BSC[68]+params.BSC[69]+params.BSC[70]+params.BSC[71]+params.BSC[72]+params.BSC[73]+params.BSC[74]+params.BSC[75]+params.BSC[76]+params.BSC[77]+params.BSC[78]+params.BSC[79]+params.BSC[80]+params.BSC[81]+params.BSC[82]+params.BSC[83]+params.BSC[84]+params.BSC[85]+params.BSC[86]+params.BSC[87]+params.BSC[88]+params.BSC[89]+params.BSC[90]+params.BSC[91]+params.BSC[92]+params.BSC[93]+params.BSC[94]+params.BSC[95]+params.BSC[96]+params.BSC[97]+params.BSC[98]+params.BSC[99]+params.BSC[100]+params.BSC[101]+params.BSC[102]+params.BSC[103]+params.BSC[104]+params.BSC[105]+params.BSC[106]+params.BSC[107]+params.BSC[108]+params.BSC[109]+params.BSC[110]+params.BSC[111]+params.BSC[112]+params.BSC[113]+params.BSC[114]+params.BSC[115]+params.BSC[116]+params.BSC[117]+params.BSC[118]+params.BSC[119]+params.BSC[120]+params.BSC[121]+params.BSC[122]+params.BSC[123]+params.BSC[124]+params.BSC[125]+params.BSC[126]+params.BSC[127]+params.BSC[128]+params.BSC[129]+params.BSC[130]+params.BSC[131]+params.BSC[132]+params.BSC[133]+params.BSC[134]+params.BSC[135]+params.BSC[136]+params.BSC[137]+params.BSC[138]+params.BSC[139]+params.BSC[140]+params.BSC[141]+params.BSC[142]+params.BSC[143]+params.BSC[144]+params.BSC[145]+params.BSC[146]+params.BSC[147]+params.BSC[148]+params.BSC[149]+params.BSC[150]+params.BSC[151]+params.BSC[152]+params.BSC[153]+params.BSC[154]+params.BSC[155]+params.BSC[156]+params.BSC[157]+params.BSC[158]+params.BSC[159]+params.BSC[160]+params.BSC[161]+params.BSC[162]+params.BSC[163]+params.BSC[164]+params.BSC[165]+params.BSC[166]+params.BSC[167]+params.BSC[168]+params.BSC[169]+params.BSC[170]+params.BSC[171]+params.BSC[172]+params.BSC[173]+params.BSC[174]+params.BSC[175]+params.BSC[176]+params.BSC[177]+params.BSC[178]+params.BSC[179]+params.BSC[180]+params.BSC[181]+params.BSC[182]+params.BSC[183]+params.BSC[184]+params.BSC[185]+params.BSC[186]+params.BSC[187]+params.BSC[188]+params.BSC[189]+params.BSC[190]+params.BSC[191]+params.BSC[192]+params.BSC[193]+params.BSC[194]+params.BSC[195]+params.BSC[196]+params.BSC[197]+params.BSC[198]+params.BSC[199]+params.BSC[200]+params.BSC[201]+params.BSC[202]+params.BSC[203]+params.BSC[204]+params.BSC[205]+params.BSC[206]+params.BSC[207]+params.BSC[208]+params.BSC[209]+params.BSC[210]+params.BSC[211]+params.BSC[212]+params.BSC[213]+params.BSC[214]+params.BSC[215]+params.BSC[216]+params.BSC[217]+params.BSC[218]+params.BSC[219]+params.BSC[220]+params.BSC[221]+params.BSC[222]+params.BSC[223]+params.BSC[224]+params.BSC[225]+params.BSC[226]+params.BSC[227]+params.BSC[228]+params.BSC[229]+params.BSC[230]+params.BSC[231]+params.BSC[232]+params.BSC[233]+params.BSC[234]+params.BSC[235]+params.BSC[236]+params.BSC[237]+params.BSC[238]+params.BSC[239]+params.BSC[240]+params.BSC[241]+params.BSC[242]+params.BSC[243]+params.BSC[244]+params.BSC[245]+params.BSC[246]+params.BSC[247]+params.BSC[248]+params.BSC[249]+params.BSC[250]+params.BSC[251]+params.BSC[252]+params.BSC[253]+params.BSC[254]+params.BSC[255]+params.BSC[256]+params.BSC[257]+params.BSC[258]+params.BSC[259]+params.BSC[260]+params.BSC[261]+params.BSC[262]+params.BSC[263]+params.BSC[264]+params.BSC[265]+params.BSC[266]+params.BSC[267]+params.BSC[268]+params.BSC[269]+params.BSC[270]+params.BSC[271]+params.BSC[272]+params.BSC[273]+params.BSC[274]+params.BSC[275]+params.BSC[276]+params.BSC[277]+params.BSC[278]+params.BSC[279]+params.BSC[280]+params.BSC[281]+params.BSC[282]+params.BSC[283]+params.BSC[284]+params.BSC[285]+params.BSC[286]+params.BSC[287]+params.BSC[288]+params.BSC[289]+params.BSC[290]+params.BSC[291]+params.BSC[292]+params.BSC[293]+params.BSC[294]+params.BSC[295]+params.BSC[296]+params.BSC[297]+params.BSC[298]+params.BSC[299]+params.BSC[300]+params.BSC[301]+params.BSC[302]+params.BSC[303]+params.BSC[304]+params.BSC[305]+params.BSC[306]+params.BSC[307]+params.BSC[308]+params.BSC[309]+params.BSC[310]+params.BSC[311]+params.BSC[312]+params.BSC[313]+params.BSC[314]+params.BSC[315]+params.BSC[316]+params.BSC[317]+params.BSC[318]+params.BSC[319]+params.BSC[320]+params.BSC[321]+params.BSC[322]+params.BSC[323]+params.BSC[324]+params.BSC[325]+params.BSC[326]+params.BSC[327]+params.BSC[328]+params.BSC[329]+params.BSC[330]+params.BSC[331]+params.BSC[332]+params.BSC[333]+params.BSC[334]+params.BSC[335]+params.BSC[336]+params.BSC[337]+params.BSC[338]+params.BSC[339]+params.BSC[340]+params.BSC[341]+params.BSC[342]+params.BSC[343]+params.BSC[344]+params.BSC[345]+params.BSC[346]+params.BSC[347]+params.BSC[348]+params.BSC[349]+params.BSC[350]+params.BSC[351]+params.BSC[352]+params.BSC[353]+params.BSC[354]+params.BSC[355]+params.BSC[356]+params.BSC[357]+params.BSC[358]+params.BSC[359]+params.BSC[360]+params.BSC[361]+params.BSC[362]+params.BSC[363]+params.BSC[364]+params.BSC[365]+params.BSC[366]+params.BSC[367]+params.BSC[368]+params.BSC[369]+params.BSC[370]+params.BSC[371]+params.BSC[372]+params.BSC[373]+params.BSC[374]+params.BSC[375]+params.BSC[376]+params.BSC[377]+params.BSC[378]+params.BSC[379]+params.BSC[380]+params.BSC[381]+params.BSC[382]+params.BSC[383]+params.BSC[384]+params.BSC[385]+params.BSC[386]+params.BSC[387]+params.BSC[388]+params.BSC[389]+params.BSC[390]+params.BSC[391]+params.BSC[392]+params.BSC[393]+params.BSC[394]+params.BSC[395]+params.BSC[396]+params.BSC[397]+params.BSC[398]+params.BSC[399]+params.BSC[400]+params.BSC[401]+params.BSC[402]+params.BSC[403]+params.BSC[404]+params.BSC[405]+params.BSC[406]+params.BSC[407]+params.BSC[408]+params.BSC[409]+params.BSC[410]+params.BSC[411]+params.BSC[412]+params.BSC[413]+params.BSC[414]+params.BSC[415]+params.BSC[416]+params.BSC[417]+params.BSC[418]+params.BSC[419]+params.BSC[420]+params.BSC[421]+params.BSC[422]+params.BSC[423]+params.BSC[424]+params.BSC[425]+params.BSC[426]+params.BSC[427]+params.BSC[428]+params.BSC[429]+params.BSC[430]+params.BSC[431]+params.BSC[432]+params.BSC[433]+params.BSC[434]+params.BSC[435]+params.BSC[436]+params.BSC[437]+params.BSC[438]+params.BSC[439]+params.BSC[440]+params.BSC[441]+params.BSC[442]+params.BSC[443]+params.BSC[444]+params.BSC[445]+params.BSC[446]+params.BSC[447]+params.BSC[448]+params.BSC[449]+params.BSC[450]+params.BSC[451]+params.BSC[452]+params.BSC[453]+params.BSC[454]+params.BSC[455]+params.BSC[456]+params.BSC[457]+params.BSC[458]+params.BSC[459]+params.BSC[460]+params.BSC[461]+params.BSC[462]+params.BSC[463]+params.BSC[464]+params.BSC[465]+params.BSC[466]+params.BSC[467]+params.BSC[468]+params.BSC[469]+params.BSC[470]+params.BSC[471]+params.BSC[472]+params.BSC[473]+params.BSC[474]+params.BSC[475]+params.BSC[476]+params.BSC[477]+params.BSC[478]+params.BSC[479]+params.BSC[480]+params.BSC[481]+params.BSC[482]+params.BSC[483]+params.BSC[484]+params.BSC[485]+params.BSC[486]+params.BSC[487]+params.BSC[488]+params.BSC[489]+params.BSC[490]+params.BSC[491]+params.BSC[492]+params.BSC[493]+params.BSC[494]+params.BSC[495]+params.BSC[496]+params.BSC[497]+params.BSC[498]+params.BSC[499])+params.lambda_1SC[0]+params.lambda_1SC[1]+params.lambda_1SC[2]+params.lambda_1SC[3]+params.lambda_1SC[4]+params.lambda_1SC[5]+params.lambda_1SC[6]+params.lambda_1SC[7]+params.lambda_1SC[8]+params.lambda_1SC[9]+params.lambda_1SC[10]+params.lambda_1SC[11]+params.lambda_1SC[12]+params.lambda_1SC[13]+params.lambda_1SC[14]+params.lambda_1SC[15]+params.lambda_1SC[16]+params.lambda_1SC[17]+params.lambda_1SC[18]+params.lambda_1SC[19]+params.lambda_1SC[20]+params.lambda_1SC[21]+params.lambda_1SC[22]+params.lambda_1SC[23]+params.lambda_1SC[24]+params.lambda_1SC[25]+params.lambda_1SC[26]+params.lambda_1SC[27]+params.lambda_1SC[28]+params.lambda_1SC[29]+params.lambda_1SC[30]+params.lambda_1SC[31]+params.lambda_1SC[32]+params.lambda_1SC[33]+params.lambda_1SC[34]+params.lambda_1SC[35]+params.lambda_1SC[36]+params.lambda_1SC[37]+params.lambda_1SC[38]+params.lambda_1SC[39]+params.lambda_1SC[40]+params.lambda_1SC[41]+params.lambda_1SC[42]+params.lambda_1SC[43]+params.lambda_1SC[44]+params.lambda_1SC[45]+params.lambda_1SC[46]+params.lambda_1SC[47]+params.lambda_1SC[48]+params.lambda_1SC[49]+params.lambda_1SC[50]+params.lambda_1SC[51]+params.lambda_1SC[52]+params.lambda_1SC[53]+params.lambda_1SC[54]+params.lambda_1SC[55]+params.lambda_1SC[56]+params.lambda_1SC[57]+params.lambda_1SC[58]+params.lambda_1SC[59]+params.lambda_1SC[60]+params.lambda_1SC[61]+params.lambda_1SC[62]+params.lambda_1SC[63]+params.lambda_1SC[64]+params.lambda_1SC[65]+params.lambda_1SC[66]+params.lambda_1SC[67]+params.lambda_1SC[68]+params.lambda_1SC[69]+params.lambda_1SC[70]+params.lambda_1SC[71]+params.lambda_1SC[72]+params.lambda_1SC[73]+params.lambda_1SC[74]+params.lambda_1SC[75]+params.lambda_1SC[76]+params.lambda_1SC[77]+params.lambda_1SC[78]+params.lambda_1SC[79]+params.lambda_1SC[80]+params.lambda_1SC[81]+params.lambda_1SC[82]+params.lambda_1SC[83]+params.lambda_1SC[84]+params.lambda_1SC[85]+params.lambda_1SC[86]+params.lambda_1SC[87]+params.lambda_1SC[88]+params.lambda_1SC[89]+params.lambda_1SC[90]+params.lambda_1SC[91]+params.lambda_1SC[92]+params.lambda_1SC[93]+params.lambda_1SC[94]+params.lambda_1SC[95]+params.lambda_1SC[96]+params.lambda_1SC[97]+params.lambda_1SC[98]+params.lambda_1SC[99]+params.lambda_1SC[100]+params.lambda_1SC[101]+params.lambda_1SC[102]+params.lambda_1SC[103]+params.lambda_1SC[104]+params.lambda_1SC[105]+params.lambda_1SC[106]+params.lambda_1SC[107]+params.lambda_1SC[108]+params.lambda_1SC[109]+params.lambda_1SC[110]+params.lambda_1SC[111]+params.lambda_1SC[112]+params.lambda_1SC[113]+params.lambda_1SC[114]+params.lambda_1SC[115]+params.lambda_1SC[116]+params.lambda_1SC[117]+params.lambda_1SC[118]+params.lambda_1SC[119]+params.lambda_1SC[120]+params.lambda_1SC[121]+params.lambda_1SC[122]+params.lambda_1SC[123]+params.lambda_1SC[124]+params.lambda_1SC[125]+params.lambda_1SC[126]+params.lambda_1SC[127]+params.lambda_1SC[128]+params.lambda_1SC[129]+params.lambda_1SC[130]+params.lambda_1SC[131]+params.lambda_1SC[132]+params.lambda_1SC[133]+params.lambda_1SC[134]+params.lambda_1SC[135]+params.lambda_1SC[136]+params.lambda_1SC[137]+params.lambda_1SC[138]+params.lambda_1SC[139]+params.lambda_1SC[140]+params.lambda_1SC[141]+params.lambda_1SC[142]+params.lambda_1SC[143]+params.lambda_1SC[144]+params.lambda_1SC[145]+params.lambda_1SC[146]+params.lambda_1SC[147]+params.lambda_1SC[148]+params.lambda_1SC[149]+params.lambda_1SC[150]+params.lambda_1SC[151]+params.lambda_1SC[152]+params.lambda_1SC[153]+params.lambda_1SC[154]+params.lambda_1SC[155]+params.lambda_1SC[156]+params.lambda_1SC[157]+params.lambda_1SC[158]+params.lambda_1SC[159]+params.lambda_1SC[160]+params.lambda_1SC[161]+params.lambda_1SC[162]+params.lambda_1SC[163]+params.lambda_1SC[164]+params.lambda_1SC[165]+params.lambda_1SC[166]+params.lambda_1SC[167]+params.lambda_1SC[168]+params.lambda_1SC[169]+params.lambda_1SC[170]+params.lambda_1SC[171]+params.lambda_1SC[172]+params.lambda_1SC[173]+params.lambda_1SC[174]+params.lambda_1SC[175]+params.lambda_1SC[176]+params.lambda_1SC[177]+params.lambda_1SC[178]+params.lambda_1SC[179]+params.lambda_1SC[180]+params.lambda_1SC[181]+params.lambda_1SC[182]+params.lambda_1SC[183]+params.lambda_1SC[184]+params.lambda_1SC[185]+params.lambda_1SC[186]+params.lambda_1SC[187]+params.lambda_1SC[188]+params.lambda_1SC[189]+params.lambda_1SC[190]+params.lambda_1SC[191]+params.lambda_1SC[192]+params.lambda_1SC[193]+params.lambda_1SC[194]+params.lambda_1SC[195]+params.lambda_1SC[196]+params.lambda_1SC[197]+params.lambda_1SC[198]+params.lambda_1SC[199]+params.lambda_1SC[200]+params.lambda_1SC[201]+params.lambda_1SC[202]+params.lambda_1SC[203]+params.lambda_1SC[204]+params.lambda_1SC[205]+params.lambda_1SC[206]+params.lambda_1SC[207]+params.lambda_1SC[208]+params.lambda_1SC[209]+params.lambda_1SC[210]+params.lambda_1SC[211]+params.lambda_1SC[212]+params.lambda_1SC[213]+params.lambda_1SC[214]+params.lambda_1SC[215]+params.lambda_1SC[216]+params.lambda_1SC[217]+params.lambda_1SC[218]+params.lambda_1SC[219]+params.lambda_1SC[220]+params.lambda_1SC[221]+params.lambda_1SC[222]+params.lambda_1SC[223]+params.lambda_1SC[224]+params.lambda_1SC[225]+params.lambda_1SC[226]+params.lambda_1SC[227]+params.lambda_1SC[228]+params.lambda_1SC[229]+params.lambda_1SC[230]+params.lambda_1SC[231]+params.lambda_1SC[232]+params.lambda_1SC[233]+params.lambda_1SC[234]+params.lambda_1SC[235]+params.lambda_1SC[236]+params.lambda_1SC[237]+params.lambda_1SC[238]+params.lambda_1SC[239]+params.lambda_1SC[240]+params.lambda_1SC[241]+params.lambda_1SC[242]+params.lambda_1SC[243]+params.lambda_1SC[244]+params.lambda_1SC[245]+params.lambda_1SC[246]+params.lambda_1SC[247]+params.lambda_1SC[248]+params.lambda_1SC[249]+params.lambda_1SC[250]+params.lambda_1SC[251]+params.lambda_1SC[252]+params.lambda_1SC[253]+params.lambda_1SC[254]+params.lambda_1SC[255]+params.lambda_1SC[256]+params.lambda_1SC[257]+params.lambda_1SC[258]+params.lambda_1SC[259]+params.lambda_1SC[260]+params.lambda_1SC[261]+params.lambda_1SC[262]+params.lambda_1SC[263]+params.lambda_1SC[264]+params.lambda_1SC[265]+params.lambda_1SC[266]+params.lambda_1SC[267]+params.lambda_1SC[268]+params.lambda_1SC[269]+params.lambda_1SC[270]+params.lambda_1SC[271]+params.lambda_1SC[272]+params.lambda_1SC[273]+params.lambda_1SC[274]+params.lambda_1SC[275]+params.lambda_1SC[276]+params.lambda_1SC[277]+params.lambda_1SC[278]+params.lambda_1SC[279]+params.lambda_1SC[280]+params.lambda_1SC[281]+params.lambda_1SC[282]+params.lambda_1SC[283]+params.lambda_1SC[284]+params.lambda_1SC[285]+params.lambda_1SC[286]+params.lambda_1SC[287]+params.lambda_1SC[288]+params.lambda_1SC[289]+params.lambda_1SC[290]+params.lambda_1SC[291]+params.lambda_1SC[292]+params.lambda_1SC[293]+params.lambda_1SC[294]+params.lambda_1SC[295]+params.lambda_1SC[296]+params.lambda_1SC[297]+params.lambda_1SC[298]+params.lambda_1SC[299]+params.lambda_1SC[300]+params.lambda_1SC[301]+params.lambda_1SC[302]+params.lambda_1SC[303]+params.lambda_1SC[304]+params.lambda_1SC[305]+params.lambda_1SC[306]+params.lambda_1SC[307]+params.lambda_1SC[308]+params.lambda_1SC[309]+params.lambda_1SC[310]+params.lambda_1SC[311]+params.lambda_1SC[312]+params.lambda_1SC[313]+params.lambda_1SC[314]+params.lambda_1SC[315]+params.lambda_1SC[316]+params.lambda_1SC[317]+params.lambda_1SC[318]+params.lambda_1SC[319]+params.lambda_1SC[320]+params.lambda_1SC[321]+params.lambda_1SC[322]+params.lambda_1SC[323]+params.lambda_1SC[324]+params.lambda_1SC[325]+params.lambda_1SC[326]+params.lambda_1SC[327]+params.lambda_1SC[328]+params.lambda_1SC[329]+params.lambda_1SC[330]+params.lambda_1SC[331]+params.lambda_1SC[332]+params.lambda_1SC[333]+params.lambda_1SC[334]+params.lambda_1SC[335]+params.lambda_1SC[336]+params.lambda_1SC[337]+params.lambda_1SC[338]+params.lambda_1SC[339]+params.lambda_1SC[340]+params.lambda_1SC[341]+params.lambda_1SC[342]+params.lambda_1SC[343]+params.lambda_1SC[344]+params.lambda_1SC[345]+params.lambda_1SC[346]+params.lambda_1SC[347]+params.lambda_1SC[348]+params.lambda_1SC[349]+params.lambda_1SC[350]+params.lambda_1SC[351]+params.lambda_1SC[352]+params.lambda_1SC[353]+params.lambda_1SC[354]+params.lambda_1SC[355]+params.lambda_1SC[356]+params.lambda_1SC[357]+params.lambda_1SC[358]+params.lambda_1SC[359]+params.lambda_1SC[360]+params.lambda_1SC[361]+params.lambda_1SC[362]+params.lambda_1SC[363]+params.lambda_1SC[364]+params.lambda_1SC[365]+params.lambda_1SC[366]+params.lambda_1SC[367]+params.lambda_1SC[368]+params.lambda_1SC[369]+params.lambda_1SC[370]+params.lambda_1SC[371]+params.lambda_1SC[372]+params.lambda_1SC[373]+params.lambda_1SC[374]+params.lambda_1SC[375]+params.lambda_1SC[376]+params.lambda_1SC[377]+params.lambda_1SC[378]+params.lambda_1SC[379]+params.lambda_1SC[380]+params.lambda_1SC[381]+params.lambda_1SC[382]+params.lambda_1SC[383]+params.lambda_1SC[384]+params.lambda_1SC[385]+params.lambda_1SC[386]+params.lambda_1SC[387]+params.lambda_1SC[388]+params.lambda_1SC[389]+params.lambda_1SC[390]+params.lambda_1SC[391]+params.lambda_1SC[392]+params.lambda_1SC[393]+params.lambda_1SC[394]+params.lambda_1SC[395]+params.lambda_1SC[396]+params.lambda_1SC[397]+params.lambda_1SC[398]+params.lambda_1SC[399]+params.lambda_1SC[400]+params.lambda_1SC[401]+params.lambda_1SC[402]+params.lambda_1SC[403]+params.lambda_1SC[404]+params.lambda_1SC[405]+params.lambda_1SC[406]+params.lambda_1SC[407]+params.lambda_1SC[408]+params.lambda_1SC[409]+params.lambda_1SC[410]+params.lambda_1SC[411]+params.lambda_1SC[412]+params.lambda_1SC[413]+params.lambda_1SC[414]+params.lambda_1SC[415]+params.lambda_1SC[416]+params.lambda_1SC[417]+params.lambda_1SC[418]+params.lambda_1SC[419]+params.lambda_1SC[420]+params.lambda_1SC[421]+params.lambda_1SC[422]+params.lambda_1SC[423]+params.lambda_1SC[424]+params.lambda_1SC[425]+params.lambda_1SC[426]+params.lambda_1SC[427]+params.lambda_1SC[428]+params.lambda_1SC[429]+params.lambda_1SC[430]+params.lambda_1SC[431]+params.lambda_1SC[432]+params.lambda_1SC[433]+params.lambda_1SC[434]+params.lambda_1SC[435]+params.lambda_1SC[436]+params.lambda_1SC[437]+params.lambda_1SC[438]+params.lambda_1SC[439]+params.lambda_1SC[440]+params.lambda_1SC[441]+params.lambda_1SC[442]+params.lambda_1SC[443]+params.lambda_1SC[444]+params.lambda_1SC[445]+params.lambda_1SC[446]+params.lambda_1SC[447]+params.lambda_1SC[448]+params.lambda_1SC[449]+params.lambda_1SC[450]+params.lambda_1SC[451]+params.lambda_1SC[452]+params.lambda_1SC[453]+params.lambda_1SC[454]+params.lambda_1SC[455]+params.lambda_1SC[456]+params.lambda_1SC[457]+params.lambda_1SC[458]+params.lambda_1SC[459]+params.lambda_1SC[460]+params.lambda_1SC[461]+params.lambda_1SC[462]+params.lambda_1SC[463]+params.lambda_1SC[464]+params.lambda_1SC[465]+params.lambda_1SC[466]+params.lambda_1SC[467]+params.lambda_1SC[468]+params.lambda_1SC[469]+params.lambda_1SC[470]+params.lambda_1SC[471]+params.lambda_1SC[472]+params.lambda_1SC[473]+params.lambda_1SC[474]+params.lambda_1SC[475]+params.lambda_1SC[476]+params.lambda_1SC[477]+params.lambda_1SC[478]+params.lambda_1SC[479]+params.lambda_1SC[480]+params.lambda_1SC[481]+params.lambda_1SC[482]+params.lambda_1SC[483]+params.lambda_1SC[484]+params.lambda_1SC[485]+params.lambda_1SC[486]+params.lambda_1SC[487]+params.lambda_1SC[488]+params.lambda_1SC[489]+params.lambda_1SC[490]+params.lambda_1SC[491]+params.lambda_1SC[492]+params.lambda_1SC[493]+params.lambda_1SC[494]+params.lambda_1SC[495]+params.lambda_1SC[496]+params.lambda_1SC[497]+params.lambda_1SC[498]+params.lambda_1SC[499]+params.gamma[0]*params.B[0]+params.lambda_1[0]+2*(-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*work.frac_121674190848;
  work.q[1] = -2*params.PgNextNu[0]*work.frac_479733190656+params.gamma[0]*params.D[0]+params.lambda_2[0];
  work.q[2] = 2*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])*work.frac_121674190848;
}
void GensolverFirstBase::fillh(void) {
  work.h[0] = -params.PgMin[0];
  work.h[1] = params.PgMax[0];
  work.h[2] = -params.RgMin[0];
  work.h[3] = params.RgMax[0];
  work.h[4] = -(params.RgMin[0]+params.PgPrev[0]);
  work.h[5] = -(-params.PgPrev[0]-params.RgMax[0]);
}
void GensolverFirstBase::fillb(void) {
}
void GensolverFirstBase::pre_ops(void) {
  work.frac_479733190656 = params.beta[0];
  work.frac_479733190656 /= 2;
  work.frac_500585390080 = params.betaInner[0];
  work.frac_500585390080 /= 2;
  work.frac_121674190848 = params.rho[0];
  work.frac_121674190848 /= 2;
  work.quad_608032620544[0] = params.PgNu[0]*work.frac_479733190656*params.PgNu[0];
  work.quad_643747205120[0] = params.PgNextNu[0]*work.frac_479733190656*params.PgNextNu[0];
  work.quad_311816339456[0] = params.PgNuInner[0]*work.frac_500585390080*params.PgNuInner[0];
  work.quad_155015135232[0] = (-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*work.frac_121674190848*(-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0]);
  work.quad_497818112000[0] = (-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])*work.frac_121674190848*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0]);
}

//ldl.cpp
/* Produced by CVXGEN, 2019-02-05 18:01:13 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: ldl.c. */
/* Description: Basic test harness for solver.c. */
/* Be sure to place ldl_solve first, so storage schemes are defined by it. */
void GensolverFirstBase::ldl_solve(double *target, double *var) {
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
  work.v[7] = target[9]-work.L[0]*work.v[1];
  work.v[8] = target[10]-work.L[1]*work.v[2];
  work.v[9] = target[13]-work.L[2]*work.v[5];
  work.v[10] = target[14]-work.L[3]*work.v[6];
  work.v[11] = target[0]-work.L[4]*work.v[7]-work.L[5]*work.v[8]-work.L[6]*work.v[9]-work.L[7]*work.v[10];
  work.v[12] = target[1];
  work.v[13] = target[11]-work.L[8]*work.v[3]-work.L[9]*work.v[11]-work.L[10]*work.v[12];
  work.v[14] = target[12]-work.L[11]*work.v[4]-work.L[12]*work.v[11]-work.L[13]*work.v[12]-work.L[14]*work.v[13];
  /* Diagonal scaling. Assume correctness of work.d_inv. */
  for (i = 0; i < 15; i++)
    work.v[i] *= work.d_inv[i];
  /* Back substitution */
  work.v[13] -= work.L[14]*work.v[14];
  work.v[12] -= work.L[10]*work.v[13]+work.L[13]*work.v[14];
  work.v[11] -= work.L[9]*work.v[13]+work.L[12]*work.v[14];
  work.v[10] -= work.L[7]*work.v[11];
  work.v[9] -= work.L[6]*work.v[11];
  work.v[8] -= work.L[5]*work.v[11];
  work.v[7] -= work.L[4]*work.v[11];
  work.v[6] -= work.L[3]*work.v[10];
  work.v[5] -= work.L[2]*work.v[9];
  work.v[4] -= work.L[11]*work.v[14];
  work.v[3] -= work.L[8]*work.v[13];
  work.v[2] -= work.L[1]*work.v[8];
  work.v[1] -= work.L[0]*work.v[7];
  /* Unpermute the result, from v to var. */
  var[0] = work.v[11];
  var[1] = work.v[12];
  var[2] = work.v[0];
  var[3] = work.v[1];
  var[4] = work.v[2];
  var[5] = work.v[3];
  var[6] = work.v[4];
  var[7] = work.v[5];
  var[8] = work.v[6];
  var[9] = work.v[7];
  var[10] = work.v[8];
  var[11] = work.v[13];
  var[12] = work.v[14];
  var[13] = work.v[9];
  var[14] = work.v[10];
#ifndef ZERO_LIBRARY_MODE
  if (settings.debug) {
    printf("Squared norm for solution is %.8g.\n", check_residual(target, var));
  }
#endif
}
void GensolverFirstBase::ldl_factor(void) {
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
  work.L[0] = (work.KKT[2])*work.d_inv[1];
  work.v[2] = work.KKT[3];
  work.d[2] = work.v[2];
  if (work.d[2] < 0)
    work.d[2] = settings.kkt_reg;
  else
    work.d[2] += settings.kkt_reg;
  work.d_inv[2] = 1/work.d[2];
  work.L[1] = (work.KKT[4])*work.d_inv[2];
  work.v[3] = work.KKT[5];
  work.d[3] = work.v[3];
  if (work.d[3] < 0)
    work.d[3] = settings.kkt_reg;
  else
    work.d[3] += settings.kkt_reg;
  work.d_inv[3] = 1/work.d[3];
  work.L[8] = (work.KKT[6])*work.d_inv[3];
  work.v[4] = work.KKT[7];
  work.d[4] = work.v[4];
  if (work.d[4] < 0)
    work.d[4] = settings.kkt_reg;
  else
    work.d[4] += settings.kkt_reg;
  work.d_inv[4] = 1/work.d[4];
  work.L[11] = (work.KKT[8])*work.d_inv[4];
  work.v[5] = work.KKT[9];
  work.d[5] = work.v[5];
  if (work.d[5] < 0)
    work.d[5] = settings.kkt_reg;
  else
    work.d[5] += settings.kkt_reg;
  work.d_inv[5] = 1/work.d[5];
  work.L[2] = (work.KKT[10])*work.d_inv[5];
  work.v[6] = work.KKT[11];
  work.d[6] = work.v[6];
  if (work.d[6] < 0)
    work.d[6] = settings.kkt_reg;
  else
    work.d[6] += settings.kkt_reg;
  work.d_inv[6] = 1/work.d[6];
  work.L[3] = (work.KKT[12])*work.d_inv[6];
  work.v[1] = work.L[0]*work.d[1];
  work.v[7] = work.KKT[13]-work.L[0]*work.v[1];
  work.d[7] = work.v[7];
  if (work.d[7] > 0)
    work.d[7] = -settings.kkt_reg;
  else
    work.d[7] -= settings.kkt_reg;
  work.d_inv[7] = 1/work.d[7];
  work.L[4] = (work.KKT[14])*work.d_inv[7];
  work.v[2] = work.L[1]*work.d[2];
  work.v[8] = work.KKT[15]-work.L[1]*work.v[2];
  work.d[8] = work.v[8];
  if (work.d[8] > 0)
    work.d[8] = -settings.kkt_reg;
  else
    work.d[8] -= settings.kkt_reg;
  work.d_inv[8] = 1/work.d[8];
  work.L[5] = (work.KKT[16])*work.d_inv[8];
  work.v[5] = work.L[2]*work.d[5];
  work.v[9] = work.KKT[17]-work.L[2]*work.v[5];
  work.d[9] = work.v[9];
  if (work.d[9] > 0)
    work.d[9] = -settings.kkt_reg;
  else
    work.d[9] -= settings.kkt_reg;
  work.d_inv[9] = 1/work.d[9];
  work.L[6] = (work.KKT[18])*work.d_inv[9];
  work.v[6] = work.L[3]*work.d[6];
  work.v[10] = work.KKT[19]-work.L[3]*work.v[6];
  work.d[10] = work.v[10];
  if (work.d[10] > 0)
    work.d[10] = -settings.kkt_reg;
  else
    work.d[10] -= settings.kkt_reg;
  work.d_inv[10] = 1/work.d[10];
  work.L[7] = (work.KKT[20])*work.d_inv[10];
  work.v[7] = work.L[4]*work.d[7];
  work.v[8] = work.L[5]*work.d[8];
  work.v[9] = work.L[6]*work.d[9];
  work.v[10] = work.L[7]*work.d[10];
  work.v[11] = work.KKT[21]-work.L[4]*work.v[7]-work.L[5]*work.v[8]-work.L[6]*work.v[9]-work.L[7]*work.v[10];
  work.d[11] = work.v[11];
  if (work.d[11] < 0)
    work.d[11] = settings.kkt_reg;
  else
    work.d[11] += settings.kkt_reg;
  work.d_inv[11] = 1/work.d[11];
  work.L[9] = (work.KKT[22])*work.d_inv[11];
  work.L[12] = (work.KKT[23])*work.d_inv[11];
  work.v[12] = work.KKT[24];
  work.d[12] = work.v[12];
  if (work.d[12] < 0)
    work.d[12] = settings.kkt_reg;
  else
    work.d[12] += settings.kkt_reg;
  work.d_inv[12] = 1/work.d[12];
  work.L[10] = (work.KKT[25])*work.d_inv[12];
  work.L[13] = (work.KKT[26])*work.d_inv[12];
  work.v[3] = work.L[8]*work.d[3];
  work.v[11] = work.L[9]*work.d[11];
  work.v[12] = work.L[10]*work.d[12];
  work.v[13] = work.KKT[27]-work.L[8]*work.v[3]-work.L[9]*work.v[11]-work.L[10]*work.v[12];
  work.d[13] = work.v[13];
  if (work.d[13] > 0)
    work.d[13] = -settings.kkt_reg;
  else
    work.d[13] -= settings.kkt_reg;
  work.d_inv[13] = 1/work.d[13];
  work.L[14] = (-work.L[12]*work.v[11]-work.L[13]*work.v[12])*work.d_inv[13];
  work.v[4] = work.L[11]*work.d[4];
  work.v[11] = work.L[12]*work.d[11];
  work.v[12] = work.L[13]*work.d[12];
  work.v[13] = work.L[14]*work.d[13];
  work.v[14] = work.KKT[28]-work.L[11]*work.v[4]-work.L[12]*work.v[11]-work.L[13]*work.v[12]-work.L[14]*work.v[13];
  work.d[14] = work.v[14];
  if (work.d[14] > 0)
    work.d[14] = -settings.kkt_reg;
  else
    work.d[14] -= settings.kkt_reg;
  work.d_inv[14] = 1/work.d[14];
#ifndef ZERO_LIBRARY_MODE
  if (settings.debug) {
    printf("Squared Frobenius for factorization is %.8g.\n", check_factorization());
  }
#endif
}
double GensolverFirstBase::check_factorization(void) {
  /* Returns the squared Frobenius norm of A - L*D*L'. */
  double temp, residual;
  /* Only check the lower triangle. */
  residual = 0;
  temp = work.KKT[21]-1*work.d[11]*1-work.L[4]*work.d[7]*work.L[4]-work.L[5]*work.d[8]*work.L[5]-work.L[6]*work.d[9]*work.L[6]-work.L[7]*work.d[10]*work.L[7];
  residual += temp*temp;
  temp = work.KKT[24]-1*work.d[12]*1;
  residual += temp*temp;
  temp = work.KKT[0]-1*work.d[0]*1;
  residual += temp*temp;
  temp = work.KKT[1]-1*work.d[1]*1;
  residual += temp*temp;
  temp = work.KKT[3]-1*work.d[2]*1;
  residual += temp*temp;
  temp = work.KKT[5]-1*work.d[3]*1;
  residual += temp*temp;
  temp = work.KKT[7]-1*work.d[4]*1;
  residual += temp*temp;
  temp = work.KKT[9]-1*work.d[5]*1;
  residual += temp*temp;
  temp = work.KKT[11]-1*work.d[6]*1;
  residual += temp*temp;
  temp = work.KKT[2]-work.L[0]*work.d[1]*1;
  residual += temp*temp;
  temp = work.KKT[4]-work.L[1]*work.d[2]*1;
  residual += temp*temp;
  temp = work.KKT[6]-work.L[8]*work.d[3]*1;
  residual += temp*temp;
  temp = work.KKT[8]-work.L[11]*work.d[4]*1;
  residual += temp*temp;
  temp = work.KKT[10]-work.L[2]*work.d[5]*1;
  residual += temp*temp;
  temp = work.KKT[12]-work.L[3]*work.d[6]*1;
  residual += temp*temp;
  temp = work.KKT[13]-work.L[0]*work.d[1]*work.L[0]-1*work.d[7]*1;
  residual += temp*temp;
  temp = work.KKT[15]-work.L[1]*work.d[2]*work.L[1]-1*work.d[8]*1;
  residual += temp*temp;
  temp = work.KKT[27]-work.L[8]*work.d[3]*work.L[8]-1*work.d[13]*1-work.L[9]*work.d[11]*work.L[9]-work.L[10]*work.d[12]*work.L[10];
  residual += temp*temp;
  temp = work.KKT[28]-work.L[11]*work.d[4]*work.L[11]-1*work.d[14]*1-work.L[12]*work.d[11]*work.L[12]-work.L[13]*work.d[12]*work.L[13]-work.L[14]*work.d[13]*work.L[14];
  residual += temp*temp;
  temp = work.KKT[17]-work.L[2]*work.d[5]*work.L[2]-1*work.d[9]*1;
  residual += temp*temp;
  temp = work.KKT[19]-work.L[3]*work.d[6]*work.L[3]-1*work.d[10]*1;
  residual += temp*temp;
  temp = work.KKT[14]-1*work.d[7]*work.L[4];
  residual += temp*temp;
  temp = work.KKT[16]-1*work.d[8]*work.L[5];
  residual += temp*temp;
  temp = work.KKT[22]-work.L[9]*work.d[11]*1;
  residual += temp*temp;
  temp = work.KKT[25]-work.L[10]*work.d[12]*1;
  residual += temp*temp;
  temp = work.KKT[23]-work.L[12]*work.d[11]*1;
  residual += temp*temp;
  temp = work.KKT[26]-work.L[13]*work.d[12]*1;
  residual += temp*temp;
  temp = work.KKT[18]-1*work.d[9]*work.L[6];
  residual += temp*temp;
  temp = work.KKT[20]-1*work.d[10]*work.L[7];
  residual += temp*temp;
  return residual;
}
void GensolverFirstBase::matrix_multiply(double *result, double *source) {
  /* Finds result = A*source. */
  result[0] = work.KKT[21]*source[0]+work.KKT[14]*source[9]+work.KKT[16]*source[10]+work.KKT[22]*source[11]+work.KKT[23]*source[12]+work.KKT[18]*source[13]+work.KKT[20]*source[14];
  result[1] = work.KKT[24]*source[1]+work.KKT[25]*source[11]+work.KKT[26]*source[12];
  result[2] = work.KKT[0]*source[2];
  result[3] = work.KKT[1]*source[3]+work.KKT[2]*source[9];
  result[4] = work.KKT[3]*source[4]+work.KKT[4]*source[10];
  result[5] = work.KKT[5]*source[5]+work.KKT[6]*source[11];
  result[6] = work.KKT[7]*source[6]+work.KKT[8]*source[12];
  result[7] = work.KKT[9]*source[7]+work.KKT[10]*source[13];
  result[8] = work.KKT[11]*source[8]+work.KKT[12]*source[14];
  result[9] = work.KKT[2]*source[3]+work.KKT[13]*source[9]+work.KKT[14]*source[0];
  result[10] = work.KKT[4]*source[4]+work.KKT[15]*source[10]+work.KKT[16]*source[0];
  result[11] = work.KKT[6]*source[5]+work.KKT[27]*source[11]+work.KKT[22]*source[0]+work.KKT[25]*source[1];
  result[12] = work.KKT[8]*source[6]+work.KKT[28]*source[12]+work.KKT[23]*source[0]+work.KKT[26]*source[1];
  result[13] = work.KKT[10]*source[7]+work.KKT[17]*source[13]+work.KKT[18]*source[0];
  result[14] = work.KKT[12]*source[8]+work.KKT[19]*source[14]+work.KKT[20]*source[0];
}
double GensolverFirstBase::check_residual(double *target, double *multiplicand) {
  /* Returns the squared 2-norm of lhs - A*rhs. */
  /* Reuses v to find the residual. */
  int i;
  double residual;
  residual = 0;
  matrix_multiply(work.v, multiplicand);
  for (i = 0; i < 3; i++) {
    residual += (target[i] - work.v[i])*(target[i] - work.v[i]);
  }
  return residual;
}
void GensolverFirstBase::fill_KKT(void) {
  work.KKT[21] = 2*(params.c2[0]+work.frac_479733190656+work.frac_500585390080+work.frac_121674190848);
  work.KKT[24] = 2*work.frac_479733190656;
  work.KKT[0] = 2*work.frac_121674190848;
  work.KKT[1] = work.s_inv_z[0];
  work.KKT[3] = work.s_inv_z[1];
  work.KKT[5] = work.s_inv_z[2];
  work.KKT[7] = work.s_inv_z[3];
  work.KKT[9] = work.s_inv_z[4];
  work.KKT[11] = work.s_inv_z[5];
  work.KKT[2] = 1;
  work.KKT[4] = 1;
  work.KKT[6] = 1;
  work.KKT[8] = 1;
  work.KKT[10] = 1;
  work.KKT[12] = 1;
  work.KKT[13] = work.block_33[0];
  work.KKT[15] = work.block_33[0];
  work.KKT[27] = work.block_33[0];
  work.KKT[28] = work.block_33[0];
  work.KKT[17] = work.block_33[0];
  work.KKT[19] = work.block_33[0];
  work.KKT[14] = -1;
  work.KKT[16] = 1;
  work.KKT[22] = 1;
  work.KKT[25] = -1;
  work.KKT[23] = -1;
  work.KKT[26] = 1;
  work.KKT[18] = -1;
  work.KKT[20] = 1;
}

//testsolver.cpp
/* Produced by CVXGEN, 2019-02-05 18:01:17 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */
GensolverFirstBase::GensolverFirstBase( int numberOfCases, double cost2, double cost1, double cost0, double PgenMax, double PgenMin, double RgenMax, double RgenMin, double Beta, double InnerBeta, double Gamma, double innerGamma, double PgenPrevious ) {
  params.caseNumber = numberOfCases;
  params.c2[0] = cost2;
  params.c1[0] = cost1;
  params.c0[0] = cost0;
  params.beta[0] = Beta;
  params.PgNu[0] = 0.0;
  params.PgNextNu[0] = 0.0;
  params.betaInner[0] = InnerBeta;
  params.PgNuInner[0] = 0.0;
  params.gamma[0] = Gamma;
  params.gammaSC[0] = innerGamma;
  params.B[0] = 0.0;
  params.D[0] = 0.0;
  params.lambda_1[0] = 0.0;
  params.lambda_2[0] = 0.0;
  for ( int i = 0; i < 500; ++i ) {
   params.BSC[i] = 0.0;
   params.lambda_1SC[i] = 0.0;
  }
  params.rho[0] = 1.0;
  params.PgMin[0] = PgenMin;
  params.PgMax[0] = PgenMax;
  params.RgMin[0] = RgenMin;
  params.RgMax[0] = RgenMax;
  params.PgPrev[0] = PgenPrevious;
  params.Pg_N_init[0] = 0.0;
  params.Pg_N_avg[0] = 0.0;
  params.ug_N[0] = 0.0;
  params.Vg_N_avg[0] = 0.0;
  params.Thetag_N_avg[0] = 0.0;
  params.vg_N[0] = 0.0;
  global_seed = 1;
  Piterate = NULL;
  PgNextiterate = NULL;
  Thiterate = NULL;
  set_defaults();
  setup_indexing();
}
GensolverFirstBase::~GensolverFirstBase() {
}
#define NUMTESTS 0
void GensolverFirstBase::mainsolve(int outerAPPIt, int APPItCount, double gsRho, double Pgenprev, double Pgenavg, double Powerprice, double Angpriceavg, double Angavg, double Angprice, double PgenAPP, double PgenAPPInner, double PgenNextAPP, double BAPP, double DAPP, double LambAPP1, double LambAPP2, double BSC[], double lambdaSC[]) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif
  set_defaults();
  setup_indexing();
  load_default_data(outerAPPIt, APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgenNextAPP, BAPP, DAPP, LambAPP1, LambAPP2, BSC, lambdaSC);
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
void GensolverFirstBase::load_default_data(int outerAPPIt, int APPItCount, double gsRho, double Pgenprev, double Pgenavg, double Powerprice, double Angpriceavg, double Angavg, double Angprice, double PgenAPP, double PgenAPPInner, double PgenNextAPP, double BAPP, double DAPP, double LambAPP1, double LambAPP2, double BSC[], double lambdaSC[]) {
  if (outerAPPIt==1) {
   if (APPItCount==1) {
    params.rho[0] = gsRho;
    params.PgNu[0] = 0.0;
    params.PgNextNu[0] = 0.0;
    params.PgNuInner[0] = PgenAPPInner;
    params.B[0] = 0.0;
    params.D[0] = 0.0;
    params.lambda_1[0] = 0.0;
    params.lambda_2[0] = 0.0;
    for (int i=0; i<params.caseNumber; ++i) {
     params.BSC[i] = 0.0;
     params.lambda_1SC[i] = 0.0;
    }
   }
   else {
    params.rho[0] = gsRho;
    params.PgNu[0] = PgenAPP;
    params.PgNextNu[0] = 0.0;
    params.PgNuInner[0] = PgenAPPInner;
    params.B[0] = 0.0;
    params.D[0] = 0.0;
    params.lambda_1[0] = 0.0;
    params.lambda_2[0] = 0.0;
    for (int i=0; i<params.caseNumber; ++i) {
     params.BSC[i] = BSC[i];
     params.lambda_1SC[i] = lambdaSC[i];
    }
   }
  }
  else {
   params.rho[0] = gsRho;
   params.PgNu[0] = PgenAPP;
   params.PgNextNu[0] = PgenNextAPP;
   params.PgNuInner[0] = PgenAPPInner;
   params.B[0] = BAPP;
   params.D[0] = DAPP;
   params.lambda_1[0] = LambAPP1;
   params.lambda_2[0] = LambAPP2;
   for ( int i = 0; i < params.caseNumber; ++i ) {
    params.BSC[i] = BSC[i];
    params.lambda_1SC[i] = lambdaSC[i];
   }
  }
  params.Pg_N_init[0] = Pgenprev;
  params.Pg_N_avg[0] = Pgenavg;
  params.ug_N[0] = Powerprice;
  params.Vg_N_avg[0] = Angpriceavg;
  params.Thetag_N_avg[0] = Angavg;
  params.vg_N[0] = Angprice;
  set_defaults();
  setup_indexing();
}

double GensolverFirstBase::getPgPrev(void) {
 return params.PgPrev[0];
}

double GensolverFirstBase::getPSol(void) {
 Piterate = (vars.Pg);
 return *Piterate;
}

double GensolverFirstBase::getPNextSol(void) {
 PgNextiterate = (vars.PgNext);
 return *PgNextiterate;
}

double GensolverFirstBase::getObj(void) {
 return (params.c2[0])*(*Piterate)*(*Piterate)+(params.c1[0])*(*Piterate)+(params.c0[0]);
}

double *GensolverFirstBase::getThetaPtr(void){
Thiterate = (vars.Thetag);
return Thiterate;
}

double GensolverFirstBase::getRMax() // start of getPMax function
{
	return params.RgMax[0];
} // end of getMax

double GensolverFirstBase::getRMin() // start of getPMin function
{
	return params.RgMin[0];
} // end of getPMin

double GensolverFirstBase::getPMax() // start of getPMax function
{
	return params.PgMax[0];
} // end of getMax

double GensolverFirstBase::getPMin() // start of getPMin function
{
	return params.PgMin[0];
} // end of getPMin

double GensolverFirstBase::getQuadCoeff() // start of getQuadCoeff function
{
	return params.c2[0]; 
} // end of getQuadCoeff

double GensolverFirstBase::getLinCoeff() // start of getLinCoeff function
{
	return params.c1[0]; 
} // end of getLinCoeff

double GensolverFirstBase::getConstCoeff() // start of getConstCoeff function
{
	return params.c0[0]; 
} // end of getConstCoeff

double GensolverFirstBase::getBeta(){return params.beta[0];}
double GensolverFirstBase::getGamma(){return params.gamma[0];}
double GensolverFirstBase::getIntGamma(){return params.gammaSC[0];}
