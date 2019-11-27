/* Produced by CVXGEN, 2016-05-03 19:26:42 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */
#include "solver.h"
void multbymA(double *lhs, double *rhs) {
}
void multbymAT(double *lhs, double *rhs) {
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
  lhs[13] = 0;
}
void multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1);
  lhs[1] = -rhs[0]*(1);
  lhs[2] = -rhs[0]*(1)-rhs[1]*(-1);
  lhs[3] = -rhs[0]*(-1)-rhs[1]*(1);
  lhs[4] = -rhs[0]*(-1)-rhs[2]*(1);
  lhs[5] = -rhs[0]*(1)-rhs[2]*(-1);
}
void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1)-rhs[1]*(1)-rhs[2]*(1)-rhs[3]*(-1)-rhs[4]*(-1)-rhs[5]*(1);
  lhs[1] = -rhs[2]*(-1)-rhs[3]*(1);
  lhs[2] = -rhs[4]*(1)-rhs[5]*(-1);
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
  lhs[13] = 0;
}
void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2*(params.c2[0]+work.frac_479733190656+work.frac_121674190848*work.quad_877793640448[0]));
  lhs[1] = rhs[1]*(2*work.frac_479733190656);
  lhs[2] = rhs[2]*(2*work.frac_479733190656);
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
  lhs[13] = rhs[13]*(2*work.frac_121674190848);
}
void fillq(void) {
  work.q[0] = params.c1[0]-2*params.PgNu[0]*work.frac_479733190656+params.gamma[0]*params.B[0]+params.lambda_1[0]-params.lambda_4[0]+2*work.frac_121674190848*((-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*params.ones[0]+(-params.Pg_N_init[1]+params.Pg_N_avg[1]+params.ug_N[1])*params.ones[1]+(-params.Pg_N_init[2]+params.Pg_N_avg[2]+params.ug_N[2])*params.ones[2]+(-params.Pg_N_init[3]+params.Pg_N_avg[3]+params.ug_N[3])*params.ones[3]+(-params.Pg_N_init[4]+params.Pg_N_avg[4]+params.ug_N[4])*params.ones[4]+(-params.Pg_N_init[5]+params.Pg_N_avg[5]+params.ug_N[5])*params.ones[5]+(-params.Pg_N_init[6]+params.Pg_N_avg[6]+params.ug_N[6])*params.ones[6]+(-params.Pg_N_init[7]+params.Pg_N_avg[7]+params.ug_N[7])*params.ones[7]+(-params.Pg_N_init[8]+params.Pg_N_avg[8]+params.ug_N[8])*params.ones[8]+(-params.Pg_N_init[9]+params.Pg_N_avg[9]+params.ug_N[9])*params.ones[9]+(-params.Pg_N_init[10]+params.Pg_N_avg[10]+params.ug_N[10])*params.ones[10]);
  work.q[1] = -2*params.PgNextNu[0]*work.frac_479733190656+params.gamma[0]*params.D[0]+params.lambda_2[0];
  work.q[2] = -2*params.PgPrevNu[0]*work.frac_479733190656+params.gamma[0]*params.A[0]-params.lambda_3[0];
  work.q[3] = 2*work.frac_121674190848*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0]);
  work.q[4] = 2*work.frac_121674190848*(-params.Vg_N_avg[1]-params.Thetag_N_avg[1]+params.vg_N[1]);
  work.q[5] = 2*work.frac_121674190848*(-params.Vg_N_avg[2]-params.Thetag_N_avg[2]+params.vg_N[2]);
  work.q[6] = 2*work.frac_121674190848*(-params.Vg_N_avg[3]-params.Thetag_N_avg[3]+params.vg_N[3]);
  work.q[7] = 2*work.frac_121674190848*(-params.Vg_N_avg[4]-params.Thetag_N_avg[4]+params.vg_N[4]);
  work.q[8] = 2*work.frac_121674190848*(-params.Vg_N_avg[5]-params.Thetag_N_avg[5]+params.vg_N[5]);
  work.q[9] = 2*work.frac_121674190848*(-params.Vg_N_avg[6]-params.Thetag_N_avg[6]+params.vg_N[6]);
  work.q[10] = 2*work.frac_121674190848*(-params.Vg_N_avg[7]-params.Thetag_N_avg[7]+params.vg_N[7]);
  work.q[11] = 2*work.frac_121674190848*(-params.Vg_N_avg[8]-params.Thetag_N_avg[8]+params.vg_N[8]);
  work.q[12] = 2*work.frac_121674190848*(-params.Vg_N_avg[9]-params.Thetag_N_avg[9]+params.vg_N[9]);
  work.q[13] = 2*work.frac_121674190848*(-params.Vg_N_avg[10]-params.Thetag_N_avg[10]+params.vg_N[10]);
}
void fillh(void) {
  work.h[0] = -params.PgMin[0];
  work.h[1] = params.PgMax[0];
  work.h[2] = -params.RgMin[0];
  work.h[3] = params.RgMax[0];
  work.h[4] = -params.RgMin[0];
  work.h[5] = params.RgMax[0];
}
void fillb(void) {
}
void pre_ops(void) {
  work.frac_479733190656 = params.beta[0];
  work.frac_479733190656 /= 2;
  work.frac_121674190848 = params.rho[0];
  work.frac_121674190848 /= 2;
  work.quad_877793640448[0] = params.ones[0]*params.ones[0]+params.ones[1]*params.ones[1]+params.ones[2]*params.ones[2]+params.ones[3]*params.ones[3]+params.ones[4]*params.ones[4]+params.ones[5]*params.ones[5]+params.ones[6]*params.ones[6]+params.ones[7]*params.ones[7]+params.ones[8]*params.ones[8]+params.ones[9]*params.ones[9]+params.ones[10]*params.ones[10];
  work.quad_399913824256[0] = params.PgPrevNu[0]*work.frac_479733190656*params.PgPrevNu[0];
  work.quad_608032620544[0] = params.PgNu[0]*work.frac_479733190656*params.PgNu[0];
  work.quad_643747205120[0] = params.PgNextNu[0]*work.frac_479733190656*params.PgNextNu[0];
  work.quad_324888432640[0] = ((-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*(-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])+(-params.Pg_N_init[1]+params.Pg_N_avg[1]+params.ug_N[1])*(-params.Pg_N_init[1]+params.Pg_N_avg[1]+params.ug_N[1])+(-params.Pg_N_init[2]+params.Pg_N_avg[2]+params.ug_N[2])*(-params.Pg_N_init[2]+params.Pg_N_avg[2]+params.ug_N[2])+(-params.Pg_N_init[3]+params.Pg_N_avg[3]+params.ug_N[3])*(-params.Pg_N_init[3]+params.Pg_N_avg[3]+params.ug_N[3])+(-params.Pg_N_init[4]+params.Pg_N_avg[4]+params.ug_N[4])*(-params.Pg_N_init[4]+params.Pg_N_avg[4]+params.ug_N[4])+(-params.Pg_N_init[5]+params.Pg_N_avg[5]+params.ug_N[5])*(-params.Pg_N_init[5]+params.Pg_N_avg[5]+params.ug_N[5])+(-params.Pg_N_init[6]+params.Pg_N_avg[6]+params.ug_N[6])*(-params.Pg_N_init[6]+params.Pg_N_avg[6]+params.ug_N[6])+(-params.Pg_N_init[7]+params.Pg_N_avg[7]+params.ug_N[7])*(-params.Pg_N_init[7]+params.Pg_N_avg[7]+params.ug_N[7])+(-params.Pg_N_init[8]+params.Pg_N_avg[8]+params.ug_N[8])*(-params.Pg_N_init[8]+params.Pg_N_avg[8]+params.ug_N[8])+(-params.Pg_N_init[9]+params.Pg_N_avg[9]+params.ug_N[9])*(-params.Pg_N_init[9]+params.Pg_N_avg[9]+params.ug_N[9])+(-params.Pg_N_init[10]+params.Pg_N_avg[10]+params.ug_N[10])*(-params.Pg_N_init[10]+params.Pg_N_avg[10]+params.ug_N[10]));
  work.quad_331507691520[0] = ((-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])+(-params.Vg_N_avg[1]-params.Thetag_N_avg[1]+params.vg_N[1])*(-params.Vg_N_avg[1]-params.Thetag_N_avg[1]+params.vg_N[1])+(-params.Vg_N_avg[2]-params.Thetag_N_avg[2]+params.vg_N[2])*(-params.Vg_N_avg[2]-params.Thetag_N_avg[2]+params.vg_N[2])+(-params.Vg_N_avg[3]-params.Thetag_N_avg[3]+params.vg_N[3])*(-params.Vg_N_avg[3]-params.Thetag_N_avg[3]+params.vg_N[3])+(-params.Vg_N_avg[4]-params.Thetag_N_avg[4]+params.vg_N[4])*(-params.Vg_N_avg[4]-params.Thetag_N_avg[4]+params.vg_N[4])+(-params.Vg_N_avg[5]-params.Thetag_N_avg[5]+params.vg_N[5])*(-params.Vg_N_avg[5]-params.Thetag_N_avg[5]+params.vg_N[5])+(-params.Vg_N_avg[6]-params.Thetag_N_avg[6]+params.vg_N[6])*(-params.Vg_N_avg[6]-params.Thetag_N_avg[6]+params.vg_N[6])+(-params.Vg_N_avg[7]-params.Thetag_N_avg[7]+params.vg_N[7])*(-params.Vg_N_avg[7]-params.Thetag_N_avg[7]+params.vg_N[7])+(-params.Vg_N_avg[8]-params.Thetag_N_avg[8]+params.vg_N[8])*(-params.Vg_N_avg[8]-params.Thetag_N_avg[8]+params.vg_N[8])+(-params.Vg_N_avg[9]-params.Thetag_N_avg[9]+params.vg_N[9])*(-params.Vg_N_avg[9]-params.Thetag_N_avg[9]+params.vg_N[9])+(-params.Vg_N_avg[10]-params.Thetag_N_avg[10]+params.vg_N[10])*(-params.Vg_N_avg[10]-params.Thetag_N_avg[10]+params.vg_N[10]));
}
