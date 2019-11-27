/* Produced by CVXGEN, 2016-05-03 19:26:42 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */
#include "solver.h"
Vars vars;
Params params;
Workspace work;
Settings settings;
#define NUMTESTS 0
int main(int argc, char **argv) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif
  set_defaults();
  setup_indexing();
  load_default_data();
  /* Solve problem instance for the record. */
  settings.verbose = 1;
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
  return 0;
}
void load_default_data(void) {
  params.c2[0] = 1.101595805149151;
  params.c1[0] = 1.4162956452362097;
  params.c0[0] = 0.5818094778258887;
  params.beta[0] = 1.021655210395326;
  params.PgPrevNu[0] = 1.5717878173906188;
  params.PgNu[0] = 1.5851723557337523;
  params.PgNextNu[0] = -1.497658758144655;
  params.gamma[0] = 0.4144857562763735;
  params.A[0] = -1.7941311867966805;
  params.B[0] = -0.23676062539745413;
  params.D[0] = -1.8804951564857322;
  params.lambda_1[0] = -0.17266710242115568;
  params.lambda_2[0] = 0.596576190459043;
  params.lambda_3[0] = -0.8860508694080989;
  params.lambda_4[0] = 0.7050196079205251;
  params.rho[0] = 1.1817256348327017;
  params.ones[0] = -1.9040724704913385;
  params.ones[1] = 0.23541635196352795;
  params.ones[2] = -0.9629902123701384;
  params.ones[3] = -0.3395952119597214;
  params.ones[4] = -0.865899672914725;
  params.ones[5] = 0.7725516732519853;
  params.ones[6] = -0.23818512931704205;
  params.ones[7] = -1.372529046100147;
  params.ones[8] = 0.17859607212737894;
  params.ones[9] = 1.1212590580454682;
  params.ones[10] = -0.774545870495281;
  params.Pg_N_init[0] = -1.1121684642712744;
  params.Pg_N_init[1] = -0.44811496977740495;
  params.Pg_N_init[2] = 1.7455345994417217;
  params.Pg_N_init[3] = 1.9039816898917352;
  params.Pg_N_init[4] = 0.6895347036512547;
  params.Pg_N_init[5] = 1.6113364341535923;
  params.Pg_N_init[6] = 1.383003485172717;
  params.Pg_N_init[7] = -0.48802383468444344;
  params.Pg_N_init[8] = -1.631131964513103;
  params.Pg_N_init[9] = 0.6136436100941447;
  params.Pg_N_init[10] = 0.2313630495538037;
  params.Pg_N_avg[0] = -0.5537409477496875;
  params.Pg_N_avg[1] = -1.0997819806406723;
  params.Pg_N_avg[2] = -0.3739203344950055;
  params.Pg_N_avg[3] = -0.12423900520332376;
  params.Pg_N_avg[4] = -0.923057686995755;
  params.Pg_N_avg[5] = -0.8328289030982696;
  params.Pg_N_avg[6] = -0.16925440270808823;
  params.Pg_N_avg[7] = 1.442135651787706;
  params.Pg_N_avg[8] = 0.34501161787128565;
  params.Pg_N_avg[9] = -0.8660485502711608;
  params.Pg_N_avg[10] = -0.8880899735055947;
  params.ug_N[0] = -0.1815116979122129;
  params.ug_N[1] = -1.17835862158005;
  params.ug_N[2] = -1.1944851558277074;
  params.ug_N[3] = 0.05614023926976763;
  params.ug_N[4] = -1.6510825248767813;
  params.ug_N[5] = -0.06565787059365391;
  params.ug_N[6] = -0.5512951504486665;
  params.ug_N[7] = 0.8307464872626844;
  params.ug_N[8] = 0.9869848924080182;
  params.ug_N[9] = 0.7643716874230573;
  params.ug_N[10] = 0.7567216550196565;
  params.Vg_N_avg[0] = -0.5055995034042868;
  params.Vg_N_avg[1] = 0.6725392189410702;
  params.Vg_N_avg[2] = -0.6406053441727284;
  params.Vg_N_avg[3] = 0.29117547947550015;
  params.Vg_N_avg[4] = -0.6967713677405021;
  params.Vg_N_avg[5] = -0.21941980294587182;
  params.Vg_N_avg[6] = -1.753884276680243;
  params.Vg_N_avg[7] = -1.0292983112626475;
  params.Vg_N_avg[8] = 1.8864104246942706;
  params.Vg_N_avg[9] = -1.077663182579704;
  params.Vg_N_avg[10] = 0.7659100437893209;
  params.Thetag_N_avg[0] = 0.6019074328549583;
  params.Thetag_N_avg[1] = 0.8957565577499285;
  params.Thetag_N_avg[2] = -0.09964555746227477;
  params.Thetag_N_avg[3] = 0.38665509840745127;
  params.Thetag_N_avg[4] = -1.7321223042686946;
  params.Thetag_N_avg[5] = -1.7097514487110663;
  params.Thetag_N_avg[6] = -1.2040958948116867;
  params.Thetag_N_avg[7] = -1.3925560119658358;
  params.Thetag_N_avg[8] = -1.5995826216742213;
  params.Thetag_N_avg[9] = -1.4828245415645833;
  params.Thetag_N_avg[10] = 0.21311092723061398;
  params.vg_N[0] = -1.248740700304487;
  params.vg_N[1] = 1.808404972124833;
  params.vg_N[2] = 0.7264471152297065;
  params.vg_N[3] = 0.16407869343908477;
  params.vg_N[4] = 0.8287224032315907;
  params.vg_N[5] = -0.9444533161899464;
  params.vg_N[6] = 1.7069027370149112;
  params.vg_N[7] = 1.3567722311998827;
  params.vg_N[8] = 0.9052779937121489;
  params.vg_N[9] = -0.07904017565835986;
  params.vg_N[10] = 1.3684127435065871;
  params.PgMin[0] = 1.4895046468487185;
  params.PgMax[0] = 1.320651812799225;
  params.RgMin[0] = -1.8279505340118756;
  params.RgMax[0] = 1.2673311275751495;
}
