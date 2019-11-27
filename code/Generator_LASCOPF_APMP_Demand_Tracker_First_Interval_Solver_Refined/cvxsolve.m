% Produced by CVXGEN, 2019-02-05 18:01:13 -0500.
% CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2017 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: cvxsolve.m.
% Description: Solution file, via cvx, for use with sample.m.
function [vars, status] = cvxsolve(params, settings)
B = params.B;
BSC = params.BSC;
D = params.D;
PgMax = params.PgMax;
PgMin = params.PgMin;
PgNextNu = params.PgNextNu;
PgNu = params.PgNu;
PgNuInner = params.PgNuInner;
PgPrev = params.PgPrev;
Pg_N_avg = params.Pg_N_avg;
Pg_N_init = params.Pg_N_init;
RgMax = params.RgMax;
RgMin = params.RgMin;
Thetag_N_avg = params.Thetag_N_avg;
Vg_N_avg = params.Vg_N_avg;
beta = params.beta;
betaInner = params.betaInner;
c0 = params.c0;
c1 = params.c1;
c2 = params.c2;
gamma = params.gamma;
gammaSC = params.gammaSC;
lambda_1 = params.lambda_1;
lambda_1SC = params.lambda_1SC;
lambda_2 = params.lambda_2;
rho = params.rho;
ug_N = params.ug_N;
vg_N = params.vg_N;
cvx_begin
  % Caution: automatically generated by cvxgen. May be incorrect.
  variable Pg;
  variable PgNext;
  variable Thetag;

  minimize(c2*square(Pg) + c1*Pg + c0 + (beta/2)*(square(Pg - PgNu) + square(PgNext - PgNextNu)) + (betaInner/2)*square(Pg - PgNuInner) + gammaSC*sum(BSC*Pg) + sum(lambda_1SC*Pg) + gamma*(B*Pg + D*PgNext) + lambda_1*Pg + lambda_2*PgNext + (rho/2)*(square(Pg - Pg_N_init + Pg_N_avg + ug_N) + square(Thetag - Vg_N_avg - Thetag_N_avg + vg_N)));
  subject to
    PgMin <= Pg;
    Pg <= PgMax;
    RgMin <= PgNext - Pg;
    PgNext - Pg <= RgMax;
    RgMin <= Pg - PgPrev;
    Pg - PgPrev <= RgMax;
cvx_end
vars.Pg = Pg;
vars.PgNext = PgNext;
vars.Thetag = Thetag;
status.cvx_status = cvx_status;
% Provide a drop-in replacement for csolve.
status.optval = cvx_optval;
status.converged = strcmp(cvx_status, 'Solved');
