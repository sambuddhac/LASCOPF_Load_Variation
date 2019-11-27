% csolve  Solves a custom quadratic program very rapidly.
%
% [vars, status] = csolve(params, settings)
%
% solves the convex optimization problem
%
%   minimize(c2*square(Pg) + c1*Pg + c0 + (beta/2)*(square(Pg - PgNu) + square(PgNext - PgNextNu)) + (betaInner/2)*square(Pg - PgNuInner) + gammaSC*sum(BSC*Pg) + sum(lambda_1SC*Pg) + gamma*(B*Pg + D*PgNext) + lambda_1*Pg + lambda_2*PgNext + (rho/2)*(square(Pg - Pg_N_init + Pg_N_avg + ug_N) + square(Thetag - Vg_N_avg - Thetag_N_avg + vg_N)))
%   subject to
%     PgMin <= Pg
%     Pg <= PgMax
%     RgMin <= PgNext - Pg
%     PgNext - Pg <= RgMax
%     RgMin <= Pg - PgPrev
%     Pg - PgPrev <= RgMax
%
% with variables
%       Pg   1 x 1
%   PgNext   1 x 1
%   Thetag   1 x 1
%
% and parameters
%        B   1 x 1
%      BSC 500 x 1
%        D   1 x 1
%    PgMax   1 x 1    positive
%    PgMin   1 x 1    positive
% PgNextNu   1 x 1
%     PgNu   1 x 1
% PgNuInner   1 x 1
%   PgPrev   1 x 1    positive
% Pg_N_avg   1 x 1
% Pg_N_init   1 x 1
%    RgMax   1 x 1    positive
%    RgMin   1 x 1    negative
% Thetag_N_avg   1 x 1
% Vg_N_avg   1 x 1
%     beta   1 x 1    positive
% betaInner   1 x 1    positive
%       c0   1 x 1    positive
%       c1   1 x 1    positive
%       c2   1 x 1    positive
%    gamma   1 x 1    positive
%  gammaSC   1 x 1    positive
% lambda_1   1 x 1
% lambda_1SC 500 x 1
% lambda_2   1 x 1
%      rho   1 x 1    positive
%     ug_N   1 x 1
%     vg_N   1 x 1
%
% Note:
%   - Check status.converged, which will be 1 if optimization succeeded.
%   - You don't have to specify settings if you don't want to.
%   - To hide output, use settings.verbose = 0.
%   - To change iterations, use settings.max_iters = 20.
%   - You may wish to compare with cvxsolve to check the solver is correct.
%
% Specify params.B, ..., params.vg_N, then run
%   [vars, status] = csolve(params, settings)
% Produced by CVXGEN, 2019-02-05 18:01:13 -0500.
% CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2017 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
