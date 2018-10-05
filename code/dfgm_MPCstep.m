function [ output ] = dfgm_MPCstep(input, time_)

% DFGM_MPCSTEP Solve relative QP with the generalized dual fast gradient
%              method used in the simulations of Kouzoupis2015 (First Order
%              Methods in Embedded Nonlinear Model Predictive Control)

WARMSTART = 0;
ITER      = 10000;
TOL       = 1e-3;

persistent LAMBDA;

qp = build_relative_qp(input);
N  = length(qp.r);
NX = size(qp.B{1},1);
NU = size(qp.B{1},2);

if time_ == 0
    LAMBDA = zeros((N+1)*NX,1);
else
    if WARMSTART == 2
        % shift
        LAMBDA = [LAMBDA(NX+1:end); LAMBDA(end-NX+1:end)];
        LAMBDA(1:NX) = -LAMBDA(1:NX); % why?
    end
end

%% SET UP OBJECTIVE

H = []; 
f = [];
for ii = 1:N
    H = blkdiag(H, qp.Q{ii}, qp.R{ii});
    f = [f; qp.q{ii}; qp.r{ii}];
end
H = blkdiag(H, qp.Q{N+1});
f = [f; qp.q{N+1}];

%% SETUP DYNAMICS

Am  = horzcat(qp.A{:});
Bm  = horzcat(qp.B{:});
beq = [qp.lbx{1}; -vertcat(qp.c{:})];

if max(max(abs(qp.lbx{1}-qp.ubx{1}))) > 1e-10
    warning('bounds on first stage do not include x0 constraint')
    keyboard
end

%% SET UP LOWER AND UPPER BOUNDS

lb = [];
ub = [];
for ii = 1:N
    lb = [lb; qp.lbx{ii}; qp.lbu{ii}];
    ub = [ub; qp.ubx{ii}; qp.ubu{ii}];
end
lb = [lb; qp.lbx{N+1};];
ub = [ub; qp.ubx{N+1};];

% lb(1:NX) = -inf*ones(NX,1);
% ub(1:NX) = +inf*ones(NX,1);

% lb(lb < -1e6) = -1e6;
% ub(ub > 1e6)  = 1e6;

%% SET UP OPTIONS

OPT.maximumIterations       = ITER;
OPT.tolerance               = TOL;
OPT.calculateAllMultipliers = 0; % or 1 to calculate KKT's
OPT.useExternalLibraries    = 1; % CONTROLLED AT COMPILE TIME! LEAVE TO 1, THERE IS A BUG FOR 0!
OPT.terminationCondition    = 1; % CONTROLLED AT COMPILE TIME!
OPT.warmStart               = double(WARMSTART > 0);

%% SOLVE RELATIVE QP

[sol, ~, timeElapsed, it, LAMBDA, ~] = mexedDGM(H, f, Am, Bm, beq, ub, lb, OPT, LAMBDA);

sol_xu = reshape([sol; nan(NU,1)], NX+NU, N+1);

dfgm_delta_x = sol_xu(1:NX,:);
dfgm_delta_u = sol_xu(NX+1:end,1:N);

dfgm_x = input.x' + dfgm_delta_x;
dfgm_u = input.u' + dfgm_delta_u;

%% SAVE RESULTS

output.x   = dfgm_x';
output.u   = dfgm_u';

output.info.QP_time = timeElapsed;
output.info.nIterations = it;
output.info.status = nan;

end

