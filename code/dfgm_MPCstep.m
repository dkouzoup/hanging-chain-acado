function [ output ] = dfgm_MPCstep(input, time_)

%DFGM_MPCSTEP SOLVE QP WITH DFGM

%% extract data from input

NX           =  size(input.x,2);
NU           =  size(input.u,2);
N            =  size(input.u,1);
M            =  (NX/3-1)/2;
Q            =  input.W(1:NX,1:NX);
R            =  input.W(NX+1:end,NX+1:end);
QN           =  input.WN;
state_traj   =  input.x';
control_traj =  input.u';
X_ref        =  [input.y(:,1:NX);input.yN]';
U_ref        =  input.y(:,NX+1:end)';
lbu          =  input.LBU;
ubu          =  input.UBU;
WALL         =  input.WALL;

lbx          =  -inf*ones(NX,1);
ubx          =  inf*ones(NX,1);
lbx(2:3:3*(M+1)) = WALL;

%% build relative QP data

SensX = cell(N,1);
SensU = cell(N,1);
SensC = cell(N,1);
GradX = cell(N+1,1);
GradU = cell(N,1);

INT_NAME = sprintf('integrate_chain_M%d_fiordos',M+2);
for ii = 1:N
    sim_input.x = state_traj(:,ii);
    sim_input.u = control_traj(:,ii);
    [states,out] = eval([ INT_NAME,'(sim_input)']);
    if out.errorCode
        error('integrator is failed')
    end
    SensX{ii,1} = states.sensX;
    SensU{ii,1} = states.sensU;
    SensC{ii,1} = states.value - state_traj(:,ii+1);
    GradX{ii,1} = Q*(sim_input.x-X_ref(:,ii));
    GradU{ii,1} = R*(sim_input.u-U_ref(:,ii));
end
GradX{N+1,1} = QN*(state_traj(:,N+1)-X_ref(:,N+1));

% calculate relative bounds
delta_lbu = cell(N+1,1);
delta_ubu = cell(N+1,1);
delta_lbx = cell(N+1,1);
delta_ubx = cell(N+1,1);

for ii = 0:N-1
    delta_lbu{ii+1} = lbu - control_traj(:,ii+1);
    delta_ubu{ii+1} = ubu - control_traj(:,ii+1);
    delta_lbx{ii+2} = lbx - state_traj(:,ii+2);
    delta_ubx{ii+2} = ubx - state_traj(:,ii+2);
end
delta_lbx{1} = input.x0' - state_traj(:,1);
delta_ubx{1} = input.x0' - state_traj(:,1);

%% SET UP OBJECTIVE

H = []; 
f = [];
for k = 1:N
    H = blkdiag(H, Q, R);
    f = [f; GradX{k}; GradU{k}];
end
H = blkdiag(H, QN);
f = [f; GradX{N+1}];

%% SETUP DYNAMICS

Am  = horzcat(SensX{:});
Bm  = horzcat(SensU{:});
beq = [delta_lbx{1}; -vertcat(SensC{:})];

%% SET UP LOWER AND UPPER BOUNDS

lb = [];
ub = [];
for k = 1:N
    lb = [lb; delta_lbx{k}; delta_lbu{k}];
    ub = [ub; delta_ubx{k}; delta_ubu{k}];
end
lb = [lb; delta_lbx{N+1};];
ub = [ub; delta_ubx{N+1};];

% lb(1:NX) = -inf*ones(NX,1);
% ub(1:NX) = +inf*ones(NX,1);
% lb(lb < -1e6) = -1e6;
% ub(ub > 1e6)  = 1e6;

%% SET UP OPTIONS

OPT.maximumIterations       = 10000;
OPT.tolerance               = 1e-3;
OPT.calculateAllMultipliers = 1; % or 1 to calculate KKT's
OPT.useExternalLibraries    = 1; % CONTROLLED AT COMPILE TIME! LEAVE TO 1, THERE IS A BUG FOR 0!
OPT.terminationCondition    = 1; % CONTROLLED AT COMPILE TIME!
OPT.warmStart               = 0; % TODO

%% SOLVE
[sol, ~, timeElapsed, it, LAMBDA, MU] = mexedDGM(H, f, Am, Bm, beq, ub, lb, OPT);

sol_xu = reshape([sol; nan(NU,1)], NX+NU, N+1);

dfgm_delta_x = sol_xu(1:NX,:);
dfgm_delta_u = sol_xu(NX+1:end,1:N);


dfgm_x = state_traj + dfgm_delta_x;
dfgm_u = control_traj + dfgm_delta_u;

%% LOG

if 0
    disp('')
    acado_x   = input.acado_sol.x';
    acado_u   = input.acado_sol.u';
    
    acado_delta_x = acado_x - state_traj;
    acado_delta_u = acado_u - control_traj;
    
    for ii = 2:N
        err = max(abs(acado_delta_x(:, ii) - SensX{ii-1}*acado_delta_x(:,ii-1) - SensU{ii-1}*acado_delta_u(:,ii-1) - SensC{ii-1}));
        disp(['acado error in dynamics' num2str(err)])
    end
    
    for ii = 2:N
        err = max(abs(dfgm_delta_x(:, ii) - SensX{ii-1}*dfgm_delta_x(:,ii-1) - SensU{ii-1}*dfgm_delta_u(:,ii-1) - SensC{ii-1}));
        disp(['dfgm error in dynamics' num2str(err)])
    end
    
    disp(['error in x between acado and dfgm ' num2str(max(max(abs(acado_x-dfgm_x))))]);
    disp(['error in u between acado and dfgm ' num2str(max(max(abs(acado_u-dfgm_u))))]);
    disp('')
    keyboard
end 

output.x   = dfgm_x';
output.u   = dfgm_u';

output.info.QP_time = timeElapsed;
output.info.nIterations = it;
output.info.status = nan;

end

