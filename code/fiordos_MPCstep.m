function output = fiordos_MPCstep( input, time_)

persistent sol_x sol_la;

if time_ == 0
    sol_x = []; 
    sol_la = [];
end

% CAUTION!! WEIGHTS HARDCODED DURING CODE GENERATION

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

%% build fiordos data

Aeu  = [];
Aex1 = [];
Aex2 = [];
for ii = 1:N
    Aeu  = blkdiag(Aeu, -SensU{ii});
    Aex1 = blkdiag(Aex1, eye(NX));
    if ii < N
        Aex2 = blkdiag(Aex2, -SensX{ii+1});
    end
end

Aex = Aex1 + [zeros(NX, size(Aex1,2)); Aex2 zeros(size(Aex2,1) ,NX)];

be = vertcat(SensC{:});
be(1:NX) = be(1:NX) + SensX{1}*delta_lbx{1};

g = [vertcat(GradU{:}); vertcat(GradX{2:end})];

%% solve QP with fiordos

mparams = struct();
msetgs  = struct();

mparams.Ae = [Aeu Aex];
mparams.be = be;
mparams.g  = g;
for ii = 1:N
    % NOTE: fiordos gives wrong results with matlab inf!
    % TODO: can I skip constraints with -+ inf bounds?
    delta_lbx{ii}(isinf(delta_lbx{ii})) = -1e8;
    delta_ubx{ii}(isinf(delta_ubx{ii})) = +1e8;
    eval(['mparams.X' num2str(ii) '.l = delta_lbu{' num2str(ii) '};']);
    eval(['mparams.X' num2str(ii) '.u = delta_ubu{' num2str(ii) '};']);
    eval(['mparams.X' num2str(N+ii) '.l = delta_lbx{' num2str(ii+1) '};']);
    eval(['mparams.X' num2str(N+ii) '.u = delta_ubx{' num2str(ii+1) '};']);
end

if time_ > 0 && input.warmstart
    msetgs.approach.apprInitX  = sol_x;                  
    msetgs.approach.apprInitLa = sol_la;
end

tic;
mres = fiordos_mpc_mex(mparams, msetgs);
tmp_cputime = toc;

disp(['fiordos returned status ' num2str(mres.exitflag) ' in ' num2str(mres.iter) ' iterations'])

fiordos_delta_x = [delta_lbx{1} reshape(mres.x(N*NU+1:end), NX, N)];
fiordos_delta_u = reshape(mres.x(1:N*NU), NU, N);

fiordos_x = state_traj + fiordos_delta_x;
fiordos_u = control_traj + fiordos_delta_u;

if input.warmstart
    sol_x  = mres.x;
    sol_la = mres.la;
end

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
        err = max(abs(fiordos_delta_x(:, ii) - SensX{ii-1}*fiordos_delta_x(:,ii-1) - SensU{ii-1}*fiordos_delta_u(:,ii-1) - SensC{ii-1}));
        disp(['fiordos error in dynamics' num2str(err)])
    end
    
    disp(['error in x between acado and fiordos ' num2str(max(max(abs(acado_x-fiordos_x))))]);
    disp(['error in u between acado and fiordos ' num2str(max(max(abs(acado_u-fiordos_u))))]);
    disp('')
    keyboard
end 

output.x   = fiordos_x';
output.u   = fiordos_u';

if isfield(mres, 'cputime')
    output.info.QP_time = mres.cputime;
else
    output.info.QP_time = nan;    
end

output.info.nIterations = mres.iter;
output.info.status = mres.exitflag;

end
