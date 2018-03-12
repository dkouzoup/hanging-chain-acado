function output = acados_MPCstep( input )

import acados.*
persistent qp;

persistent indx;
    
if isempty(indx)
    indx = 0;
else
    indx = indx+1;
end

% persistent logs;
% logs_length = 5;
% if isempty(logs)
%     logs = cell(logs_length,1);
% end
% logs{indx+1} = input;
% if indx+1 == logs_length
%     save('logs', 'logs')
%     keyboard
% end


if isempty(qp)
    first_iter = 1;
else
    first_iter = 0;
end

solver = input.solver;

SKIP_ACADOS_SOLVE = 0;
CODE_GENERATE_QPS = 0;

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

INT_NAME = sprintf('integrate_chain_M%d_acados',M+2);
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

%% code generate OCP_QP

if CODE_GENERATE_QPS
    
    % info to build filename
    qpdata.solver  = input.acado_solver;
    qpdata.nmasses = input.nmasses; 
    qpdata.warmstart = input.warmstart;
    
    % optimal solution
    qpdata.acado_iter = input.acado_sol.info.nIterations;
    qpdata.acado_sol  = [input.acado_sol.u' - control_traj; input.acado_sol.x(1:end-1,:)' - state_traj(:,1:end-1)];
    qpdata.acado_sol  = [qpdata.acado_sol(:); input.acado_sol.x(end,:)' - state_traj(:,end)];
        
    % dims
    qpdata.N  = N;
    qpdata.nx = NX*ones(1, N+1);
    qpdata.nu = [NU*ones(1, N) 0];
    
    % objective
    
    %     TODO: pass objective also as cell arays
    
    %     qpdata.Q = cell(N+1, 1);
    %     qpdata.R = cell(N+1, 1);
    %     qpdata.S = cell(N+1, 1);
    %     qpdata.q = cell(N+1, 1);
    %     qpdata.r = cell(N+1, 1);
    %
    %     for ii = 1:N
    %         qpdata.Q{ii} = Q;
    %         qpdata.R{ii} = R;
    %         qpdata.S{ii} = zeros(qpdata.nu(ii), qpdata.nx(ii));  % or transposed?
    %         qpdata.q{ii} = q;
    %         qpdata.r{ii} = r;
    %     end
    %     qpdata.Q{N+1} = QN;
    
    qpdata.Qv = [repmat(Q(:), N, 1); QN(:)];
    qpdata.Rv = [repmat(R(:), N, 1)];
    qpdata.Sv = repmat(zeros(NX*NU,1), N, 1);
    qpdata.q  = vertcat(GradX{:});
    qpdata.r  = vertcat(GradU{:});
    
    % dynamics
    
    qpdata.A = SensX;
    qpdata.B = SensU;
    qpdata.b = SensC;
    
    % constraints
    
    qpdata.lbu = delta_lbu;
    qpdata.ubu = delta_ubu;
    qpdata.lbx = delta_lbx;
    qpdata.ubx = delta_ubx;
    
    code_generate_ocp_qp(qpdata, indx, '../../acados/examples/c/ocp_qp_bugs/');
    
end

%% solve with acados

if ~SKIP_ACADOS_SOLVE
    
    if first_iter
        qp = ocp_qp(N, NX, NU);
        
        qp.set('Q', Q);
        qp.set('R', R);
        qp.set('Q', N, QN);
        
        % max iter
        opts = struct('max_iter', 1000);
        
        % warmstart (for qpOASES, QORE and qpDUNES)
        if ~contains(solver, 'hpipm') && ~contains(solver, 'hpmpc')
            opts.warm_start = input.warmstart;
            if strcmp(solver, 'qore')
                %opts.warm_strategy = 0;  % 0 or 1
            end
        end
        
        % partial condensing
        if ~strcmp(solver, 'qpoases') && ~strcmp(solver, 'qore')&& ~strcmp(solver, 'condensing_hpipm')
            if ~isinf(input.N2) && input.N2 ~= N
                opts.N2 = input.N2;
                if strcmp(solver, 'qpdunes')
                    opts.clipping = 0;
                end
            else
                opts.N2 = N;
                if strcmp(solver, 'qpdunes')
                    opts.clipping = 1;
                end
            end
        end
        
        % TODO: fix segfaults in HPMPC!!!
        if strcmp(solver, 'hpmpc')
            opts.mu0 = 1e6;
            % TODO: Why setting this (same as default arg) breaks solver?
            % opts.tol = 1e-12;
        end
        opts
        qp.initialize_solver(solver, opts);
    end
    
    % set dynamics
    for ii = 0:N-1
        qp.set('A', ii, SensX{ii+1});
        qp.set('B', ii, SensU{ii+1});
        qp.set('b', ii, SensC{ii+1});
        qp.set('q', ii, GradX{ii+1});
        qp.set('r', ii, GradU{ii+1});
    end
    qp.set('q', N, GradX{N+1});
    
    for ii = 0:N-1
        qp.set('lbu', ii, delta_lbu{ii+1});
        qp.set('ubu', ii, delta_ubu{ii+1});
        qp.set('lbx', ii+1, delta_lbx{ii+2});
        qp.set('ubx', ii+1, delta_ubx{ii+2});
    end
    
    
    qp.set('lbx', 0, delta_lbx{1});
    qp.set('ubx', 0, delta_ubx{1});
    
    if first_iter
        qp.initialize_solver(solver, opts);
    else
        % keyboard
    end
    
    % solve QP
    acados_output = qp.solve();
    acados_output.info.num_iter
    
    Delta_xopt = acados_output.states();
    Delta_uopt = acados_output.controls();

    output.x   = (state_traj   + [Delta_xopt{:}])';
    output.u   = (control_traj + [Delta_uopt{:}])';
    output.info.QP_time = acados_output.info.total_time;
    output.info.nIterations = acados_output.info.num_iter;

else
    output.x = state_traj*NaN;
    output.u = control_traj*NaN;
    output.info.QP_time = NaN;
    output.info.nIterations = NaN;
end

end
