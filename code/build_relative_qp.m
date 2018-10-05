function [qp] = build_relative_qp(input)

% BUILD_RELATIVE_QP derive the relative QP that ACADO solves, based on the
%   current linearization trajectory

%% extract data from input struct

NX           =  size(input.x,2);
NU           =  size(input.u,2);
N            =  size(input.u,1);
Q            =  input.W(1:NX,1:NX);
R            =  input.W(NX+1:NX+NU,NX+1:NX+NU);
QN           =  input.WN;
state_traj   =  input.x';
control_traj =  input.u';
X_ref        =  [input.y(:,1:NX);input.yN]';
U_ref        =  input.y(:,NX+1:end)';
lbu          =  input.LBU;
ubu          =  input.UBU;

M            =  (NX/3-1)/2;
WALL         =  input.WALL;
lbx          =  -inf*ones(NX,1);
ubx          =  inf*ones(NX,1);

lbx(2:3:3*(M+1)) = WALL;

%% build relative QP

SensX = cell(N,1);
SensU = cell(N,1);
SensC = cell(N,1);
GradX = cell(N+1,1);
GradU = cell(N,1);

INT_NAME = sprintf('integrate_chain_M%d_tmp', M+2);
for ii = 1:N
    sim_input.x = state_traj(:,ii);
    sim_input.u = control_traj(:,ii);
    [states,out] = eval([ INT_NAME,'(sim_input)']);
    if out.errorCode
        error('integrator failed')
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

%% save result to struct

qp.Q = cell(N+1,1);
qp.R = cell(N,1);

for ii = 1:N
    qp.Q{ii} = Q;
    qp.R{ii} = R;
end
qp.Q{N+1} = QN;

qp.q = GradX;
qp.r = GradU;

qp.A = SensX;
qp.B = SensU;
qp.c = SensC;

qp.lbx = delta_lbx;
qp.lbu = delta_lbu;
qp.ubx = delta_ubx;
qp.ubu = delta_ubu;

end

