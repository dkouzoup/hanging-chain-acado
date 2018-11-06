function [primal_res, dual_res, slack_res] = check_solution_accuracy_acado_hpmpc(input, output)

% CHECK_SOLUTION_ACCURACY_ACADO_HPMPC When using acado-dev and ONLY for the
%                                     chain of masses example, calculate
%                                     primal and dual residuals of solution

%% build relative QP that ACADO solves

qp = build_relative_qp(input);
N  = length(qp.r);
NX = size(qp.B{1},1);
NU = size(qp.B{1},2);
M  = (NX/3 - 1)/2;

%% write QP in sparse form

H = [];
f = [];
for ii = 1:N
    H = blkdiag(H, qp.Q{ii}, qp.R{ii});
    f = [f; qp.q{ii}; qp.r{ii}];
end
H = blkdiag(H, qp.Q{N+1});
f = [f; qp.q{N+1}];

Aeq = zeros(N*NX+NX,N*(NX+NU)+NX);

Aeq(1:NX,1:NX) = eye(NX);
for k = 1:N
    Aeq(NX+(k-1)*NX+1:NX+k*NX,(k-1)*(NX+NU)+1:k*(NX+NU)+NX) = [qp.A{k} qp.B{k} -eye(NX)];
end
Aeq(1:NX,1:NX) = eye(NX);

beq = [qp.lbx{1}; -vertcat(qp.c{:})];

Aineq = eye(size(H,1)); 
lb    = [];
ub    = [];

for ii = 1:N
    lb = [lb; qp.lbx{ii}; qp.lbu{ii}];
    ub = [ub; qp.ubx{ii}; qp.ubu{ii}];
end
lb = [lb; qp.lbx{N+1}];
ub = [ub; qp.ubx{N+1}];

%% get (relative) solution in sparse form

delta_x = output.x - input.x;
delta_u = [output.u - input.u; nan(1, NU)];

delta_z = [delta_x delta_u]';
delta_z = delta_z(1:end-NU)';

%% calculate primal residual

primal_res_eq   = max(abs(Aeq*delta_z - beq));
primal_res_ineq = max([lb - delta_z; delta_z - ub]);
primal_res      = max(primal_res_eq, primal_res_ineq);

%% sort multipliers

hpmpc_mu = reshape(output.mu, 2*(NX+NU), N+1);
hpmpc_mu(:,1) = [];

NBU  = NU;
NBX  = M+1; % except for k = 0 where NBX = 0
IDXB = 2:3:NX/2+1;

hpmpc_lbu = [output.mu(1:NBU) hpmpc_mu(1:NBU, :)];
hpmpc_ubu = [output.mu(NBU+1:2*NBU) hpmpc_mu(NBU+NBX+1:2*NBU+NBX, :)];
hpmpc_lbx = zeros(NX, N+1); hpmpc_lbx(IDXB,2:end) = hpmpc_mu(NBU+1:NBX+NBU,:);
hpmpc_ubx = zeros(NX, N+1);

MU_1 = [hpmpc_ubx; hpmpc_ubu];
MU_1 = MU_1(NX+1:end-NU)';
MU_2 = [hpmpc_lbx; hpmpc_lbu];
MU_2 = MU_2(NX+1:end-NU)';

%% calculate dual residual

% x0 variable is eliminated in HPMPC
Aeq = Aeq(NX+1:end,NX+1:end);
beq = beq(NX+1:end);
H   = H(NX+1:end, NX+1:end);
f   = f(NX+1:end);
dz  = delta_z(NX+1:end);

tmp = H*dz + f + Aeq'*output.lam + MU_1 - MU_2;

dual_res = max(abs(tmp));

%% calculate residual for complementary slackness

% TODO
slack_res = nan;

end

