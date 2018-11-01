function output = osqp_MPCstep(input, opts, time_)

% OSQP_MPCSTEP Solve sparse relative QP with OSQP solver. Ordering of
%              variables is: [u_{0}, ..., u_{N-1}, x_{1}, ..., x{N}]

%% BUILD RELATIVE QP

p  = input.prob;
qp = build_relative_qp(input);
N  = length(qp.r);
NX = size(qp.B{1},1);
NU = size(qp.B{1},2);

%% SET UP DYNAMICS 

if max(max(abs(qp.lbx{1}-qp.ubx{1}))) > 1e-10
    warning('bounds on first stage do not include x0 constraint')
    keyboard
end

[Ae, be] = build_equality_constraints(qp.A, qp.B, qp.c, qp.lbx{1});

% convert Ae to vector of non-zero elements (zeros inside dynamics don't count)
Aes = build_equality_constraints_sparsity(N, NX, NU);
Aes(Aes == 0) = nan;

Ae = Ae.*Aes;

%% SET UP OBJECTIVE

g = [vertcat(qp.r{:}); vertcat(qp.q{2:end})];

%% SET UP BOUNDS

lb = [vertcat(qp.lbu{:}); vertcat(qp.lbx{2:end})];
ub = [vertcat(qp.ubu{:}); vertcat(qp.ubx{2:end})];

ub(isinf(lb)) = [];
lb(isinf(lb)) = [];

Ais = build_inequality_constraints_sparsity(N, NX, NU);
Ais(Ais == 0) = nan;

%% UPDATE OSQP DATA

A  = [Ae; Ais];
Ax = A(:);
Ax(isnan(Ax)) = [];

l = [be; lb];
u = [be; ub];

p.update('q', g, 'l', l, 'u', u, 'Ax', Ax);

%% SOLVE RELATIVE QP

res = p.solve();

%% STORE ABSOLUTE SOLUTION AND INFO

osqp_delta_x = [qp.lbx{1} reshape(res.x(N*NU+1:end), NX, N)];
osqp_delta_u = reshape(res.x(1:N*NU), NU, N);

osqp_x = input.x' + osqp_delta_x;
osqp_u = input.u' + osqp_delta_u;

output.x   = osqp_x';
output.u   = osqp_u';

% TODO: return matrices instead of vectors, re-order mu
% output.lam = res.y(1:size(Ae,1));
% output.mu  = res.y(size(Ae,1)+1:end);

% TODO: update time is needed here but not returned by solver (yet)
if opts.with_setup_time
    output.info.cpuTime = res.info.setup_time + res.info.solve_time;
else
    output.info.cpuTime = res.info.solve_time;
end

output.info.objValue = res.info.obj_val;
output.info.nIterations = res.info.iter;

if ~strcmp(res.info.status, 'solved')
    warning('OSQP did not solve the problem!')
    output.info.status = -1;
    keyboard
else
    output.info.status = 0;
end

%% CALCULATE RESIDUALS

% primal res
A = sparse([build_equality_constraints_sparsity(N, NX, NU);
    build_inequality_constraints_sparsity(N, NX, NU)]);

A(A~=0) = Ax;

pri_res = max([A*res.x - u; l - A*res.x]); % typically <= res.info.pri_res
output.info.primal_res = pri_res;

% dual res
Q  = input.W(1:NX,1:NX);
R  = input.W(NX+1:NX+NU,NX+1:NX+NU);
QN = input.WN;
H  = sparse(blkdiag(kron(eye(N), R),  kron(eye(N-1), Q), QN));

dua_res = max(abs(H*res.x + g + A'*res.y));

% check consistency
res_tol = 1e-10;

if abs(res.info.dua_res - dua_res) > res_tol
    [res.info.dua_res, dua_res]
    warning('inconsistent residuals')
    keyboard
end

output.info.primal_res = pri_res; 
output.info.dual_res   = dua_res; 

end

