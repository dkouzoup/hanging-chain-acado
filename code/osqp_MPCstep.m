function [ output ] = osqp_MPCstep(input, opts, time_)

% OSQP_MPCSTEP Solve sparse relative QP with OSQP

p  = input.prob;
qp = build_relative_qp(input);
N  = length(qp.r);
NX = size(qp.B{1},1);
NU = size(qp.B{1},2);

% TODO: CONSTRAINTS

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

% ub(ub > 1e8) = 1e8;

Ais = build_inequality_constraints_sparsity(N, NX, NU);
Ais(Ais == 0) = nan;

%% UPDATE DATA

A  = [Ae; Ais];
Ax = A(:);
Ax(isnan(Ax)) = [];

l = [be; lb];
u = [be; ub];

p.update('q', g, 'l', l, 'u', u, 'Ax', Ax);

%% SOLVE RELATIVE QP
res = p.solve();

%% STORE RESULTS

osqp_delta_x = [qp.lbx{1} reshape(res.x(N*NU+1:end), NX, N)];
osqp_delta_u = reshape(res.x(1:N*NU), NU, N);

osqp_x = input.x' + osqp_delta_x;
osqp_u = input.u' + osqp_delta_u;

output.x = osqp_x';
output.u = osqp_u';

output.info.QP_time = res.info.run_time; % TODO: this one?
output.info.nIterations = res.info.iter;
output.info.status = res.info.status;

if ~strcmp(res.info.status, 'solved')
    warning('OSQP did not solve the problem!')
    keyboard
end

end

