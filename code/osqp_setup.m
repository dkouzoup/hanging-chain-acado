function prob = osqp_setup(N, NX, NU, W, WN, opts)

% OSQP_SETUP Set up OSQP solver

% set up objective
Q  = W(1:NX,1:NX);
R  = W(NX+1:NX+NU,NX+1:NX+NU);
QN = WN;
H  = sparse(blkdiag(kron(eye(N), R),  kron(eye(N-1), Q), QN));
g  = ones(size(H,1), 1);

% set up bounds
Ai = sparse(build_inequality_constraints_sparsity(N, NX, NU));

% set up equalities
Ae = sparse(build_equality_constraints_sparsity(N, NX, NU));

% merge constraints
A  = [Ae; Ai];
l  = -1*ones(size(A,1),1);
u  = +1*ones(size(A,1),1);

prob = osqp;
prob.setup(H, g, A, l, u, 'warm_start', logical(opts.warmstart), ...
    'check_termination', opts.check_ter, 'max_iter', opts.maxit, ...
    'verbose', 0, 'eps_abs', opts.abstol, 'eps_rel', opts.reltol);

end

