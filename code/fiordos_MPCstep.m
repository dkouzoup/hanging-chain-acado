function output = fiordos_MPCstep(input, opts, time_)

persistent sol_x sol_la;

if time_ == 0
    sol_x  = []; 
    sol_la = [];
end

qp = build_relative_qp(input);
N  = length(qp.r);
NX = size(qp.B{1},1);
NU = size(qp.B{1},2);

%% SET UP DYNAMICS

Aeu  = [];
Aex1 = [];
Aex2 = [];
for ii = 1:N
    Aeu  = blkdiag(Aeu, -qp.B{ii});
    Aex1 = blkdiag(Aex1, eye(NX));
    if ii < N
        Aex2 = blkdiag(Aex2, -qp.A{ii+1});
    end
end

Aex = Aex1 + [zeros(NX, size(Aex1,2)); Aex2 zeros(size(Aex2,1) ,NX)];

be = vertcat(qp.c{:});

% x0 embedding
be(1:NX) = be(1:NX) + qp.A{1}*qp.lbx{1};

if max(max(abs(qp.lbx{1}-qp.ubx{1}))) > 1e-10
    warning('bounds on first stage do not include x0 constraint')
    keyboard
end

%% SET UP OBJECTIVE

g = [vertcat(qp.r{:}); vertcat(qp.q{2:end})];

%% SET UP BOUNDS

for ii = 1:N
    % NOTE: fiordos gives wrong results with matlab inf!
    % TODO: can we skip constraints with -+ inf bounds?
    qp.lbu{ii}(isinf(qp.lbu{ii})) = -opts.infval;
    qp.ubu{ii}(isinf(qp.ubu{ii})) = +opts.infval;
    qp.lbx{ii}(isinf(qp.lbx{ii})) = -opts.infval;
    qp.ubx{ii}(isinf(qp.ubx{ii})) = +opts.infval;        
end
qp.lbx{N+1}(isinf(qp.lbx{N+1})) = -opts.infval;
qp.ubx{N+1}(isinf(qp.ubx{N+1})) = +opts.infval;

%% SOLVE RELATIVE QP

mparams = struct();
msetgs  = struct();

mparams.Ae = [Aeu Aex];
mparams.be = be;
mparams.g  = g;

for ii = 1:N
    eval(['mparams.X' num2str(ii) '.l = qp.lbu{' num2str(ii) '};']);
    eval(['mparams.X' num2str(ii) '.u = qp.ubu{' num2str(ii) '};']);
    eval(['mparams.X' num2str(N+ii) '.l = qp.lbx{' num2str(ii+1) '};']);
    eval(['mparams.X' num2str(N+ii) '.u = qp.ubx{' num2str(ii+1) '};']);
end

% NOTE: Hessian is hard-coded to allow clipping
% mparams.H = diag(blkdiag(kron(eye(N), R),  kron(eye(N-1), Q), QN));

if time_ > 0 && opts.warmstart
    if strcmp(opts.approach, 'primal-dual')
        msetgs.approach.apprInitX  = sol_x;
        msetgs.approach.apprInitLa = sol_la;
    else
        msetgs.algoInner.init = sol_x;
        msetgs.algoOuter.init = sol_la;   
    end
end

tic;
mres = fiordos_mpc_mex(mparams, msetgs);
tmp_cputime = toc;

disp(['fiordos returned status ' num2str(mres.exitflag) ' in ' num2str(mres.iter) ' iterations'])

fiordos_delta_x = [qp.lbx{1} reshape(mres.x(N*NU+1:end), NX, N)];
fiordos_delta_u = reshape(mres.x(1:N*NU), NU, N);

fiordos_x = input.x' + fiordos_delta_x;
fiordos_u = input.u' + fiordos_delta_u;

if opts.warmstart
    sol_x  = mres.x;
    sol_la = mres.la;
end

if 0
    disp('')
    acado_x   = input.acado_sol.x';
    acado_u   = input.acado_sol.u';
    
    acado_delta_x = acado_x - input.x';
    acado_delta_u = acado_u - input.u';
    
    for ii = 2:N
        err = max(abs(acado_delta_x(:, ii) - qp.A{ii-1}*acado_delta_x(:,ii-1) - qp.B{ii-1}*acado_delta_u(:,ii-1) - qp.c{ii-1}));
        disp(['acado error in dynamics' num2str(err)])
    end
    
    for ii = 2:N
        err = max(abs(fiordos_delta_x(:, ii) - qp.A{ii-1}*fiordos_delta_x(:,ii-1) - qp.B{ii-1}*fiordos_delta_u(:,ii-1) - qp.c{ii-1}));
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
    if 0
        fprintf('\n\nTIMING IN C CODE %3.2f TIMES FASTER\n\n', tmp_cputime/mres.cputime);
    end
else
    output.info.QP_time = nan;    
end

output.info.nIterations = mres.iter;
output.info.status = mres.exitflag;

end
