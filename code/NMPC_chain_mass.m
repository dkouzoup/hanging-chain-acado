function logged_data = NMPC_chain_mass(sim_opts)

%% INITIALIZE

if ~exist('sim_opts','var')
    clear all
    BATCH = false;
else
    BATCH = true;
end

clear global % for ACADO workspace
clearvars -except sim_opts BATCH; clc; close all; % sim_opts contain options to overwrite current script
close all
clc

addpath([pwd filesep 'utils'])

%% SIMULATION OPTIONS

SIM_EXPORT  = 1;            % export code for ACADO simulator

SIM_COMPILE = 1;            % compile exported code for ACADO simulator

MPC_EXPORT  = 1;            % export code for ACADO solver

MPC_COMPILE = 1;            % compile exported code for ACADO solver

NRUNS       = 5;            % run closed-loop simulation NRUNS times and store minimum timings (to minimize OS interference)

SOLVER      = 'HPMPC_B10'; % 'qpDUNES_BXX' (with XX block size, 0 for clipping), 'qpOASES_N2', 'qpOASES_e_N3', 'qpOASES_e_N2', 'qpOASES_N3', 'FORCES', 'HPMPC', 'dfgm', 'osqp', 'fiordos'

WARMSTART   = 0;            % applicable for qpOASES/qpDUNES

VISUAL      = 0;            % set to 1 to visualize chain of masses (only for the first out of the NRUNS simulations)

WALL        = -0.1;         % wall position (re-export if changed)

Ts          = 0.2;          % sampling time [s]

Tf          = 5;            % final simulation time [s]

N           = 50;           % prediction horizon

NMASS       = 4;            % number of masses in chain (data available from 3 to 6 masses)

To          = 5*Ts;         % how many seconds to overwrite uMPC

Uo          = [-1;1;1];     % value to overwrite uMPC for To sec

CHECK_AGAINST_REF_SOL = 1;  % if 1, exports and compiles reference solver (qpOASES, CN2 by default)

SOL_TOL = 1e-6;             % maximum accepted 2-norm of the deviation of the solution from the reference solution
                            % (only used if CHECK_AGAINST_REF_SOL = 1).

ACADOS_BENCHMARK = true;    % change some settings to match what is currently supported in acados (ERK, no shifting, RTI timings instead of QP timings)


%% Load simulation options and overwrite local ones

if exist('sim_opts','var')
    names  = fieldnames(sim_opts);
    for ii = 1:length(names)
        eval([names{ii} ' = sim_opts.' names{ii} ';'])
    end
end

%% Initialization

if ACADOS_BENCHMARK
    USE_EXPLICIT_RK = true; % currently only ERK supported in acados-matlab
    NUM_STEPS       = 2;    % NEEDS TO BE HARDCODED IN ACADOS DEFAULT OPTS ATM (one step - the default option - becomes infeasible in closed loop)
else
    USE_EXPLICIT_RK = false;
    NUM_STEPS       = 2;
end

% dimensions
M  = NMASS - 2;     % number of intermediate masses
NX = (2*M + 1)*3;   % differential states
NU = 3;             % control inputs

% MPC weights (hard-coded for some solvers)
W  = blkdiag(25*eye(3), 25*eye(3*M), 1*eye(3*M), 0.01*eye(NU));
WN = blkdiag(25*eye(3), 25*eye(3*M), 1*eye(3*M));
% W  = blkdiag(15*eye(3), 10*eye(length(x)), eye(length(v)), 0.05*eye(3));
% WN = blkdiag(15*eye(3), 10*eye(length(x)), eye(length(v)));

% solver-specific options
% TODO move ACADO opts also here
switch SOLVER

    case 'fiordos'

        opts.approach  = 'dual'; % 'dual' or 'primal-dual'
        opts.tol       = 1e-2;
        opts.maxit     = 100000;
        opts.export    = MPC_EXPORT;
        opts.compile   = MPC_COMPILE;
        opts.infval    = 1e8;
        opts.warmstart = 1; % 0: no warmstart, 1: same solution

    case 'dfgm'

        opts.warmstart = 2; % 0: no warmstart, 1: same solution, 2: shifted solution
        opts.maxit     = 100000;
        opts.tol       = 1e-2;
        opts.criterion = 1; % 1: solver's own criterion, 2: compare with acado solution (not working yet properly)
        opts.useExternalLibraries = 0;

    case 'osqp'
        % leave empty for default
        opts.warmstart       = 1;
        opts.maxit           = 1000;
        opts.check_ter       = [];
        opts.tol             = [];
        opts.with_setup_time = true;
end


% extract block size from solver name
QPCONDENSINGSTEPS = extract_block_size(SOLVER, N);

DifferentialState xEnd(3,1);                           % 3-dimensional position of end point (M+1)
eval(['DifferentialState x(',num2str(M*3),',1);'])     % 3-dimensional position of masses 1, 2, ..., M
eval(['DifferentialState v(',num2str(M*3),',1);'])     % 3-dimensional velocity of masses 1, 2, ..., M
Control u(3,1);

x0 = [0; 0; 0];
L  = 0.033;
D  = 1.0;
m  = 0.03;

%% Differential Equation

A1 = zeros(3,3);
B1 = eye(3,3);

ode_rhs = acado.Expression(zeros(3*M, 1));
ode_rhs = chain_dynamics(x, v, xEnd, L, D, m, M, x0, ode_rhs);

ode = [ dot(x) == ode_rhs(1:3*M); ...
        dot(v) == ode_rhs(3*M+1:2*3*M)];

% compute rest position with fsolve

xEnd_ref = [1; 0; 0]; % position of actuated mass

[fsolve_x, ~, fsolve_flag] = fsolve(@(x)chain_dynamics(x, zeros(3*M, 1), xEnd_ref, L, D, m, M, x0, zeros(3*M,1)), (1:3*M).', optimoptions('fsolve','Algorithm','Levenberg-Marquardt', 'Display', 'off'));

if fsolve_flag ~= 1
    error(['Calculation of rest position failed. fsolve returned status ' num2str(fsolve_flag) '. Check xEnd_ref']);
end

fsolve_ref = [xEnd_ref; fsolve_x; zeros(3*M, 1)];

% build ode as casadi expression for acados
if contains(SOLVER, 'acados')
    import casadi.*
    try
        x_casadi   = SX.sym('x', 3*M);
        v_casadi   = SX.sym('v', 3*M);
        xE_casadi  = SX.sym('xEnd', 3);
        u_casadi   = SX.sym('u', 3);
        ode_casadi = SX.zeros(3*M, 1);
        ode_casadi = chain_dynamics(x_casadi, v_casadi, xE_casadi, L, D, m, M, x0, ode_casadi);
        ode_casadi = [u_casadi; ode_casadi];
        ode_ca_fun = Function('ode_fun', {[xE_casadi; x_casadi; v_casadi], u_casadi}, {ode_casadi});
    catch
        error('casadi-matlab is probably not in your path');
    end

     acados_nlp = acados_setup(NMASS, ode_ca_fun, WALL, N, Ts, W, WN, fsolve_ref, SOLVER, WARMSTART);
end

%% SIMexport

acadoSet('problemname', 'sim');

sim = acado.SIMexport(Ts);

if USE_EXPLICIT_RK
    sim.setModel([u; ode_rhs(1:3*M); ode_rhs(3*M+1:2*3*M)]);
    sim.set( 'INTEGRATOR_TYPE',            'INT_RK4' );
else
    sim.setLinearInput(A1, B1);
    sim.setModel(ode);
    sim.set( 'INTEGRATOR_TYPE',        'INT_IRK_GL2' );
end
sim.set( 'NUM_INTEGRATOR_STEPS',        NUM_STEPS    );

if SIM_EXPORT
    sim.exportCode( 'export_SIM' );
end

sim_name = sprintf('integrate_chain_M%d',NMASS);
if SIM_COMPILE
    cd export_SIM
    sim_path = '../';
    make_acado_integrator([sim_path sim_name])
    if ~is_acado_solver(SOLVER) || contains(SOLVER, 'HPMPC')
        % NOTE: using the same integrator for simulation gives somehow wrong results
        make_acado_integrator([sim_path sim_name '_tmp'])
    end
    cd ..
end

%% MPCexport

acadoSet('problemname', 'mpc');

ocp = acado.OCP( 0.0, N*Ts, N );

rf = [xEnd; x; v; u];
S  = acado.BMatrix(eye(NX+NU));

ocp.minimizeLSQ( S, rf );

rfN = [xEnd; x; v];
SN  = acado.BMatrix(eye(NX));

ocp.minimizeLSQEndTerm( SN, rfN );

ocp.subjectTo( WALL <= [x([2:3:end]); xEnd(2)]); % state constraint on positions
ocp.subjectTo( -1 <= u <= 1 ); % bounds on controls

if USE_EXPLICIT_RK
    ocp.setModel([u; ode_rhs(1:3*M); ode_rhs(3*M+1:2*3*M)]);
else
    ocp.setLinearInput(A1,B1);
    ocp.setModel(ode);
end

mpc = acado.OCPexport( ocp );

mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'       );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING'  );

if USE_EXPLICIT_RK
    mpc.set( 'INTEGRATOR_TYPE',             'INT_RK4'        );
else
    mpc.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL2'    );
end

mpc.set( 'NUM_INTEGRATOR_STEPS',        NUM_STEPS*N          );

if strcmp(SOLVER,'qpOASES_N3')

    mpc.set( 'QP_SOLVER',               'QP_QPOASES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'CONDENSING'         );
    if WARMSTART
        mpc.set( 'HOTSTART_QP',         'YES'                );
    end

elseif strcmp(SOLVER,'qpOASES_N2')

    mpc.set( 'QP_SOLVER',               'QP_QPOASES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'FULL_CONDENSING_N2' );
    if WARMSTART
        mpc.set( 'HOTSTART_QP',         'YES'                );
    end

elseif strcmp(SOLVER,'qpOASES_e_N3')

    mpc.set( 'QP_SOLVER',               'QP_QPOASES3'        );
    mpc.set( 'SPARSE_QP_SOLUTION',      'CONDENSING'         );
    if WARMSTART
        mpc.set( 'HOTSTART_QP',         'YES'                );
    end

elseif strcmp(SOLVER,'qpOASES_e_N2')

    mpc.set( 'QP_SOLVER',               'QP_QPOASES3'        );
    mpc.set( 'SPARSE_QP_SOLUTION',      'FULL_CONDENSING_N2' );
    if WARMSTART
        mpc.set( 'HOTSTART_QP',         'YES'                );
    end

elseif contains(SOLVER,'qpDUNES') && QPCONDENSINGSTEPS <= 1

    mpc.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'SPARSE_SOLVER'      );
    mpc.set( 'HOTSTART_QP',             'YES'                ); % probably not needed
    if WARMSTART == 0
        warning('qpDUNES with cold start not implemented in ACADO');
        keyboard
    end

elseif contains(SOLVER,'qpDUNES') && QPCONDENSINGSTEPS > 1

    mpc.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'BLOCK_CONDENSING_N2');
    mpc.set( 'CONDENSING_BLOCK_SIZE',    QPCONDENSINGSTEPS   );
    mpc.set( 'HOTSTART_QP',             'YES'                ); % probably not needed

    if WARMSTART == 0
        warning('qpDUNES with cold start not implemented in ACADO');
        keyboard
    end

elseif strcmp(SOLVER,'FORCES_BC')

    error('interface of FORCES with block condensing broken');

    if QPCONDENSINGSTEPS == 1
        mpc.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER');
    else
        mpc.set( 'SPARSE_QP_SOLUTION',  'BLOCK_CONDENSING_N2');
        mpc.set( 'CONDENSING_BLOCK_SIZE', QPCONDENSINGSTEPS  );
    end

    mpc.set( 'QP_SOLVER',               'QP_FORCES'          );


elseif strcmp(SOLVER,'FORCES')

    mpc.set( 'QP_SOLVER',               'QP_FORCES'          );
    mpc.set( 'SPARSE_QP_SOLUTION',      'SPARSE_SOLVER'      );


elseif contains(SOLVER,'HPMPC')

    mpc.set( 'QP_SOLVER',               'QP_HPMPC'           );
    mpc.set( 'SPARSE_QP_SOLUTION',      'SPARSE_SOLVER'      );

    if QPCONDENSINGSTEPS > 1
        mpc.set( 'CONDENSING_BLOCK_SIZE',  QPCONDENSINGSTEPS );
    end

else
    if ~strcmp(SOLVER, 'dfgm') && ~strcmp(SOLVER, 'osqp') && ~strcmp(SOLVER, 'fiordos') && ~contains(SOLVER, 'acados')
        error('SPECIFIED SOLVER DOES NOT EXIST')
    end
end

mpc.set( 'MAX_NUM_QP_ITERATIONS', 1000 );
mpc.set( 'PRINTLEVEL', 'LOW' );

if exist('export_MPC', 'dir') && MPC_EXPORT == 1 && MPC_COMPILE == 1
    rmdir('export_MPC', 's')
end

if MPC_EXPORT
    if is_acado_solver(SOLVER)
        mpc.exportCode( 'export_MPC' );
    end
    if CHECK_AGAINST_REF_SOL
        mpc.set( 'QP_SOLVER',               'QP_QPOASES'         );
        mpc.set( 'SPARSE_QP_SOLUTION',      'FULL_CONDENSING_N2' );
        mpc.exportCode( 'export_ref_MPC' );
    end
end

if MPC_COMPILE
    if is_acado_solver(SOLVER)
        cd export_MPC

        % add solver directory in export_MPC folder
        if contains(SOLVER,'qpDUNES')
            if QPCONDENSINGSTEPS == 0
                % use clipping
                waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'qpDUNES'], 'qpdunes'));
            else
                % use qpDUNES+qpOASES
                waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'qpDUNES-dev'], 'qpdunes'));
            end
        end
        if contains(SOLVER,'qpOASES_N')
            waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'acado-dev' filesep 'external_packages' filesep 'qpoases'], 'qpoases'));
        end
        if contains(SOLVER,'qpOASES_e_N')
            waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'acado-dev' filesep 'external_packages' filesep 'qpoases3'], 'qpoases3'));
        end
        if contains(SOLVER,'HPMPC')
            % TODO: either support HPMPC_old too (by adding submodule) or remove option in acado template
            waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'blasfeo'], 'blasfeo'));
            waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'hpmpc'], 'hpmpc'));
        end


        if contains(SOLVER, 'FORCES')
            % needed to give time to overwrite things
            % keyboard
            pause(10)
        end

        make_acado_solver('../acado_MPCstep')
        cd ..
    end
    if CHECK_AGAINST_REF_SOL
        cd export_ref_MPC
        waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'acado-dev' filesep 'external_packages' filesep 'qpoases'], 'qpoases'));
        make_acado_solver('../acado_ref_MPCstep')
        cd ..
    end
end


if strcmp(SOLVER, 'fiordos')
   fiordos_code_generate(N, NX, NU, W, WN, opts);
end

if strcmp(SOLVER, 'dfgm')
   dfgm_compile(N, NX, NU, opts);
end

%% Closed loop simulations

logged_outputs = cell(1, Tf/Ts);
solve_qp_times = nan(Tf/Ts, NRUNS);
simulate_times = nan(Tf/Ts, NRUNS);
solve_qp_iters = nan(Tf/Ts, NRUNS);
primal_feas    = nan(Tf/Ts, NRUNS);
dual_feas      = nan(Tf/Ts, NRUNS);
ref_sol_err    = nan(Tf/Ts, NRUNS);
ref_val_err    = nan(Tf/Ts, NRUNS);

for iRUNS = 1:NRUNS

    if strcmp(SOLVER, 'osqp')
        % TODO: FIX PROBLEM WITH INCONSISTENT NUMBER OF ITERATIONS ACROSS RUNS
        osqp_prob = osqp_setup(N, NX, NU, W, WN, opts);
    end

    X0   = fsolve_ref;
    Xref = repmat(fsolve_ref.',N+1,1);
    Uref = zeros(N,NU);

    % initialize input struct
    input.x    = repmat(X0.',N+1,1);
    input.u    = Uref;
    input.y    = [Xref(1:N,:) Uref];
    input.yN   = fsolve_ref.';
    input.W    = W;
    input.WN   = WN;
    input.WALL = WALL;
    input.LBU  = -ones(NU,1);
    input.UBU  = ones(NU,1);

    disp('------------------------------------------------------------------')
    disp('               Simulation Loop'                                    )
    disp('------------------------------------------------------------------')

    iter = 0;
    time = 0;

    controls_MPC = [];
    state_sim    = X0.';

    if VISUAL && iRUNS == 1
        visualize;
    end

    % clear mex memory for acado solver
    clear mex %#ok<CLMEX>

    while time(end) < Tf

        % update initial condition
        input.x0 = state_sim(end,:);

        % RTI step
        % TODO: merge to one call using function pointers and wrapping acado function
        if is_acado_solver(SOLVER)

            output = acado_MPCstep(input);

            % field name changed in future ACADO versions
            if isfield(output.info, 'QP_iter')
                output.info.nIterations = output.info.QP_iter;
            else
                output.info.nIterations = output.info.nIterations;
            end

            if contains(SOLVER, 'HPMPC')
                [output.info.primal_res, output.info.dual_res] = check_solution_accuracy_acado_hpmpc(input, output);
            end

        elseif contains(SOLVER, 'acados')

            input.nlp  = acados_nlp;
            input.time = time(end);
            output = acados_MPCstep(input);

        else

            switch SOLVER

                case 'dfgm'
                    output = dfgm_MPCstep(input, opts, time(end));
                case 'fiordos'
                    output = fiordos_MPCstep(input, opts, time(end));
                case 'osqp'
                    input.prob = osqp_prob;
                    output = osqp_MPCstep(input, opts, time(end));
            end
        end
        % TODO: check also slack_res

        % TODO: correct objVal of non acado solvers with constant term
        if CHECK_AGAINST_REF_SOL
            ref_output = acado_ref_MPCstep(input);
            sol_err    = max(norm(output.x(:) - ref_output.x(:), Inf), norm(output.u(:) - ref_output.u(:), Inf));
            val_err    = abs(output.info.objValue - ref_output.info.objValue)/max(1, ref_output.info.objValue);

            if sol_err > SOL_TOL || val_err > SOL_TOL
                warning(['failed to meet accuracy of ', num2str(SOL_TOL), '( sol_err = ', ...
                    num2str(sol_err), ', val_err = ', num2str(val_err), ')' ])
            end

            ref_sol_err(iter+1, iRUNS) = sol_err;
            ref_val_err(iter+1, iRUNS) = val_err;

            % save(['ws_' SOLVER '_N' num2str(length(time))], 'input', 'output', 'ref_output', 'sol_err');
        end

        if output.info.status ~= 1 &&  output.info.status ~= 0
            if time(end) < To
                warning('MPC STEP FAILED WHILE DISTURBANCE IS BEING APPLIED!')
            else
                warning('MPC STEP FAILED!')
            end
            output.info.cpuTime = NaN;
            keyboard
        end


        % LOG DATA
        if iRUNS == 1
            logged_outputs{iter+1} = output;
        end

        solve_qp_times(iter+1, iRUNS) = output.info.cpuTime;
        solve_qp_iters(iter+1, iRUNS) = output.info.nIterations;

        if isfield(output.info, 'simTime')
            simulate_times(iter+1, iRUNS) = output.info.simTime;
        else
            simulate_times(iter+1, iRUNS) = 0;
        end

        if isfield(output.info, 'primal_res')
            primal_feas(iter+1, iRUNS) = output.info.primal_res;
            dual_feas(iter+1, iRUNS)   = output.info.dual_res;
        else
            primal_feas(iter+1, iRUNS) = nan;
            dual_feas(iter+1, iRUNS)   = nan;
        end

        fprintf('%s:\t\t %d it\t %f ms\n', SOLVER, output.info.nIterations, 1000*output.info.cpuTime);

        % Save the MPC step
        controls_MPC = [controls_MPC; output.u(1,:)];

        % Optionally shift trajectories
        if ~ACADOS_BENCHMARK
            input.x = [output.x(2:end,:); output.x(end,:)];
            input.u = [output.u(2:end,:); output.u(end,:)];
        else
            % default init. in acados
            input.x = output.x;
            input.u = output.u;
        end

        % Simulate system
        sim_input.x = state_sim(end,:).';
        sim_input.u = output.u(1,:).';

        if time(end) < To
            sim_input.u = Uo;
            % do not take these instances into account for timings
            solve_qp_times(iter+1, iRUNS) = NaN;
            simulate_times(iter+1, iRUNS) = NaN;
        end

        [states, outputs] = eval([ sim_name,'(sim_input)']);

        state_sim = [state_sim; states.value'];
        iter      = iter + 1;
        nextTime  = iter*Ts;
        time      = [time nextTime];

        % check number of active constraints
        [nbx, nbu] = active_constraints(output, WALL, 1);
        disp(['nbx:' num2str(nbx) ' nbu:' num2str(nbu)])

        disp(['current time: ' num2str(nextTime) '   ' char(9) ' (RTI step: ' num2str(output.info.cpuTime*1e3) ' ms)'])

        if VISUAL && iRUNS == 1
            visualize;
        end
    end

    % delete last time instant
    time(end) = [];
end

if exist('sim_opts','var')
    figure
    plot(time,1000*min(solve_qp_times, [], 2),'lineWidth',2)
    title('Timings of RTI scheme in closed loop','fontSize',16,'fontWeight','Normal')
    xlabel('Time [s]','fontSize',16)
    ylabel('CPU time [ms]','fontSize',16)
    leg = legend(SOLVER);
    leg.FontSize = 16;
    set(gca,'fontSize',16)
end

% store simulation information
logged_data.wall        = WALL;
logged_data.Ts          = Ts;
logged_data.N           = N;
logged_data.Nmass       = NMASS;
logged_data.nruns       = NRUNS;
logged_data.outputs     = logged_outputs;
logged_data.Nblock      = QPCONDENSINGSTEPS;
logged_data.solver      = SOLVER;
logged_data.cputime     = solve_qp_times;
logged_data.simtime     = simulate_times;
logged_data.iters       = solve_qp_iters;
logged_data.primal_feas = primal_feas;
logged_data.dual_feas   = dual_feas;

if CHECK_AGAINST_REF_SOL
    logged_data.ref_sol_err = ref_sol_err;
    logged_data.ref_val_err = ref_val_err;
end

if ~BATCH

    if CHECK_AGAINST_REF_SOL
        disp(' ');
        disp(['MAX ERROR IN SOLUTION: ' num2str(max(max(logged_data.ref_sol_err)))]);
    end

    disp('ERRORS:')
    [max(max(primal_feas)) max(max(dual_feas)) max(max(ref_sol_err))]

    % consistency checks
    if ~check_consistency(solve_qp_iters)
        warning('inconsistent iterations between runs')
    end
    if ~any(any(isnan(primal_feas)))
        if ~check_consistency(primal_feas) || ~check_consistency(dual_feas)
            warning('inconsistent solution accuracy between runs')
        end
    end
    keyboard
end

end
