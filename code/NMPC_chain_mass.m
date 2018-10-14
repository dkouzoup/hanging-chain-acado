function logged_data = NMPC_chain_mass(sim_opts)

%% INITIALIZE

if ~exist('sim_opts','var')
    clear all
end

clear global % for ACADO workspace
clearvars -except sim_opts; clc; close all; % sim_opts contain options to overwrite current script
close all
clc

addpath([pwd filesep 'utils'])

%% SIMULATION OPTIONS

SIM_EXPORT    = 0;              % export code for ACADO simulator

SIM_COMPILE   = 0;              % compile exported code for ACADO simulator

MPC_EXPORT    = 0;              % export code for ACADO solver

MPC_COMPILE   = 0;              % compile exported code for ACADO solver

NRUNS         = 5;              % run closed-loop simulation NRUNS times and store minimum timings (to minimize OS interference)

ACADOSOLVER   = 'qpOASES_N2';   % 'qpDUNES_BXX' (with XX block size, 0 for clipping), 'qpOASES_N2', 'qpOASES_e_N3', 'qpOASES_e_N2', 'qpOASES_N3', 'FORCES', 'HPMPC'

WARMSTART     = 0;              % applicable for qpOASES/qpDUNES

VISUAL        = 0;              % set to 1 to visualize chain of masses (only for the first out of the NRUNS simulations)

WALL          = -0.1;           % wall position (re-export if changed)

Ts            = 0.2;            % sampling time [s]

Tf            = 5;              % final simulation time [s]

N             = 10;             % prediction horizon

NMASS         = 3;              % number of masses in chain (data available from 3 to 6 masses)

To            = 5*Ts;           % how many seconds to overwrite uMPC

Uo            = [-1;1;1];       % value to overwrite uMPC for To sec

DETAILED_TIME = 0;              % if 1, time preparation/feedback step: ONLY WORKS FOR FORCES

CHECK_AGAINST_REF_SOL = 0;      % if 1, exports and compiles reference solver (qpOASES, CN2 by default)

SOL_TOL = 1e-6;                 % maximum accepted 2-norm of the deviation of the solution from the reference solution
                                % (only used if CHECK_AGAINST_REF_SOL = 1).

SECONDARY_SOLVER = 'osqp';   % empty, 'fiordos', 'dfgm' or 'osqp'


% dimensions
M  = NMASS - 2;     % number of intermediate masses
NX = (2*M + 1)*3;   % differential states
NU = 3;             % control inputs

% MPC weights (may be hard-coded for some solvers)
W  = blkdiag(25*eye(3), 25*eye(3*M), 1*eye(3*M), 0.01*eye(NU));
WN = blkdiag(25*eye(3), 25*eye(3*M), 1*eye(3*M));

if ~isempty(SECONDARY_SOLVER)
   
    switch SECONDARY_SOLVER
        
        case 'fiordos'
            
            sec_opts.approach  = 'dual'; % 'dual' or 'primal-dual'
            sec_opts.tol       = 1e-2;
            sec_opts.maxit     = 100000;
            sec_opts.export    = 1;
            sec_opts.compile   = 1;
            sec_opts.infval    = 1e8;
            sec_opts.warmstart = 1; % 0: no warmstart, 1: same solution
            
        case 'dfgm'
            
            sec_opts.warmstart = 2; % 0: no warmstart, 1: same solution, 2: shifted solution
            sec_opts.maxit     = 100000;
            sec_opts.tol       = 1e-2;
            
        case 'osqp'

            sec_opts.warmstart = 1;
            sec_opts.maxit     = 1000;
            sec_opts.check_ter = 5;
    end
    
end


%% Load simulation options and overwrite local ones

if exist('sim_opts','var')
    names  = fieldnames(sim_opts);
    for ii = 1:length(names)
        eval([names{ii} ' = sim_opts.' names{ii} ';'])
    end
end

if DETAILED_TIME == 1 && ~strcmp(ACADOSOLVER, 'FORCES')
    error('detailed timings only implemented for FORCES solver');
end

%% Initialization

% extract block size from solver name
if contains(ACADOSOLVER,'qpDUNES')
    bpos = strfind(ACADOSOLVER,'B');
    QPCONDENSINGSTEPS = str2double(ACADOSOLVER(bpos+1:end));
    if QPCONDENSINGSTEPS ~= 0 &&(QPCONDENSINGSTEPS < 1 || mod(N,QPCONDENSINGSTEPS)~= 0)
        error('Invalid block size for given horizon length N.')
    end
elseif contains(ACADOSOLVER,'HPMPC')
    bpos = strfind(ACADOSOLVER,'B');
    QPCONDENSINGSTEPS = str2double(ACADOSOLVER(bpos+1:end));
    if ~(QPCONDENSINGSTEPS >=0)
        error('Invalid block size (block size has to be >= 0)')
    end
else
    QPCONDENSINGSTEPS = [];
end

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

%% SIMexport

acadoSet('problemname', 'sim');

sim = acado.SIMexport(Ts);
sim.setLinearInput(A1, B1);
sim.setModel(ode);
sim.set( 'INTEGRATOR_TYPE',        'INT_IRK_GL2' );
sim.set( 'NUM_INTEGRATOR_STEPS',        2        );

if SIM_EXPORT
    sim.exportCode( 'export_SIM' );
end

sim_name = sprintf('integrate_chain_M%d',NMASS);
if SIM_COMPILE
    cd export_SIM
    sim_path = '../';
    make_acado_integrator([sim_path sim_name])
    if ~isempty(SECONDARY_SOLVER)
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

ocp.setLinearInput(A1,B1);
ocp.setModel(ode);

mpc = acado.OCPexport( ocp );

mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'       );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING'  );

if strcmp(ACADOSOLVER,'qpOASES_N3')

    mpc.set( 'QP_SOLVER',               'QP_QPOASES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'CONDENSING'         );
    if WARMSTART
        mpc.set( 'HOTSTART_QP',         'YES'                );
    end

elseif strcmp(ACADOSOLVER,'qpOASES_N2')

    mpc.set( 'QP_SOLVER',               'QP_QPOASES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'FULL_CONDENSING_N2' );
    if WARMSTART
        mpc.set( 'HOTSTART_QP',         'YES'                );
    end

elseif strcmp(ACADOSOLVER,'qpOASES_e_N3')

    mpc.set( 'QP_SOLVER',               'QP_QPOASES3'        );
    mpc.set( 'SPARSE_QP_SOLUTION',      'CONDENSING'         );
    if WARMSTART
        mpc.set( 'HOTSTART_QP',         'YES'                );
    end

elseif strcmp(ACADOSOLVER,'qpOASES_e_N2')

    mpc.set( 'QP_SOLVER',               'QP_QPOASES3'        );
    mpc.set( 'SPARSE_QP_SOLUTION',      'FULL_CONDENSING_N2' );
    if WARMSTART
        mpc.set( 'HOTSTART_QP',         'YES'                );
    end

elseif contains(ACADOSOLVER,'qpDUNES') && QPCONDENSINGSTEPS <= 1

    mpc.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'SPARSE_SOLVER'      );
    if WARMSTART % TODO: IS THIS USED ANYWHERE IN ACADO?!
        mpc.set( 'HOTSTART_QP',         'YES'                );
    end

elseif contains(ACADOSOLVER,'qpDUNES') && QPCONDENSINGSTEPS > 1

    mpc.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'BLOCK_CONDENSING_N2');
    mpc.set( 'CONDENSING_BLOCK_SIZE',    QPCONDENSINGSTEPS   );
    if WARMSTART == 0
        warning('qpDUNES with cold start not implemented in ACADO');
        keyboard
    end
elseif strcmp(ACADOSOLVER,'FORCES_BC') % TODO FIX

    if QPCONDENSINGSTEPS == 1
        mpc.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER');
    else
        mpc.set( 'SPARSE_QP_SOLUTION',  'BLOCK_CONDENSING_N2');
        mpc.set( 'CONDENSING_BLOCK_SIZE', QPCONDENSINGSTEPS  );
    end

    mpc.set( 'QP_SOLVER',               'QP_FORCES'          );


elseif strcmp(ACADOSOLVER,'FORCES')

    mpc.set( 'QP_SOLVER',               'QP_FORCES'          );
    mpc.set( 'SPARSE_QP_SOLUTION',      'SPARSE_SOLVER'      );


elseif contains(ACADOSOLVER,'HPMPC')

    mpc.set( 'QP_SOLVER',               'QP_HPMPC'           );
    mpc.set( 'SPARSE_QP_SOLUTION',      'SPARSE_SOLVER'      );

    if QPCONDENSINGSTEPS > 1
        mpc.set( 'CONDENSING_BLOCK_SIZE',  QPCONDENSINGSTEPS );
    end

else

    error('SPECIFIED SOLVER DOES NOT EXIST')

end

mpc.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL2'        );
mpc.set( 'NUM_INTEGRATOR_STEPS',        2*N                  );

mpc.set( 'MAX_NUM_QP_ITERATIONS', 1000 );
mpc.set( 'PRINTLEVEL', 'LOW' );

if exist('export_MPC', 'dir') && MPC_EXPORT == 1 && MPC_COMPILE == 1
    rmdir('export_MPC', 's')
end

if MPC_EXPORT
    mpc.exportCode( 'export_MPC' );
    if CHECK_AGAINST_REF_SOL
        mpc.set( 'QP_SOLVER',               'QP_QPOASES'         );
        mpc.set( 'SPARSE_QP_SOLUTION',      'FULL_CONDENSING_N2' );
        mpc.exportCode( 'export_ref_MPC' );
    end
end

if MPC_COMPILE
    cd export_MPC

    % add solver directory in export_MPC folder
    if contains(ACADOSOLVER,'qpDUNES')
        if QPCONDENSINGSTEPS == 0
            % use clipping
            waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'qpDUNES'], 'qpdunes'));
        else
            % use qpDUNES+qpOASES
            waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'qpDUNES-dev'], 'qpdunes'));
        end
    end
    if contains(ACADOSOLVER,'qpOASES_N')
        waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'acado-dev' filesep 'external_packages' filesep 'qpoases'], 'qpoases'));
    end
    if contains(ACADOSOLVER,'qpOASES_e_N')
        waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'acado-dev' filesep 'external_packages' filesep 'qpoases3'], 'qpoases3'));
    end
    if contains(ACADOSOLVER,'HPMPC')
        % TODO: either support HPMPC_old too (by adding submodule) or remove option in acado template
        waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'blasfeo'], 'blasfeo'));
        waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'hpmpc'], 'hpmpc'));
    end


    if contains(ACADOSOLVER, 'FORCES')
        % needed to give time to overwrite things
        % keyboard
        pause(10)
    end
    if DETAILED_TIME
        make_timing_forces('../acado_MPCstep')
    else
        make_acado_solver('../acado_MPCstep')
    end
    cd ..

    if CHECK_AGAINST_REF_SOL
        cd export_ref_MPC
        waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'acado-dev' filesep 'external_packages' filesep 'qpoases'], 'qpoases'));
        make_acado_solver('../acado_ref_MPCstep')
        cd ..
    end
end


if strcmp(SECONDARY_SOLVER, 'fiordos')
   fiordos_code_generate(N, NX, NU, W, WN, sec_opts); 
end

if strcmp(SECONDARY_SOLVER, 'dfgm')
   dfgm_compile(N, NX, NU); 
end

if strcmp(SECONDARY_SOLVER, 'osqp')
   osqp_prob = osqp_setup(N, NX, NU, W, WN, sec_opts); 
end

%% Closed loop simulations

% minimum timings of solver over NRUNS
minACADOtLog = [];
if ~isempty(SECONDARY_SOLVER)
    secondary_solve_qp_min_times = [];
end

for iRUNS = 1:NRUNS
    
    X0   = fsolve_ref;
    Xref = repmat(fsolve_ref.',N+1,1);
    Uref = zeros(N,NU);

    input.x = repmat(X0.',N+1,1);
    input.u  = Uref;
    input.y  = [Xref(1:N,:) Uref];
    input.yN = fsolve_ref.';

    disp('------------------------------------------------------------------')
    disp('               Simulation Loop'                                    )
    disp('------------------------------------------------------------------')

    iter = 0;
    time = 0;

    controls_MPC = [];
    state_sim    = X0.';

    ACADOtLog     = [];  % log timings of solver
    ACADOtSimLog  = [];  % log simulation timings of solver
    ACADOoutputs  = {};  % log all ACADO outputs
    ACADOnIter    = [];  % log iterations (if available)
    sol_accuracy  = [];  % log accuracy (deviation from reference solution) of the solution (if available)
    val_accuracy  = [];  % log accuracy (deviation from reference obj value) of the solution (if available)

    if ~isempty(SECONDARY_SOLVER)
        secondary_solve_qp_tmp_times = [];
        secondary_solve_qp_iters     = [];
        secondary_solve_qp_error     = [];
    end

    if DETAILED_TIME
        ACADOtprepLog = [];  % log preparation times
    end

    if VISUAL && iRUNS == 1
        visualize;
    end

    % weights
    % input.W  = blkdiag(15*eye(3), 10*eye(length(x)), eye(length(v)), 0.05*eye(3));
    % input.WN = blkdiag(15*eye(3), 10*eye(length(x)), eye(length(v)));
    input.W  = W;
    input.WN = WN;

    % clear mex memory for acado solver
    clear mex %#ok<CLMEX>

    while time(end) < Tf

        % Solve NMPC OCP with ACADO
        input.x0 = state_sim(end,:);
        output   = acado_MPCstep(input);

        if ~isempty(SECONDARY_SOLVER)
            
            sec_input = input;
            sec_input.acado_sol = output;
            sec_input.WALL    = WALL;
            sec_input.LBU     = -ones(NU,1);
            sec_input.UBU     = ones(NU,1);
            
            switch SECONDARY_SOLVER
                
                case 'fiordos'
                    sec_output = fiordos_MPCstep(sec_input, sec_opts, time(end));
                    
                case 'dfgm'
                  
                    sec_output = dfgm_MPCstep(sec_input, sec_opts, time(end));
                  
                case 'osqp'
                    sec_input.prob = osqp_prob;
                    sec_output = osqp_MPCstep(sec_input, sec_opts, time(end));
            end
            
        end

        if CHECK_AGAINST_REF_SOL
            ref_output = acado_ref_MPCstep(input);
            sol_err = max(norm(output.x - ref_output.x, Inf), norm(output.u - ref_output.u, Inf));
            val_err = abs(output.info.objValue - ref_output.info.objValue)/max(1, ref_output.info.objValue);
            if sol_err > SOL_TOL || val_err > SOL_TOL
%                 keyboard
                warning(['failed to meet accuracy of ', num2str(SOL_TOL), '( sol_err = ', ...
                    num2str(sol_err), ', val_err = ', num2str(val_err), ')' ])
            end

            sol_accuracy = [sol_accuracy sol_err];
            val_accuracy = [val_accuracy val_err];
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

        % field name changed in future ACADO versions
        if isfield(output.info, 'QP_iter')
            niter = output.info.QP_iter;
        else
            niter = output.info.nIterations;
        end

        ACADOnIter = [ACADOnIter niter];
        ACADOoutputs{end+1} = output;
        ACADOtLog = [ACADOtLog; output.info.cpuTime];
        if isfield(output.info, 'simTime')
            ACADOtSimLog = [ACADOtSimLog; output.info.simTime];
        else
            ACADOtSimLog = [ACADOtSimLog; inf];
        end

        if DETAILED_TIME
            ACADOtprepLog = [ACADOtprepLog; output.info.preparationTime];
        end

        if ~isempty(SECONDARY_SOLVER)
            secondary_solve_qp_tmp_times(end+1) = sec_output.info.QP_time;
            secondary_solve_qp_iters(end+1)     = sec_output.info.nIterations;
            secondary_solve_qp_error(end+1)     = norm([sec_output.x(:); sec_output.u(:)] - [output.x(:); output.u(:)], inf);

            fprintf('ACADO:\t\t %d it\t %f ms\n', niter, 1000*(ACADOtLog(end) - ACADOtSimLog(end)));
            fprintf('%s:\t\t %d it\t %f ms\n', SECONDARY_SOLVER, sec_output.info.nIterations, 1000*sec_output.info.QP_time);
            fprintf('solution gap:\t %5.e\n', norm([sec_output.x(:); sec_output.u(:)] - [output.x(:); output.u(:)], inf));
        end

        % Save the MPC step
        controls_MPC = [controls_MPC; output.u(1,:)];

        % Shift trajectories
        input.x = [output.x(2:end,:); output.x(end,:)];
        input.u = [output.u(2:end,:); output.u(end,:)];

        % Simulate system
        sim_input.x = state_sim(end,:).';
        sim_input.u = output.u(1,:).';

        if time(end) < To
            sim_input.u = Uo;
            % do not take these instances into account for timings
            ACADOtLog(end) = NaN;
            ACADOtSimLog(end) = NaN;
            if ~isempty(SECONDARY_SOLVER)
                secondary_solve_qp_tmp_times(end) = NaN;
            end
        end

        [states, outputs] = eval([ sim_name,'(sim_input)']);

        state_sim = [state_sim; states.value'];
        iter      = iter + 1;
        nextTime  = iter*Ts;
        time      = [time nextTime];

        % check number of active constraints
        [nbx, nbu] = active_constraints(output, WALL, 1);
        disp(['nbx:' num2str(nbx) ' nbu:' num2str(nbu)])

        if DETAILED_TIME
            disp(['current time: ' num2str(nextTime) '   ' char(9) ' (pre. step: ' num2str(output.info.preparationTime*1e3) ' ms)'])
        else
            disp(['current time: ' num2str(nextTime) '   ' char(9) ' (RTI step: ' num2str(output.info.cpuTime*1e3) ' ms)'])
        end

        if VISUAL && iRUNS == 1
            visualize;
        end

    end

    % delete last time instant
    time(end) = [];

    if iRUNS == 1
        minACADOtLog = ACADOtLog;
        minACADOtSimLog = ACADOtSimLog;
        if DETAILED_TIME
            minACADOtprepLog = ACADOtprepLog;
        end
        if ~isempty(SECONDARY_SOLVER)
            secondary_solve_qp_min_times = secondary_solve_qp_tmp_times;
        end
    else
        minACADOtLog = min(ACADOtLog, minACADOtLog);
        minACADOtSimLog = min(ACADOtSimLog, minACADOtSimLog);

        if DETAILED_TIME
            minACADOtprepLog = min(ACADOtprepLog, minACADOtprepLog);
        end
        if ~isempty(SECONDARY_SOLVER)
            secondary_solve_qp_min_times = min(secondary_solve_qp_tmp_times, secondary_solve_qp_min_times);
        end
    end
end

if exist('sim_opts','var')
    figure
    plot(time,1000*minACADOtLog,'lineWidth',2)
    title('Timings of RTI scheme in closed loop','fontSize',16,'fontWeight','Normal')
    xlabel('Time [s]','fontSize',16)
    ylabel('CPU time [ms]','fontSize',16)
    leg = legend(ACADOSOLVER);
    leg.FontSize = 16;
    set(gca,'fontSize',16)
end

% store simulation information
logged_data.wall    = WALL;
logged_data.Ts      = Ts;
logged_data.N       = N;
logged_data.Nmass   = NMASS;
logged_data.nruns   = NRUNS;
logged_data.cputime = minACADOtLog;
logged_data.simtime = minACADOtSimLog;
logged_data.iters   = ACADOnIter;
logged_data.outputs = ACADOoutputs;
logged_data.Nblock  = QPCONDENSINGSTEPS;
logged_data.sol_accuracy = sol_accuracy;
logged_data.val_accuracy = val_accuracy;

if DETAILED_TIME
   logged_data.prepTime = ACADOtprepLog;
end

if ~isempty(SECONDARY_SOLVER)
   logged_data.secondary_qptime     = secondary_solve_qp_min_times';
   logged_data.secondary_iter       = secondary_solve_qp_iters;
   logged_data.secondary_error_sol  = secondary_solve_qp_error;

   if 1
       disp(['MAX ERROR IN SOLUTION:          ' num2str(max(abs(logged_data.secondary_error_sol)))]);
       close all
%        plot(minACADOtLog-minACADOtSimLog);
%        hold on
%        plot(fiordos_solve_qp_min_times);
%        figure
       plot(minACADOtLog-minACADOtSimLog)
       hold on
       plot(secondary_solve_qp_min_times')
       title('CPU times')
       legend('acado', SECONDARY_SOLVER)
       [(minACADOtLog-minACADOtSimLog) secondary_solve_qp_min_times' logged_data.secondary_error_sol']
       logged_data.secondary_iter'
       keyboard
   end
end

end
