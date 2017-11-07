function logging = NMPC_chain_mass(sim_opts)

%% INITIALIZE

if ~exist('sim_opts','var')
    clear all
end

clear global % for ACADO workspace
clearvars -except sim_opts; clc; close all; % sim_opts contain options to overwrite current script
close all
clc

%% SIMULATION OPTIONS

SIM_EXPORT    = 1;              % export code for ACADO simulator

SIM_COMPILE   = 1;              % compile exported code for ACADO simulator

MPC_EXPORT    = 1;              % export code for ACADO solver

MPC_COMPILE   = 1;              % compile exported code for ACADO solver

NRUNS         = 5;              % run closed-loop simulation NRUNS times and store minimum timings (to minimize OS interference)

ACADOSOLVER   = 'qpDUNES_B0';   % 'qpDUNES_BXX' (with XX block size, 0 for clipping), 'qpOASES_N2', 'qpOASES_N3', 'FORCES', 'HPMPC'

VISUAL        = 1;              % set to 1 to visualize chain of masses (only for the first out of the NRUNS simulations)

WALL          = -0.05;          % wall position (re-export if changed)

Ts            = 0.1;            % sampling time [s]

Tf            = 5;              % final simulation time [s]

N             = 40;             % prediction horizon

NMASS         = 6;              % number of masses in chain (data available from 3 to 6 masses)

INITMODE      = 1;              % 1: initialize lin. at X0
                                % 2: initialize lin. at XREF
                                % NOTE: does not make a difference if SCENARIO == 2

SCENARIO      = 2;              % 1: read x0 from file and go to XREF
                                % 2: start from XREF and overwrite uMPC until To [s]

To            = 0.5;            % how many seconds to overwrite uMPC when SCENARIO == 2
Uo            = [-1;1;0];       % value to overwrite uMPC when SCENARIO == 2

RECORD_VID    = 0;              % record avi video

DETAILED_TIME = 0;              % if 1, time preparation/feedback step: ONLY WORKS FOR FORCES


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

% prepare video writer
if RECORD_VID
    % TODO why doesn't it work? (at leat on mac)
    writerObj = VideoWriter('chain_video.avi');
    writerObj.FrameRate = 5;
    open(writerObj);
end

% extract block size from solver name
if contains(ACADOSOLVER,'qpDUNES')
    bpos = strfind(ACADOSOLVER,'B');
    QPCONDENSINGSTEPS = str2double(ACADOSOLVER(bpos+1:end));
    if QPCONDENSINGSTEPS ~= 0 &&(QPCONDENSINGSTEPS < 1 || mod(N,QPCONDENSINGSTEPS)~= 0)
        error('Invalid block size for given horizon length N.')
    end
else
    QPCONDENSINGSTEPS = [];
end

M  = NMASS - 2;     % number of intermediate masses
NX = (2*M + 1)*3;   % differential states
NU = 3;             % control inputs

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

% Compute the spring forces:
g = acado.Expression([0; 0; -9.81]);
f = is(repmat(g, M, 1));
Force = [];
for i = 1:M+1
    if i == 1
        dist = is(x((i-1)*3+1:i*3) - x0);
    elseif( i <= M )
        dist = is(x((i-1)*3+1:i*3) - x((i-2)*3+1:(i-1)*3));
    else
        dist = is(xEnd - x((M-1)*3+1:end));
    end

    scale = D/m*(1-L/norm(dist));
    F = is(scale*dist);

    Force = [Force; F];

    % mass on the right
    if i < M+1
        f((i-1)*3+1:i*3) = f((i-1)*3+1:i*3) - F;
    end
    % mass on the left
    if i > 1
        f((i-2)*3+1:(i-1)*3) = f((i-2)*3+1:(i-1)*3) + F;
    end
end

ode = [ dot(x) == v; ...
        dot(v) == f ];

%% SIMexport

acadoSet('problemname', 'sim');

sim = acado.SIMexport( Ts );
sim.setLinearInput(A1,B1);
sim.setModel(ode);
sim.set( 'INTEGRATOR_TYPE',        'INT_IRK_GL2' );
sim.set( 'NUM_INTEGRATOR_STEPS',        2        );

if SIM_EXPORT
    sim.exportCode( 'export_SIM' );
end
if SIM_COMPILE
    cd export_SIM
    make_acado_integrator('../integrate_chain')
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

ocp.subjectTo( WALL <= [x([2:3:end]); xEnd(2)] <= 100 ); % constraint on y-position TODO upper bound 10?
ocp.subjectTo( -1 <= u <= 1 );                          % box constraints on controls

ocp.setLinearInput(A1,B1);
ocp.setModel(ode);

mpc = acado.OCPexport( ocp );

mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'       );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING'  );

if strcmp(ACADOSOLVER,'qpOASES_N3')

    mpc.set( 'QP_SOLVER',               'QP_QPOASES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'CONDENSING'         );

elseif strcmp(ACADOSOLVER,'qpOASES_N2')

    mpc.set( 'QP_SOLVER',               'QP_QPOASES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'FULL_CONDENSING_N2' );

elseif contains(ACADOSOLVER,'qpDUNES') && QPCONDENSINGSTEPS <= 1

    mpc.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'SPARSE_SOLVER'      );

elseif contains(ACADOSOLVER,'qpDUNES') && QPCONDENSINGSTEPS > 1

    mpc.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'BLOCK_CONDENSING_N2');
    mpc.set( 'CONDENSING_BLOCK_SIZE',    QPCONDENSINGSTEPS   );

elseif strcmp(ACADOSOLVER,'FORCES_BC') % TODO FIX

    if QPCONDENSINGSTEPS == 1
        mpc.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER');
    else
        mpc.set( 'SPARSE_QP_SOLUTION',  'BLOCK_CONDENSING_N2');
        mpc.set( 'CONDENSING_BLOCK_SIZE', QPCONDENSINGSTEPS  );
    end

    mpc.set( 'QP_SOLVER',               'QP_FORCES'          );


elseif strcmp(ACADOSOLVER,'FORCES')

    mpc.set( 'QP_SOLVER',           'QP_FORCES'              );
    mpc.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER'          );

elseif strcmp(ACADOSOLVER,'HPMPC')

    mpc.set( 'QP_SOLVER',      'QP_HPMPC'                    );
    mpc.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER'          );

else

    error('SPECIFIED SOLVER DOES NOT EXIST')

end

mpc.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL2'        );
mpc.set( 'NUM_INTEGRATOR_STEPS',        2*N                  );

mpc.set( 'MAX_NUM_QP_ITERATIONS', 2000 );
mpc.set( 'PRINTLEVEL', 'LOW' );

if exist('export_MPC', 'dir')
    rmdir('export_MPC', 's')
end
if MPC_EXPORT
    mpc.exportCode( 'export_MPC' );
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
    if contains(ACADOSOLVER,'qpOASES')
        waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'acado-dev' filesep 'external_packages' filesep 'qpoases'], 'qpoases'));
    end
    if contains(ACADOSOLVER,'HPMPC')
        % TODO: either support HPMPC_old too (by adding submodule) or remove option in acado template
        waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'blasfeo'], 'blasfeo'));
        waitfor(copyfile(['..' filesep '..' filesep 'external' filesep 'hpmpc'], 'hpmpc'));
    end
    
    
    if contains(ACADOSOLVER, 'FORCES')
        % needed to give time to overwrite things
        keyboard
    end
    if DETAILED_TIME
        make_timing_forces('../acado_MPCstep')
    else
        make_acado_solver('../acado_MPCstep')
    end
    cd ..
end

%% Closed loop simulations

minACADOtLog = []; % minimum timings of solver over NRUNS

for iRUNS = 1:NRUNS

    eval(['ref = textread(' '''' 'chain_mass' filesep 'chain_mass_model_eq_M' num2str(NMASS) '.txt' '''', ', ''''' ');']);
    ref  = [ref(end-3+1:end); ref(1:end-3)]; % fix ordering convention

    if SCENARIO == 1
        eval(['X0  = textread(' '''' 'chain_mass' filesep 'chain_mass_model_dist_M' num2str(NMASS) '.txt' '''', ', ''''' ');']);
        X0 = [X0(end-3+1:end); X0(1:end-3)]; % fix ordering convention
    elseif SCENARIO == 2
        X0 = ref;
    end

    Xref = repmat(ref.',N+1,1);
    Uref = zeros(N,NU);

    if INITMODE == 1
        input.x = repmat(X0.',N+1,1);
    elseif INITMODE == 2
        input.x = repmat(ref.',N+1,1);
    else
        error('wrong initialization flag INITMODE')
    end

    input.u  = Uref;
    input.y  = [Xref(1:N,:) Uref];
    input.yN = ref.';

    disp('------------------------------------------------------------------')
    disp('               Simulation Loop'                                    )
    disp('------------------------------------------------------------------')

    iter = 0;
    time = 0;

    controls_MPC = [];
    state_sim    = X0.';

    ACADOtLog     = [];  % log timings of solver
    ACADOoutputs  = {};  % log all ACADO outputs
    ACADOnIter    = [];  % log iterations (if available)

    if DETAILED_TIME
        ACADOtprepLog = [];  % log preparation times
    end

    if VISUAL && iRUNS == 1
        visualize;
    end

    % weights % TODO MOVE UP TO OPTIONS
    input.W  = blkdiag(15*eye(3), 10*eye(length(x)), eye(length(v)), 0.05*eye(3));
    input.WN = blkdiag(15*eye(3), 10*eye(length(x)), eye(length(v)));

    if NRUNS > 1
        clear mex %#ok<CLMEX>
    end

    while time(end) < Tf

        % Solve NMPC OCP with ACADO
        input.x0 = state_sim(end,:);
        output   = acado_MPCstep(input);

        if output.info.status ~= 1 &&  output.info.status ~= 0
            warning('MPC STEP FAILED')
            output.info.cpuTime = NaN;
            keyboard
        end

        if isfield(output.info, 'QP_iter')
            niter = output.info.QP_iter;
        else
            niter = output.info.nIterations;
        end

        ACADOnIter = [ACADOnIter niter];
        ACADOoutputs{end+1} = output;
        ACADOtLog = [ACADOtLog; output.info.cpuTime];

        if DETAILED_TIME
            ACADOtprepLog = [ACADOtprepLog; output.info.preparationTime];
        end

        % Save the MPC step
        controls_MPC = [controls_MPC; output.u(1,:)];

        % Shift trajectories
        input.x = [output.x(2:end,:); output.x(end,:)];
        input.u = [output.u(2:end,:); output.u(end,:)];

        % Simulate system
        sim_input.x = state_sim(end,:).';
        sim_input.u = output.u(1,:).';

        if SCENARIO == 2
            if time(end) < To
                sim_input.u = Uo;
            end
        end
        [states,outputs] = integrate_chain(sim_input);

        state_sim = [state_sim; states.value'];
        iter      = iter + 1;
        nextTime  = iter*Ts;
        time      = [time nextTime];

        if DETAILED_TIME
            disp(['current time: ' num2str(nextTime) '   ' char(9) ' (pre. step: ' num2str(output.info.preparationTime*1e3) ' ms)'])
        else
            disp(['current time: ' num2str(nextTime) '   ' char(9) ' (RTI step: ' num2str(output.info.cpuTime*1e3) ' ms)'])
        end

        if VISUAL && iRUNS == 1
            visualize;
            if RECORD_VID
                if ~exist('frame','var')
                    frame(1) = getframe(gcf);
                else
                    frame(end+1) = getframe(gcf);
                end
                writeVideo(writerObj, frame(end));
            end
        end

    end

    % delete last time instant
    time(end) = [];

    if iRUNS > 1
        minACADOtLog = min(ACADOtLog, minACADOtLog);
        if DETAILED_TIME
            minACADOtprepLog = min(ACADOtprepLog, minACADOtprepLog);
        end
    else
        minACADOtLog = ACADOtLog;
        if DETAILED_TIME
            minACADOtprepLog = ACADOtprepLog;
        end
    end

end

if RECORD_VID
    close(writerObj)
    % close all
    % movie(frame);
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
logging.solver  = ACADOSOLVER;
logging.wall    = WALL;
logging.Ts      = Ts;
logging.N       = N;
logging.Nmass   = NMASS;
logging.init    = INITMODE;
logging.nruns   = NRUNS;
logging.cputime = minACADOtLog;
logging.iters   = ACADOnIter;
logging.outputs = ACADOoutputs;
logging.Nblock  = QPCONDENSINGSTEPS;

if DETAILED_TIME
   logging.prepTime = ACADOtprepLog;
end

end
