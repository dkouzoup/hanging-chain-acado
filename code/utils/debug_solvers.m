function logging = debug(sim_opts)

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

NRUNS         = 1;              % run closed-loop simulation NRUNS times and store minimum timings (to minimize OS interference)

ACADOSOLVER1  = 'HPMPC';       % 'qpDUNES_BXX' (with XX block size, 0 for clipping), 'qpOASES_N2', 'qpOASES_N3', 'FORCES', 'HPMPC'
ACADOSOLVER2  = 'qpOASES_N2';   % solve with both solvers, compare solutions, use output of second solver

VISUAL        = 1;              % set to 1 to visualize chain of masses (only for the first out of the NRUNS simulations)

WALL          = -0.05;          % wall position (re-export if changed)

Ts            = 0.1;            % sampling time [s]

Tf            = 5;              % final simulation time [s]

N             = 100;             % prediction horizon

NMASS         = 4;              % number of masses in chain (data available from 3 to 6 masses)

INITMODE      = 1;              % 1: initialize lin. at X0
                                % 2: initialize lin. at XREF
                                % NOTE: does not make a difference if SCENARIO == 2

SCENARIO      = 2;              % 1: read x0 from file and go to XREF
                                % 2: start from XREF and overwrite uMPC until To [s]

To            = 0.5;            % how many seconds to overwrite uMPC when SCENARIO == 2
Uo            = [-1;1;0];       % value to overwrite uMPC when SCENARIO == 2

RECORD_VID    = 0;              % record avi video

% Load simulation options and overwrite local ones
if exist('sim_opts','var')
    names  = fieldnames(sim_opts);
    for ii = 1:length(names)
        eval([names{ii} ' = sim_opts.' names{ii} ';'])
    end
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
if contains(ACADOSOLVER1,'qpDUNES')
    bpos = strfind(ACADOSOLVER1,'B');
    QPCONDENSINGSTEPS = str2double(ACADOSOLVER1(bpos+1:end));
    if QPCONDENSINGSTEPS ~= 0 &&(QPCONDENSINGSTEPS < 1 || mod(N,QPCONDENSINGSTEPS)~= 0)
        error('Invalid block size for given horizon length N.')
    end
else
    QPCONDENSINGSTEPS1 = [];
end
if contains(ACADOSOLVER2,'qpDUNES')
    bpos = strfind(ACADOSOLVER2,'B');
    QPCONDENSINGSTEPS = str2double(ACADOSOLVER2(bpos+1:end));
    if QPCONDENSINGSTEPS ~= 0 &&(QPCONDENSINGSTEPS < 1 || mod(N,QPCONDENSINGSTEPS)~= 0)
        error('Invalid block size for given horizon length N.')
    end
else
    QPCONDENSINGSTEPS2 = [];
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

%% MPCexport (solver 1)

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

mpc1 = acado.OCPexport( ocp );

mpc1.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'       );
mpc1.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING'  );

if strcmp(ACADOSOLVER1,'qpOASES_N3')

    mpc1.set( 'QP_SOLVER',               'QP_QPOASES'         );
    mpc1.set( 'SPARSE_QP_SOLUTION',      'CONDENSING'         );

elseif strcmp(ACADOSOLVER1,'qpOASES_N2')

    mpc1.set( 'QP_SOLVER',               'QP_QPOASES'         );
    mpc1.set( 'SPARSE_QP_SOLUTION',      'FULL_CONDENSING_N2' );

elseif contains(ACADOSOLVER1,'qpDUNES') && QPCONDENSINGSTEPS <= 1

    mpc1.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc1.set( 'SPARSE_QP_SOLUTION',      'SPARSE_SOLVER'      );

elseif contains(ACADOSOLVER1,'qpDUNES') && QPCONDENSINGSTEPS > 1

    mpc1.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc1.set( 'SPARSE_QP_SOLUTION',      'BLOCK_CONDENSING_N2');
    mpc1.set( 'CONDENSING_BLOCK_SIZE',    QPCONDENSINGSTEPS   );

elseif strcmp(ACADOSOLVER1,'FORCES_BC') % TODO FIX

    if QPCONDENSINGSTEPS == 1
        mpc1.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER');
    else
        mpc1.set( 'SPARSE_QP_SOLUTION',  'BLOCK_CONDENSING_N2');
        mpc1.set( 'CONDENSING_BLOCK_SIZE', QPCONDENSINGSTEPS  );
    end

    mpc1.set( 'QP_SOLVER',               'QP_FORCES'          );


elseif strcmp(ACADOSOLVER1,'FORCES')

    mpc1.set( 'QP_SOLVER',           'QP_FORCES'              );
    mpc1.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER'          );

elseif strcmp(ACADOSOLVER1,'HPMPC')

    mpc1.set( 'QP_SOLVER',      'QP_HPMPC'                    );
    mpc1.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER'          );

else

    error('SPECIFIED SOLVER DOES NOT EXIST')

end

mpc1.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL2'        );
mpc1.set( 'NUM_INTEGRATOR_STEPS',        2*N                  );

mpc1.set( 'MAX_NUM_QP_ITERATIONS', 2000 );
mpc1.set( 'PRINTLEVEL', 'LOW' );

if MPC_EXPORT
    mpc1.exportCode( 'export_MPC1' );
end

if MPC_COMPILE
    cd export_MPC1

    % add correct version of qpDUNES (clipping or with qpOASES)
    if contains(ACADOSOLVER1,'qpDUNES')
        try
            rmdir('qpdunes','s');
        catch
            % if failed, directory not there anyway
        end
        if QPCONDENSINGSTEPS == 0
            % use clipping
            waitfor(copyfile('qpdunes_pub','qpdunes'));
        else
            % use qpDUNES+qpOASES
            waitfor(copyfile('qpdunes_dev','qpdunes'));
        end
    end
    if contains(ACADOSOLVER1, 'FORCES')
        % needed to give time to overwrite things
        keyboard
    end

    make_acado_solver('../acado_MPCstep1')
    cd ..
end

%% MPCexport (solver 2)

% acadoSet('problemname', 'mpc');
% 
% ocp = acado.OCP( 0.0, N*Ts, N );
% 
% rf = [xEnd; x; v; u];
% S  = acado.BMatrix(eye(NX+NU));
% 
% ocp.minimizeLSQ( S, rf );
% 
% rfN = [xEnd; x; v];
% SN  = acado.BMatrix(eye(NX));
% 
% ocp.minimizeLSQEndTerm( SN, rfN );
% 
% ocp.subjectTo( WALL <= [x([2:3:end]); xEnd(2)] <= 100 ); % constraint on y-position TODO upper bound 10?
% ocp.subjectTo( -1 <= u <= 1 );                          % box constraints on controls
% 
% ocp.setLinearInput(A1,B1);
% ocp.setModel(ode);

mpc2 = acado.OCPexport( ocp );

mpc2.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'       );
mpc2.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING'  );

if strcmp(ACADOSOLVER2,'qpOASES_N3')

    mpc2.set( 'QP_SOLVER',               'QP_QPOASES'         );
    mpc2.set( 'SPARSE_QP_SOLUTION',      'CONDENSING'         );

elseif strcmp(ACADOSOLVER2,'qpOASES_N2')

    mpc2.set( 'QP_SOLVER',               'QP_QPOASES'         );
    mpc2.set( 'SPARSE_QP_SOLUTION',      'FULL_CONDENSING_N2' );

elseif contains(ACADOSOLVER2,'qpDUNES') && QPCONDENSINGSTEPS <= 1

    mpc2.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc2.set( 'SPARSE_QP_SOLUTION',      'SPARSE_SOLVER'      );

elseif contains(ACADOSOLVER2,'qpDUNES') && QPCONDENSINGSTEPS > 1

    mpc2.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc2.set( 'SPARSE_QP_SOLUTION',      'BLOCK_CONDENSING_N2');
    mpc2.set( 'CONDENSING_BLOCK_SIZE',    QPCONDENSINGSTEPS   );

elseif strcmp(ACADOSOLVER2,'FORCES_BC') % TODO FIX

    if QPCONDENSINGSTEPS == 1
        mpc2.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER');
    else
        mpc2.set( 'SPARSE_QP_SOLUTION',  'BLOCK_CONDENSING_N2');
        mpc2.set( 'CONDENSING_BLOCK_SIZE', QPCONDENSINGSTEPS  );
    end

    mpc2.set( 'QP_SOLVER',               'QP_FORCES'          );


elseif strcmp(ACADOSOLVER2,'FORCES')

    mpc2.set( 'QP_SOLVER',           'QP_FORCES'              );
    mpc2.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER'          );

elseif strcmp(ACADOSOLVER2,'HPMPC')

    mpc2.set( 'QP_SOLVER',      'QP_HPMPC'                    );
    mpc2.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER'          );

else

    error('SPECIFIED SOLVER DOES NOT EXIST')

end

mpc2.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL2'        );
mpc2.set( 'NUM_INTEGRATOR_STEPS',        2*N                  );

mpc2.set( 'MAX_NUM_QP_ITERATIONS', 2000 );
mpc2.set( 'PRINTLEVEL', 'LOW' );

if MPC_EXPORT
    mpc2.exportCode( 'export_MPC2' );
end

if MPC_COMPILE
    cd export_MPC2

    % add correct version of qpDUNES (clipping or with qpOASES)
    if contains(ACADOSOLVER2,'qpDUNES')
        try
            rmdir('qpdunes','s');
        catch
            % if failed, directory not there anyway
        end
        if QPCONDENSINGSTEPS == 0
            % use clipping
            waitfor(copyfile('qpdunes_pub','qpdunes'));
        else
            % use qpDUNES+qpOASES
            waitfor(copyfile('qpdunes_dev','qpdunes'));
        end
    end
    if contains(ACADOSOLVER2, 'FORCES')
        % needed to give time to overwrite things
        keyboard
    end

    make_acado_solver('../acado_MPCstep2')
    cd ..
end

%% SIMULATIONS

minACADOtLog = []; % minimum timings of solver over NRUNS

for iRUNS = 1:NRUNS

    eval(['ref = textread(' '''' 'chain_mass_model_eq_M' num2str(NMASS) '.txt' '''', ', ''''' ');']);
    ref  = [ref(end-3+1:end); ref(1:end-3)]; % fix ordering convention

    if SCENARIO == 1
        eval(['X0  = textread(' '''' 'chain_mass_model_dist_M' num2str(NMASS) '.txt' '''', ', ''''' ');']);
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
        % run first solver
        output1 = acado_MPCstep1(input);
        % run second solver
        output2 = acado_MPCstep2(input);
        
        err_x = max(max(abs(output2.x - output1.x)));
        err_u = max(max(abs(output2.u - output1.u)));    

        if output1.info.status ~= 1 &&  output1.info.status ~= 0
            warning('SOLVER 1 FAILED')
            output1.info.cpuTime = NaN;
            keyboard
        end
        if output2.info.status ~= 1 &&  output2.info.status ~= 0
            warning('SOLVER 2 FAILED')
            output2.info.cpuTime = NaN;
            keyboard
        end

        ACADOnIter = [ACADOnIter output1.info.nIterations];
        ACADOoutputs{end+1} = output1;
        ACADOtLog = [ACADOtLog; output1.info.cpuTime];
        
        % Save the MPC step
        controls_MPC = [controls_MPC; output2.u(1,:)];

        % Shift trajectories
        input.x = [output2.x(2:end,:); output2.x(end,:)];
        input.u = [output2.u(2:end,:); output2.u(end,:)];

        % Simulate system
        sim_input.x = state_sim(end,:).';
        sim_input.u = output2.u(1,:).';

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

        disp(['current time: ' num2str(nextTime) '   ' char(9) ' (error in solutions: ' num2str(max(err_x,err_u)) ')'])

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
logging.solver  = ACADOSOLVER1;
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

end
