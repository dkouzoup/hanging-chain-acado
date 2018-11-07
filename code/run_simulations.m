
%% Initialize

clear variables;
clear global;
close all;
clc

USE_ACADO_DEV = 1; % leave to 1 if HPMPC is used in the benchmark or to measure timings

initialize(USE_ACADO_DEV)

%% Choose simulation options

% AVAILABLE SOLVERS:
%
% 'qpOASES_N3'  qpOASES with N3 condensing
% 'qpOASES_N2'  qpOASES with N2 condensing
% 'qpDUNES_B0'  qpDUNES with clipping
% 'qpDUNES_BX'  qpDUNES with qpOASES and partial condensing with block size X
% 'HPMPC_B0'    HPMPC without partial condensing
% 'HPMPC_BX'    HPMPC with (its own) partial condensing and block size X
% 'FORCES'      FORCES QP solver (if license is available)

% set_of_solvers = {'qpOASES_N3', 'qpOASES_N2', 'qpDUNES_B0', 'HPMPC_B0', 'qpDUNES_B10', 'HPMPC_B10'};
% set_of_N       = 10:10:100;

% set_of_solvers = {'dfgm', 'osqp', 'qpOASES_N2', 'HPMPC_B0'}; 
set_of_solvers = {'dfgm', 'osqp'}; 
set_of_N       = 20:20:80; % fiordos crashses for N > 50 (and NMASS = 3)

% set_of_solvers = {'acados_qpOASES_e_N2', 'acados_HPIPM_B10', 'acados_HPIPM_B10', 'qpOASES_e_N2', 'HPMPC_B10'}; 
% set_of_N       = 10:10:80; 

sim_opts.WARMSTART   = 0;
sim_opts.NMASS       = 3;
sim_opts.NRUNS       = 8;
sim_opts.MPC_EXPORT  = 1;
sim_opts.MPC_COMPILE = 1;
sim_opts.SIM_EXPORT  = 1;
sim_opts.SIM_COMPILE = 1;
sim_opts.CHECK_AGAINST_REF_SOL = 1;

sim_opts.ACADOS_BENCHMARK = false;
for solver = set_of_solvers
   if contains(solver, 'acados')
       sim_opts.ACADOS_BENCHMARK = true;
   end
end

%% Run simulations

logs = {};

for jj = 1:length(set_of_solvers)

    sim_opts.VISUAL = 1;
    sim_opts.SOLVER = set_of_solvers{jj};

    for ii = 1:length(set_of_N)

        clear mex % for DFGM

        if set_of_N(ii) > 50 && strcmp(set_of_solvers{jj}, 'fiordos')
            break
        end
        sim_opts.N = set_of_N(ii);
        logs{end+1} = NMPC_chain_mass(sim_opts);

        % plot only once per solver
        sim_opts.VISUAL = 0;

    end

    % clear all AFTER EACH SOLVER
    save('temp_data','sim_opts','logs','ii', 'jj','set_of_solvers', 'set_of_N');
    clear all; %#ok<CLALL>
    load('temp_data');

end

%% Sanity checks

for ii = 1:length(logs)

    if ~check_consistency(logs{ii}.iters)
        warning(['inconsistent iterations between runs for log{' num2str(ii) '} with ' logs{ii}.solver])
        keyboard
    end

end

%% Save results and clean up

% build log name based on current time and save results
t = clock;
t = t(1:end-1);     % remove seconds
t = mat2str(t);     % convert to string
t = t(2:end-1);     % remove [ ]
t(t == ' ') = '_';  % substitute spaces with underscore

save(['logs' filesep 'data_' t], 'logs');

delete_temp_files();
