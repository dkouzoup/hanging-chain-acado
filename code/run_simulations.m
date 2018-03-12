
%% Initialize

clear variables;
clear global;
close all;
clc

initialize()

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

% set_of_solvers = {'qpOASES_N3', 'qpOASES_N2', 'qpDUNES_B0', 'HPMPC_B0'};
% set_of_N       = 10:10:100;
set_of_solvers = {'qpOASES_e_N2'};
set_of_N       = 30; %10:10:100;

sim_opts.WARMSTART   = 0;
sim_opts.NMASS       = 4;
sim_opts.NRUNS       = 5;
sim_opts.MPC_EXPORT  = 1;
sim_opts.MPC_COMPILE = 1;
sim_opts.SIM_EXPORT  = 1;
sim_opts.SIM_COMPILE = 1;
sim_opts.CHECK_AGAINST_REF_SOL = 0;
sim_opts.CHECK_AGAINST_ACADOS  = 'qpoases';
sim_opts.SOL_TOL = 1e-3;

%% Run simulations

% % TEMP!
% if isempty(sim_opts.CHECK_AGAINST_ACADOS)
%     if sim_opts.NRUNS == 1
%         error('set more runs!')
%     end
% else
%    if sim_opts.NRUNS > 1
%         error('set one run!')
%    end
% end

logs = {};

for jj = 1:length(set_of_solvers)

    sim_opts.ACADOSOLVER = set_of_solvers{jj};

    sim_opts.VISUAL = 1;

    for ii = 1:length(set_of_N)

        sim_opts.N = set_of_N(ii);
        logs{end+1} = NMPC_chain_mass(sim_opts);

        % plot only once per solver
        sim_opts.VISUAL = 0;

    end

    % clear all AFTER EACH SOLVER
    save('temp_data','sim_opts','logs','jj','set_of_solvers','set_of_N');
    clear all; %#ok<CLALL>
    load('temp_data');

end

%% Save results and clean up

% build log name based on current time and save results
t = clock;
t = t(1:end-1);     % remove seconds
t = mat2str(t);     % convert to string
t = t(2:end-1);     % remove [ ]
t(t == ' ') = '_';  % substitute spaces with underscore

save(['logs' filesep 'data_' t], 'logs');

plot_logs(logs);

if sim_opts.CHECK_AGAINST_REF_SOL
    max_val_err = 0;
    for i = 1:length(logs)
        err = max(logs{i}.val_accuracy);
        if err > max_val_err
            max_val_err = err;
        end
    end
    display(['max err = ', num2str(max_val_err)])
end

delete_temp_files();
