
%% Initialize

clear all; close all; clc

USE_ACADO_DEV = 1; % leave to 1 if HPMPC is used in the benchmark

addpath([pwd filesep 'utils'])

% remove acado installations from path and add version of this repository
remove_acado_from_path();

curr_path = pwd;

if USE_ACADO_DEV
    cd(['..' filesep 'external' filesep 'acado-dev' filesep 'interfaces' filesep 'matlab']);
else
    cd(['..' filesep 'external' filesep 'acado' filesep 'interfaces' filesep 'matlab']);
end
make

% compile blasfeo and hpmpc
cd([curr_path filesep '..' filesep 'external' filesep 'blasfeo'])
system('make static_library')
cd([curr_path filesep '..' filesep 'external' filesep 'hpmpc'])
system('make static_library BLASFEO_PATH=$(pwd)/../blasfeo/')

cd(curr_path)

% add helper functions missing from older matlab versions
if verLessThan('matlab', 'R2016a')
    addpath([pwd filesep 'legacy'])
end

%% Choose simulation options

% AVAILABLE SOLVERS:
%
% 'qpDUNES'
% 'HPMPC'

SOLVER = 'HPMPC';

sim_opts.NRUNS       = 5;
sim_opts.MPC_EXPORT  = 1;
sim_opts.MPC_COMPILE = 1;
sim_opts.N           = 80;
sim_opts.CHECK_AGAINST_REF_SOL = 1;
sim_opts.SOL_TOL = 1e-6;

set_of_M  = [0 2:14 15:5:40];
set_of_NM = [3 4 5];

if strcmp(SOLVER, 'qpDUNES')
    set_of_M((mod(sim_opts.N, set_of_M) ~= 0) & set_of_M~=0 ) = [];
end

%% Run simulations

logs = cell(length(set_of_NM), length(set_of_M));

% delete code generated controller
[~] = rmdir('export_MPC', 's');
delete('acado_MPCstep.*');

for ii = 1:length(set_of_NM)
    
    sim_opts.VISUAL      = 1;
    sim_opts.SIM_EXPORT  = 1;
    sim_opts.SIM_COMPILE = 1;
    sim_opts.NMASS = set_of_NM(ii);

    for jj = 1:length(set_of_M)

        sim_opts.ACADOSOLVER = [SOLVER '_B' num2str(set_of_M(jj))];

        logs{ii, jj} = NMPC_chain_mass(sim_opts);

        % plot only once per NM
        sim_opts.VISUAL = 0;
        % compile sim only once per NM
        sim_opts.SIM_EXPORT  = 0;
        sim_opts.SIM_COMPILE = 0;
    
        % delete code generated controller
        [~] = rmdir('export_MPC', 's');
        delete('acado_MPCstep.*');
    end
    
    % clear all AFTER EACH SOLVER
    save('temp_data','sim_opts','logs','jj', 'set_of_M', 'set_of_NM', 'SOLVER');
    clear all; %#ok<CLALL>
    load('temp_data');
    
end

set_of_M(set_of_M == 0) = 1;

CPUTIMES = zeros(length(set_of_NM), length(set_of_M));
for ii = 1:size(CPUTIMES, 1)
    for jj = 1:size(CPUTIMES, 2)
        CPUTIMES(ii, jj) = max(logs{ii, jj}.cputime);
    end
end

plot(repmat(set_of_M, length(set_of_NM), 1)', CPUTIMES', '-o', 'linewidth', 1.5)
legend('HPMPC 3 masses', 'HPMPC 4 masses', 'HPMPC 5 masses', 'qpDUNES 3 masses', 'qpDUNES 4 masses', 'qpDUNES 5 masses')
title('qpDUNES and HPMPC with partial condensing');
ylabel('cpu time');
xlabel('block size');

SPEEDUPS = 1./(CPUTIMES./repmat(CPUTIMES(:,1), 1, length(set_of_M)));
plot(repmat(set_of_M, length(set_of_NM), 1)', SPEEDUPS', '-o', 'linewidth', 1.5)
legend('3 masses', '4 masses', '5 masses')
title(['speedup of partial condensing for ' SOLVER])
ylabel('speedup')
xlabel('block size')

%% Save results and clean up

% build log name based on current time and save results
t = clock;
t = t(1:end-1);     % remove seconds
t = mat2str(t);     % convert to string
t = t(2:end-1);     % remove [ ]
t(t == ' ') = '_';  % substitute spaces with underscore

save(['logs' filesep 'data_' t], 'logs');

if sim_opts.CHECK_AGAINST_REF_SOL
    max_val_err = 0;
    for i = 1:length(set_of_NM)
        for j = 1:length(logs)
            err = max(logs{i,j}.val_accuracy);
            if err > max_val_err
                max_val_err = err;
            end
        end
    end
    display(['max val err = ', num2str(max_val_err)])

    max_sol_err = 0;
    for i = 1:length(set_of_NM)
        for j = 1:length(logs)
            err = max(logs{i,j}.sol_accuracy);
            if err > max_sol_err
                max_sol_err = err;
            end
        end
    end
    display(['max sol err = ', num2str(max_sol_err)])
end

delete_temp_files();
