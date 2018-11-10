
clear all; close all; clc

PATH = '~/Documents/Repositories/GIT/thesis/image/qpstory/';
LA   = 'HP';
NM   = 5;

load(['acados_M' num2str(NM) '_' LA '.mat']);

% TODO: see why acado_ref fails for NM = 9 while acados not

%%

SAVEFIGS = 0;
    
if ~exist('NM', 'var')
    NM = logs{1}.Nmass;
else
    if NM ~= logs{1}.Nmass
        error('ups')
    end
end

switch NM
   
    case 3
        ylim_max = 60;
    case 4
        ylim_max = 80;
    case 5
        ylim_max = 120;
end

acados = zeros(length(logs),1);
acado  = zeros(length(logs),1);

for ii = 1:length(logs)
    if contains(logs{ii}.solver, 'acados')
        acados(ii) = 1;
    else
        acado(ii) = 1;
    end
end

xlim_max = 80;

% do NOT subtract min time from ACADO timings
FULL_RTI_TIMINGS = true;

fhandle = plot_timings(logs(find(acados)), false, 'max', false, [], [10 xlim_max], [0 ylim_max], FULL_RTI_TIMINGS);
plot_timings(logs(find(acado)), true, 'max', false, fhandle, [10 xlim_max], [0 ylim_max], FULL_RTI_TIMINGS);

if SAVEFIGS
    exportfig([PATH 'acados_M' num2str(NM) '_' LA '.pdf'])
end

fhandle = plot_timings(logs(find(acados)), false, 'max', true, [], [10 xlim_max], [0 ylim_max], FULL_RTI_TIMINGS);
plot_timings(logs(find(acado)), true, 'max', true, fhandle, [10 xlim_max], [0 ylim_max], FULL_RTI_TIMINGS);

if SAVEFIGS
    exportfig([PATH 'acados_M' num2str(NM)  '_' LA '_log.pdf'])
end

% calculate_tolerances(logs)
