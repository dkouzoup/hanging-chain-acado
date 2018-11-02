
clear all; close all; clc

PATH = '~/Documents/Repositories/GIT/thesis/image/qpstory/';

NM  = 5;

load(['acados_M' num2str(NM) '.mat']);

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
        ylim_av  = 20;
    case 4
        ylim_max = 80;
        ylim_av  = 15;
    case 5
        ylim_max = 120;
        ylim_av  = 20;
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

% do NOT subtract min time from ACADO timings
FULL_RTI_TIMINGS = true;

fhandle = plot_timings(logs(find(acados)), false, 'max', false, [], [10 80], [0 ylim_max], FULL_RTI_TIMINGS);
plot_timings(logs(find(acado)), true, 'max', false, fhandle, [10 80], [0 ylim_max], FULL_RTI_TIMINGS);

if SAVEFIGS
    exportfig([PATH 'acados_M' num2str(NM) '.pdf'])
end

fhandle = plot_timings(logs(find(acados)), false, 'max', true, [], [10 80], [0 ylim_max], FULL_RTI_TIMINGS);
plot_timings(logs(find(acado)), true, 'max', true, fhandle, [10 80], [0 ylim_max], FULL_RTI_TIMINGS);

if SAVEFIGS
    exportfig([PATH 'acados_M' num2str(NM) '_log.pdf'])
end

calculate_tolerances(logs)
