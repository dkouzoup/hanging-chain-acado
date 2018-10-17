
clear all; close all; clc

PATH = '~/Documents/Repositories/GIT/thesis/image/qpstory/';

NM = 5;

load(['logs/dfgm_osqp_m' num2str(NM) '.mat']);

%%

SAVEFIGS = 1;

if NM ~= logs{1}.Nmass
    error('ups')
end

switch NM
   
    case 3
        ylim_max = 60;
        ylim_av  = 10;
    case 4
        ylim_max = 80;
        ylim_av  = 15;
    case 5
        ylim_max = 120;
        ylim_av  = 20;
end


[~, err] = plot_secondary_logs(logs, 'max', false, [], [10 80], [0 ylim_max])

if SAVEFIGS
    exportfig([PATH 'foms_max_M' num2str(NM) '.pdf'])
end

plot_secondary_logs(logs, 'max', true, [], [10 80], [0 ylim_max])

if SAVEFIGS
    exportfig([PATH 'foms_max_log_M' num2str(NM) '.pdf'])
end

plot_secondary_logs(logs, 'av', false, [], [10 80], [0 ylim_av])

if SAVEFIGS
    exportfig([PATH 'foms_av_M' num2str(NM) '.pdf'])
end

plot_secondary_logs(logs, 'av', true, [], [10 80], [0 ylim_av])

if SAVEFIGS
    exportfig([PATH 'foms_av_log_M' num2str(NM) '.pdf'])
end