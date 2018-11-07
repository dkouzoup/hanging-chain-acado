
%% load mat files

clear all; close all; clc

PATH = '~/Documents/Repositories/GIT/thesis/image/qpstory/';

NM  = 3;

OSQP_WITH_SETUP = true;

if ~OSQP_WITH_SETUP
    error('timings not available');
end

load(['logs/foms_M' num2str(NM) '.mat']);

% TOL3 ERR
% 0.0239
% 0.0416

% TOL4 ERR
% 0.0020
% 0.0016
    
%% plot results

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

foms = zeros(length(logs),1);
soms = zeros(length(logs),1);

for ii = 1:length(logs)
    if strcmp(logs{ii}.solver, 'fiordos') || strcmp(logs{ii}.solver, 'dfgm') || strcmp(logs{ii}.solver, 'osqp')
        foms(ii) = 1;
    else
        soms(ii) = 1;
    end
end

fhandle = plot_timings(logs(find(foms)), false, 'max', false, [], [10 80], [0 ylim_max]);
plot_timings(logs(find(soms)), true, 'max', false, fhandle, [10 80], [0 ylim_max]);

if SAVEFIGS
    exportfig([PATH 'foms_M' num2str(NM) '.pdf'])
end

fhandle = plot_timings(logs(find(foms)), false, 'max', true, [], [10 80], [0 ylim_max]);
plot_timings(logs(find(soms)), true, 'max', true, fhandle, [10 80], [0 ylim_max]);

if SAVEFIGS
    exportfig([PATH 'foms_M' num2str(NM) '_log.pdf'])
end

calculate_tolerances(logs)
