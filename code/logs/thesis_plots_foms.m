
clear all; close all; clc

PATH = '~/Documents/Repositories/GIT/thesis/image/qpstory/';

NM  = 5;
TOL = 'low';

OSQP_WITH_SETUP = true

if OSQP_WITH_SETUP
    load(['logs/osqp_with_setup_M' num2str(NM) '_' TOL '.mat']);
    logs_osqp = logs;
end

load(['logs/foms_M' num2str(NM) '_' TOL '.mat']);

if OSQP_WITH_SETUP
    if ~strcmp(logs{1}.solver, 'osqp')
        error('assuming osqp is the first solver')
    end
    logs(1:length(logs_osqp)) = logs_osqp;
end

% TOL3 ERR
% 0.0239
% 0.0416

% TOL4 ERR
% 0.0020
% 0.0016
    
%%

% for ii = 1:length(logs)
%     if ~isfield(logs{ii}, 'secondary_solver')
%         logs{ii}.secondary_solver = 'none';
%     end
%     if strcmp(logs{ii}.secondary_solver, 'none')
%        logs{ii}.secondary_qptime = -1*ones(size(logs{ii}.cputime));
%        logs{ii}.secondary_error_sol = nan; 
%     end
% end

SAVEFIGS = 1;
    
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

% if SAVEFIGS
%     exportfig([PATH 'foms_max_M' num2str(NM) '_T' num2str(TOL) '.pdf'])
% end

% plot_secondary_logs(logs, 'max', true, [], [10 80], [0 ylim_max])
% 
% if SAVEFIGS
%     exportfig([PATH 'foms_max_log_M' num2str(NM) '.pdf'])
% end

% plot_secondary_logs(logs, 'av', false, [], [10 80], [0 ylim_av])
% 
% if SAVEFIGS
%     exportfig([PATH 'foms_av_M' num2str(NM) '_T' num2str(TOL) '.pdf'])
% end
% 
% % plot_secondary_logs(logs, 'av', true, [], [10 80], [0 ylim_av])
% % 
% % if SAVEFIGS
% %     exportfig([PATH 'foms_av_log_M' num2str(NM) '.pdf'])
% end