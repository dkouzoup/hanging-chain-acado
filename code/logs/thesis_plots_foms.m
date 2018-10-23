
clear all; close all; clc

PATH = '~/Documents/Repositories/GIT/thesis/image/qpstory/';

NM  = 3;
TOL = 4;

load(['logs/dfgm_osqp_hpmpc_qpoases_m' num2str(NM) '_t' num2str(TOL) '.mat']);

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
        ylim_max = 45;
        ylim_av  = 20;
    case 4
        ylim_max = 80;
        ylim_av  = 15;
    case 5
        ylim_max = 120;
        ylim_av  = 20;
end


fhandle = plot_timings(logs(1:2), true, 'max', false, [], [10 80], [0 ylim_max])
plot_timings(logs(3:end), false, 'max', false, fhandle, [10 80], [0 ylim_max])

fhandle = plot_timings(logs(1:2), true, 'av', false, [], [10 80], [0 ylim_max])
plot_timings(logs(3:end), false, 'av', false, fhandle, [10 80], [0 ylim_max])


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