
% paper plots

clear all; clc

addpath([pwd filesep '../utils']);

%% plot condensing and sparse solvers
 
close all
clear variables

SAVEFIGS = 0;

NM = 5;  % number of masses 

xlims  = [10 80];

if NM == 3
    ylims1 = [0 60];
    ylims2 = [0 60];
elseif NM == 4
    ylims1 = [0 80];
    ylims2 = [0 80];    
elseif NM == 5
    ylims1 = [0 120];
    ylims2 = [0 120];
end

% first plot 

load(['M' num2str(NM) '_most_solvers.mat'], 'logs')

plot_1 = logs;

for ii = length(plot_1):-1:1
    
    if ~contains(plot_1{ii}.solver, 'qpOASES') %&& ~strcmp(plot_1{ii}.solver, 'FORCES')
        plot_1(ii) = [];
    end
    
end

plot_logs(plot_1, [], [], xlims, ylims1);

if SAVEFIGS
    exportfig(['~/Documents/Repositories/GIT/QP_story/paper/figures/solvers_1_M' num2str(NM) '.pdf'])
end

% second plot 

f = plot_logs(plot_1, true, [], xlims, ylims2);

% plot_2 = [];
load(['M' num2str(NM) '_most_solvers.mat'], 'logs')
plot_2 = logs;
% plot_2 = [plot_2 logs];
% load(['M' num2str(NM) '_HPMPC_B' num2str(NB) '.mat'], 'logs')
% plot_2 = [plot_2 logs];

for ii = length(plot_2):-1:1
    
    if contains(plot_2{ii}.solver, 'qpOASES') %|| strcmp(plot_2{ii}.solver, 'FORCES')
        plot_2(ii) = [];
    end
    
end


plot_logs(plot_2, false, f, xlims, ylims2);

if SAVEFIGS
    exportfig(['~/Documents/Repositories/GIT/QP_story/paper/figures/solvers_2_M' num2str(NM) '.pdf'])
end


% third plot

ff = plot_logs([plot_1 plot_2], true, [], xlims, ylims2);

load(['M' num2str(NM) '_bc_solvers.mat'], 'logs')
plot_3 = logs;
plot_logs(plot_3, false, ff, xlims, ylims2);

if SAVEFIGS
   exportfig(['~/Documents/Repositories/GIT/QP_story/paper/figures/solvers_3_M' num2str(NM) '.pdf'])
end

%% plot partial condensing
 
close all
clear variables

SAVEFIGS = 1;

load([pwd filesep 'M345_HPMPC_BC'], 'logs')
plot_partial_condensing(logs)

if SAVEFIGS
    exportfig(['~/Documents/Repositories/GIT/QP_story/paper/figures/bc_hpmpc.pdf'])
end

load([pwd filesep 'M345_QPDUNES_BC'], 'logs')
plot_partial_condensing(logs)

if SAVEFIGS
    exportfig(['~/Documents/Repositories/GIT/QP_story/paper/figures/bc_qpdunes.pdf'])
end



