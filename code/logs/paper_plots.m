
% paper plots
 
clear all; clc

addpath([pwd filesep '../utils']);

%% plot condensing and sparse solvers
 
close all
clear variables

SAVEFIGS = false; % this is for personal use, do not set to true
SAVEPATH = '~/Documents/Repositories/GIT/QP_story/paper/';

NM = 5;  % number of masses 

xlims  = [10 80];

if NM == 3
    ylims1 = [0 63];
    ylims2 = [0 63];
elseif NM == 4
    ylims1 = [0 80];
    ylims2 = [0 80];    
elseif NM == 5
    ylims1 = [0 120];
    ylims2 = [0 120];
end

% first plot 

load(['M' num2str(NM) '_all_solvers.mat'], 'logs')

plot_1 = logs;

for ii = length(plot_1):-1:1
    
    if ~contains(plot_1{ii}.solver, 'qpOASES')
        plot_1(ii) = [];
    end
    
end

plot_logs(plot_1, false, false, [], xlims, ylims1);

if SAVEFIGS
    exportfig([SAVEPATH 'solvers_1_M' num2str(NM) '.pdf'])
end

plot_logs(plot_1, false, true, [], xlims, ylims1);

if SAVEFIGS
    exportfig([SAVEPATH 'solvers_1_M' num2str(NM) '_log.pdf'])
end

% second plot 

plot_2 = logs;

for ii = length(plot_2):-1:1
    
    if contains(plot_2{ii}.solver, 'qpOASES') || contains(plot_2{ii}.solver, 'B10')
        plot_2(ii) = [];
    end
    
end

f1 = plot_logs(plot_1, true, false, [], xlims, ylims2);
plot_logs(plot_2, false, false, f1, xlims, ylims2);

if SAVEFIGS
    exportfig([SAVEPATH 'solvers_2_M' num2str(NM) '.pdf'])
end

f2 = plot_logs(plot_1, true, true, [], xlims, ylims2);
plot_logs(plot_2, false, true, f2, xlims, ylims2);

if SAVEFIGS
    exportfig([SAVEPATH 'solvers_2_M' num2str(NM) '_log.pdf'])
end

% third plot

plot_3 = logs;

for ii = length(plot_3):-1:1 
    if contains(plot_3{ii}.solver, 'qpOASES') || contains(plot_3{ii}.solver, 'B0')
        plot_3(ii) = [];
    end
end

ff1 = plot_logs([plot_1 plot_2], true, false, [], xlims, ylims2);
plot_logs(plot_3, false, false, ff1, xlims, ylims2);

if SAVEFIGS
   exportfig([SAVEPATH 'solvers_3_M' num2str(NM) '.pdf'])
end

ff2 = plot_logs([plot_1 plot_2], true, true, [], xlims, ylims2);
plot_logs(plot_3, false, true, ff2, xlims, ylims2);

if SAVEFIGS
   exportfig([SAVEPATH 'solvers_3_M' num2str(NM) '_log.pdf'])
end

%% plot partial condensing
 
close all
clear variables

SAVEFIGS = false; % this is for personal use, do not set to true
SAVEPATH = '~/Documents/Repositories/GIT/QP_story/paper/';

load([pwd filesep 'BC_HPMPC'], 'logs')
plot_partial_condensing(logs)

if SAVEFIGS
    exportfig([SAVEPATH 'bc_hpmpc.pdf'])
end

load([pwd filesep 'BC_qpDUNES'], 'logs')
plot_partial_condensing(logs)

if SAVEFIGS
    exportfig([SAVEPATH 'bc_qpdunes.pdf'])
end



