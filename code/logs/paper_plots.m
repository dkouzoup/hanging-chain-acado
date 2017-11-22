
% paper plots

clear all; clc

addpath([pwd filesep '../utils']);

%% plot 1a

close all

load('data_M4_most_solvers.mat', 'loggings')

plot_1_M4 = loggings;

for ii = length(plot_1_M4):-1:1
    
    if ~contains(plot_1_M4{ii}.solver, 'qpOASES')
        plot_1_M4(ii) = [];
    end
    
end

load('data_M4_FORCES.mat', 'loggings')

plot_1_M4 = [plot_1_M4 loggings];
% plot_1_M4 = [plot_1_M4(1:10) loggings plot_1_M4(11:end)];

plot_logs(plot_1_M4);

% save
exportfig('~/Documents/Repositories/GIT/QP_story/paper/figures/solvers_1_M4.pdf')

%% plot 1b

close all

load('data_M6_most_solvers.mat', 'loggings')

plot_1_M6 = loggings;

for ii = length(plot_1_M6):-1:1
    
    if ~contains(plot_1_M6{ii}.solver, 'qpOASES')
        plot_1_M6(ii) = [];
    end
    
end

load('data_M6_FORCES.mat', 'loggings')

plot_1_M6 = [plot_1_M6 loggings];

plot_logs(plot_1_M6);

% save
exportfig('~/Documents/Repositories/GIT/QP_story/paper/figures/solvers_1_M6.pdf')

%% plot 2a

close all

f = plot_logs(plot_1_M4, true);

plot_2_M4 = [];
load('data_M4_most_solvers.mat', 'loggings')
plot_2_M4 = [plot_2_M4 loggings];
load('data_M4_HPMPC_B8.mat', 'loggings')
plot_2_M4 = [plot_2_M4 loggings];

for ii = length(plot_2_M4):-1:1
    
    if contains(plot_2_M4{ii}.solver, 'qpOASES') || strcmp(plot_2_M4{ii}.solver, 'FORCES')
        plot_2_M4(ii) = [];
    end
    
end

plot_logs(plot_2_M4, false, f);

exportfig('~/Documents/Repositories/GIT/QP_story/paper/figures/solvers_2_M4.pdf')

%% plot 2b

close all

f = plot_logs(plot_1_M6, true);

plot_2_M6 = [];
load('data_M6_most_solvers.mat', 'loggings')
plot_2_M6 = [plot_2_M6 loggings];
load('data_M6_HPMPC_B10.mat', 'loggings')
plot_2_M6 = [plot_2_M6 loggings];

for ii = length(plot_2_M6):-1:1
    
    if contains(plot_2_M6{ii}.solver, 'qpOASES') || strcmp(plot_2_M6{ii}.solver, 'FORCES')
        plot_2_M6(ii) = [];
    end
    
end

plot_logs(plot_2_M6, false, f);

exportfig('~/Documents/Repositories/GIT/QP_story/paper/figures/solvers_2_M6.pdf')


