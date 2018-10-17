function [FHANDLE, maxerror] = plot_secondary_logs(logs, MODE, LOGSCALE, FHANDLE, xlims, ylims)

% PLOT_LOGS     plot performance of QP solvers as a function of prediction
%               horizon N.
%
% INPUTS:
%
% logs          logged data from simulation (cell array)
% MODE          'max' or 'av'
% FHANDLE       pass existing figure handle to get multiple logs in on plot

if nargin < 6
    ylims = [0 130];
end

if nargin < 5
    xlims = [10 100];
end


%% default values for inputs

if nargin < 4 || isempty(FHANDLE)
    FHANDLE = figure;
end

if nargin < 3 || isempty(LOGSCALE)
    LOGSCALE = false;
end

if nargin < 2 || isempty(MODE)
    MODE = 'max';
end

alpha_faded = 0.3;
style_faded = '--';
color_faded = [0 0 0 alpha_faded];

alpha = 1.0;
style = '-';
color = [0 0 0 alpha];

%% process data

solver = 'undefined';
nexp   = length(logs);
kk     = 0; % will contain number of different solvers in log cell array

for ii = 1:nexp
    
    if ~strcmp(solver,logs{ii}.secondary_solver)
        kk = kk+1;
        solver = logs{ii}.secondary_solver;
        data(kk).x = [];
        data(kk).y = [];
        data(kk).ya = [];
        data(kk).error = [];
        data(kk).marker = set_up_marker(solver);
        data(kk).solver = set_up_solver_name(solver);
    end
    
    try
    data(kk).x = [data(kk).x logs{ii}.N];
    catch
        keyboard
    end
    y  = logs{ii}.secondary_qptime; y(isnan(y)) = [];
    ya = logs{ii}.cputime - logs{ii}.simtime; ya(isnan(ya)) = [];
    
    data(kk).y = [data(kk).y y];
    data(kk).ya = [data(kk).ya ya];
    
    data(kk).error = [data(kk).error max(logs{ii}.secondary_error_sol)];

end

for kk = 1:length(data)
    if strcmp(MODE, 'max')
        data(kk).yplot  = max(data(kk).y);
        data(kk).ayplot = max(data(kk).ya);
    elseif strcmp(MODE, 'av')
        data(kk).ayplot = sum(data(kk).ya)/size(data(kk).ya, 1);
        data(kk).yplot  = sum(data(kk).y)/size(data(kk).y, 1);        
    end
end

%% plot error

if 0
    
    figure
    
    for kk = 1:length(data)
        plot(data(kk).x, data(kk).error);
        hold on
    end
    
end

%% plot timings

figure(FHANDLE);

if ~LOGSCALE
    
    for kk = 1:length(data)
        plot(data(kk).x, 1e3*data(kk).yplot, ...
            'Marker', data(kk).marker, 'MarkerSize', 12, 'MarkerEdgeColor', [1-alpha 1-alpha 1-alpha], ...
            'Color', color, 'Linewidth',1.5, 'LineStyle', style);
        hold on
    end
    plot(data(kk).x, 1e3*data(kk).ayplot, ...
        'Marker', set_up_marker('qpOASES_N2'), 'MarkerSize', 12, 'MarkerEdgeColor', [1-alpha_faded 1-alpha_faded 1-alpha_faded], ...
        'Color', color_faded, 'Linewidth',1.5, 'LineStyle', style_faded);
    
    grid on
    
    set_up_plot(data, false);
    xlim(xlims)
    ylim(ylims)

else
        
    for kk = 1:length(data)
        loglog(data(kk).x, 1e3*data(kk).yplot, ...
            'Marker', data(kk).marker, 'MarkerSize', 12, 'MarkerEdgeColor', [1-alpha 1-alpha 1-alpha], ...
            'Color', color, 'linewidth',1.5, 'LineStyle', style);
        hold on
    end
    loglog(data(kk).x, 1e3*data(kk).ayplot, ...
        'Marker', set_up_marker('qpOASES_N2'), 'MarkerSize', 12, 'MarkerEdgeColor', [1-alpha_faded 1-alpha_faded 1-alpha_faded], ...
        'Color', color_faded, 'Linewidth',1.5, 'LineStyle', style_faded);
    
    grid on
    
    set_up_plot(data, true);
    xlim(xlims)
    ylim(ylims)
end


FHANDLE.Position = [100 300 600 500];

maxerror = nan(length(data),1);
for kk = 1:length(data)
    maxerror(kk) = max(data(kk).error);
end

end

function set_up_plot(data, LOGPLOT)

set(gca, 'fontsize',20);
xlabel('Prediction horizon $N$', 'interpreter','latex', 'fontsize',20);
ylabel('CPU time $(\mathrm{ms})$', 'interpreter','latex', 'fontsize',20);

set(gca,'TickLabelInterpreter','latex')

if ~LOGPLOT
    hLegend = findobj(gcf, 'Type', 'Legend');
    
    if isempty(hLegend)
        l = legend({data.solver 'qpOASES CN$^2$'});
        l.Interpreter = 'latex';
        l.Location = 'northwest';
    end
end

if data(1).x(end) > data(1).x(1)
    xlim([data(1).x(1) data(1).x(end)])
end

end


function marker = set_up_marker(solver)

if strcmp(solver, 'qpOASES_N2')
    
    marker = '^';
    
elseif strcmp(solver, 'osqp')
    
    marker = 's';
    
elseif strcmp(solver, 'fiordos')
    
    marker = '*';
    
elseif strcmp(solver, 'dfgm')
    
    marker = 'p';
    
else 
    marker = 'o';
end

end


function solver_name_latex = set_up_solver_name(solver)

solver_name_latex = solver;

if strcmp(solver_name_latex, 'dfgm')
    solver_name_latex = 'DFGM';
end
if strcmp(solver_name_latex, 'osqp')
    solver_name_latex = 'OSQP';
end

end