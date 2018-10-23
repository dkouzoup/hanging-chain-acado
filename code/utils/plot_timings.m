function [FHANDLE] = plot_timings(logs, FADED, MODE, LOGSCALE, FHANDLE, xlims, ylims)

% PLOT_LOGS     plot performance of QP solvers as a function of prediction
%               horizon N.
%
% INPUTS:
%
% logs          logged data from simulation (cell array)
% FADED         set to true to plot solver curves faded (boolean)
% MODE          'max' or 'av' to plot worst-case or average timings
% FHANDLE       pass existing figure handle to get multiple logs in on plot
% xlims         optional limits to x axis
% ylims         optional limits to y axis


%% default values for inputs

if nargin < 7
    ylims = [0 100];
end

if nargin < 6
    xlims = [10 100];
end

if nargin < 5 || isempty(FHANDLE)
    FHANDLE = figure;
end

if nargin < 4 || isempty(LOGSCALE)
    LOGSCALE = false;
end

if nargin < 3 || isempty(MODE)
    MODE = 'max';
end

if nargin < 2 || isempty(FADED)
    FADED = false;
end

if FADED
    alpha = 0.3;
    style = '--';
else
    alpha = 1.0;
    style = '-';
end

color = [0 0 0 alpha];

%% process data

solver = 'undefined';
nexp   = length(logs);
kk     = 0; % will contain number of different solvers in log cell array

for ii = 1:nexp
    
    if ~strcmp(solver,logs{ii}.solver)
        kk = kk+1;
        solver = logs{ii}.solver;
        data(kk).x = [];
        data(kk).y = [];
        data(kk).marker = set_up_marker(solver);
        data(kk).solver = set_up_solver_name(solver);
    end
    
    data(kk).x = [data(kk).x logs{ii}.N];
    cputimings = logs{ii}.cputime - logs{ii}.simtime;
    data(kk).y = [data(kk).y min(cputimings, [], 2)];
  
end

for kk = 1:length(data)
    if strcmp(MODE, 'max')
        data(kk).yplot  = max(data(kk).y);
    elseif strcmp(MODE, 'av')
        nnans = sum(isnan(data(kk).y(:,1)));
        y_tmp = data(kk).y(nnans+1:end,:);
        data(kk).yplot  = sum(y_tmp)/size(y_tmp,1);
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
    grid on
    
    set_up_plot(data, true);
    xlim(xlims)
    ylim(ylims)
    % title('Worst case CPU time in closed-loop','interpreter','latex', 'fontsize',20)
end


FHANDLE.Position = [100 300 600 500];

end

function set_up_plot(data, LOGPLOT)

set(gca, 'fontsize',20);
xlabel('Prediction horizon $N$', 'interpreter','latex', 'fontsize',20);
ylabel('CPU time $(\mathrm{ms})$', 'interpreter','latex', 'fontsize',20);

set(gca,'TickLabelInterpreter','latex')

if ~LOGPLOT
    hLegend = findobj(gcf, 'Type', 'Legend');
    
    if isempty(hLegend)
        l = legend(data.solver);
        l.Interpreter = 'latex';
        l.Location = 'northwest';
    else
        for ii = 1:length(data)
            hLegend.String{end-ii+1} = data(end-ii+1).solver;  
        end
    end
end

if data(1).x(end) > data(1).x(1)
    xlim([data(1).x(1) data(1).x(end)])
end

end


function marker = set_up_marker(solver)

if strcmp(solver, 'qpOASES_N2')
    
    marker = '^';
    
elseif strcmp(solver, 'qpOASES_N3')
    
    marker = 'v';

elseif strcmp(solver, 'FORCES')
    
    marker = 's';
    
elseif strcmp(solver, 'qpDUNES') || strcmp(solver, 'qpDUNES_B0')
    
    marker = 'p';
    
elseif strcmp(solver, 'HPMPC') || strcmp(solver, 'HPMPC_B0')
    
    marker = '*';
    
elseif contains(solver, 'HPMPC_B')
    
    marker = 'x';
    
elseif contains(solver, 'qpDUNES_B')
    
    marker = 'd';
    
elseif strcmp(solver, 'osqp')
    
    marker = 'o';
    
elseif strcmp(solver, 'dfgm')
    
    marker = 'h';
    
elseif strcmp(solver, 'fiordos')
    
    marker = '>';
    
else 
    marker = '<';
end

end


function solver_name_latex = set_up_solver_name(solver)

solver_name_latex = solver;

solver_name_latex(solver_name_latex == '_') = ' ';
if contains(solver_name_latex, 'qpOASES')
    solver_name_latex = replace(solver_name_latex, 'N', 'C$N^');
    solver_name_latex(end+1) = '$';
end

if strcmp(solver_name_latex, 'HPMPC B0')
    solver_name_latex = 'HPMPC';
end

if strcmp(solver_name_latex, 'qpDUNES B0')
    solver_name_latex = 'qpDUNES';
end

if contains(solver_name_latex, 'HPMPC B') || contains(solver_name_latex, 'qpDUNES B')
    solver_name_latex = [solver_name_latex(1:strfind(solver_name_latex, 'B')-1) 'PC'];
    % solver_name_latex = replace(solver_name_latex, 'B', 'B$_{');
    % solver_name_latex(end+1:end+2) = '}$';
end

if strcmp(solver_name_latex, 'dfgm')
    solver_name_latex = 'DFGM';
end

if strcmp(solver_name_latex, 'osqp')
    solver_name_latex = 'OSQP';
end

end