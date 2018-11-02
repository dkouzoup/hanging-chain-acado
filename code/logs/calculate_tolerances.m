function calculate_tolerances(logs)%, FADED, MODE, LOGSCALE, FHANDLE, xlims, ylims)

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
    ylims = [];
end

if nargin < 6
    xlims = [];
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
        data(kk).x  = [];
        data(kk).yp = [];
        data(kk).yd = [];
        data(kk).ye = [];
        data(kk).marker = set_up_marker(solver);
        data(kk).solver = set_up_solver_name(solver);
    end
    
    data(kk).x = [data(kk).x logs{ii}.N];
    if ~check_consistency(logs{ii}.primal_feas) || ~check_consistency(logs{ii}.dual_feas)
        warning(['log(' num2str(ii) ') with solver ' solver ' contains inconsistent data'])
    end
    
    data(kk).yp = [data(kk).yp max(logs{ii}.primal_feas(:,1))];
    data(kk).yd = [data(kk).yd max(logs{ii}.dual_feas(:,1))];
    data(kk).ye = [data(kk).ye max(logs{ii}.ref_sol_err(:,1))];
end

%% plot timings

% figure(FHANDLE);

for kk = 1:length(data)
    semilogy(data(kk).x, data(kk).yp, ...
        'Marker', data(kk).marker, 'MarkerSize', 12, 'MarkerEdgeColor', [1-alpha 1-alpha 1-alpha], ...
        'Color', color, 'Linewidth',1.5, 'LineStyle', '-');
    hold on
end
for kk = 1:length(data)
    semilogy(data(kk).x, data(kk).yd, ...
        'Marker', data(kk).marker, 'MarkerSize', 12, 'MarkerEdgeColor', [1-alpha 1-alpha 1-alpha], ...
        'Color', color, 'Linewidth',1.5, 'LineStyle', '--');
    hold on
end
grid on

set_up_plot(data, false);

if ~isempty(xlims)
    xlim(xlims)
end
if ~isempty(ylims)
    ylim(ylims)
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
        legend_txt = {};
        for ii = 1:length(data)
            legend_txt{ii} = [data(ii).solver ', primal residual'];
        end
        for ii = 1:length(data)
            legend_txt{ii+length(data)} = [data(ii).solver ', dual residual'];
        end
        
        l = legend(data.solver, legend_txt);
        l.Interpreter = 'latex';
        l.Location = 'northwest';
%     else
%         for ii = 1:length(data)
%             hLegend.String{end-ii+1} = data(end-ii+1).solver;  
%         end
    end
end

if data(1).x(end) > data(1).x(1)
    xlim([data(1).x(1) data(1).x(end)]);
end

% PRINT INFO

for kk = 1:length(data)
    disp('--------------------------------')
    disp(['SOLVER ' data(kk).solver])
    disp('--------------------------------')
    disp(['max primal res: ' num2str(max(data(kk).yp),'%2.2e')])
    disp(['max dual res:   ' num2str(max(data(kk).yd),'%2.2e')])
    disp(['max acado res:  ' num2str(max(data(kk).ye),'%2.2e')])
    disp('++++++++++++++++++++++++++++++++')
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