function plot_results( loggings, Tprep )

% PLOT_RESULTS plot performance of QP solvers as a function of prediction
%   horizon N.

if nargin == 1 || (nargin == 2 && isempty(Tprep))
    Tprep = 0;
end

%% process data

solver = 'undefined';
nexp   = length(loggings);
kk     = 0;

for ii = 1:nexp
    
    if ~strcmp(solver,loggings{ii}.solver)
        kk = kk+1;
        solver = loggings{ii}.solver;
        data(kk).s = solver;
        data(kk).x = [];
        data(kk).y = [];
        
        % update solver name for proper legends
        data(kk).s(data(kk).s == '_') = ' ';
        if contains(data(kk).s, 'qpOASES')
            data(kk).s = replace(data(kk).s, 'N', 'C$N^');
            data(kk).s(end+1) = '$';
        end
        if contains(data(kk).s, 'HPMPC B')
            data(kk).s = replace(data(kk).s, 'B', 'B$_{');
            data(kk).s(end+1:end+2) = '}$';
        end
        if strcmp(data(kk).s, 'qpDUNES B0')
            data(kk).s = 'qpDUNES';
        end
    end
    
    data(kk).x = [data(kk).x loggings{ii}.N];
    data(kk).y = [data(kk).y loggings{ii}.cputime];
  
end

%% plot average timings

figure

for kk = 1:length(data)
    plot(data(kk).x, 1e3*(sum(data(kk).y)/size(data(kk).y,1) - Tprep), '-o', 'linewidth',1.5);
    hold on
end

set_up_plot(data);
title('Average CPU time in closed-loop','interpreter','latex', 'fontsize',20)

%% plot worst-case timings

figure

for kk = 1:length(data)
    plot(data(kk).x, 1e3*(max(data(kk).y)-Tprep), '-o', 'linewidth',1.5);
    % loglog(data(kk).x, 1e3*max(data(kk).y), '-o', 'linewidth',1.5);
    hold on
end

set_up_plot(data);
title('Worst case CPU time in closed-loop','interpreter','latex', 'fontsize',20)

end

function set_up_plot(data)
    set(gca, 'fontsize',20);
    xlabel('Prediction horizon $N$', 'interpreter','latex', 'fontsize',20);
    ylabel('Time $(\mathrm{ms})$', 'interpreter','latex', 'fontsize',20);
    l = legend(data.s);
    l.Interpreter = 'latex';
    l.Location = 'northwest';
    if data(1).x(end) > data(1).x(1)
        xlim([data(1).x(1) data(1).x(end)])
    end
end
