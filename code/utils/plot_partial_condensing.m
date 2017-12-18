function FHANDLE = plot_partial_condensing(logs)

%% process data

solver = logs{1}.solver(1:strfind(logs{1}.solver,'_')-1);
NMASS  = size(logs, 1);
BS     = size(logs,2);

CPUTIMES  = zeros(NMASS, BS);
BLOCKSIZE = zeros(NMASS, BS);

for ii = 1:NMASS
    for jj = 1:BS
        if ~contains(logs{ii, jj}.solver, solver)
            error('log of different solver detected')
        end
        CPUTIMES(ii, jj)  = max(logs{ii, jj}.cputime - logs{ii, jj}.simtime);
        BLOCKSIZE(ii, jj) = str2double(logs{ii, jj}.solver(strfind(logs{ii, jj}.solver, '_B')+2:end));
    end
end

SPEEDUPS = 1./(CPUTIMES./repmat(CPUTIMES(:,1), 1, BS));
BLOCKSIZE(BLOCKSIZE == 0) = 1;
 
%% plot

FHANDLE = figure;

legends = {};
for ii = 1:NMASS
    MARKER   = set_up_marker(logs{ii, 1}.Nmass);
    plot(BLOCKSIZE(ii,:), SPEEDUPS(ii,:), 'Marker', MARKER, 'MarkerSize', 12, 'color', 'k', 'Linewidth',1.5, 'LineStyle', '-');
    hold on
    legends{end+1} = ['$n_{\mathrm{m}} = ' num2str(logs{ii, 1}.Nmass) '$'];
end
grid on

set(gca, 'fontsize',20);
xlabel('Block size $M$', 'interpreter','latex', 'fontsize',20);
ylabel('Speedup', 'interpreter','latex', 'fontsize',20);
set(gca,'TickLabelInterpreter','latex')
title(['Partial condensing with \texttt{' solver '}'],'interpreter','latex', 'fontsize',20);
l = legend(legends);
l.Interpreter = 'latex';
l.Location = 'northeast';

xlim([BLOCKSIZE(1,1) BLOCKSIZE(1,end)]);

WIDE = 0;

if WIDE
    FHANDLE.Position = [100 300 1200 500];
else
    FHANDLE.Position = [100 300 600 500];    
end


end


function marker = set_up_marker(nmasses)

if nmasses == 3
    marker = 'o';
elseif nmasses == 4
    marker = '>';
elseif nmasses == 5
    marker = 'h';
else
    marker = '.';
end

end
