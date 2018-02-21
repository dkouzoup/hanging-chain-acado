

clear all
clc
load('acados_qpdunes.mat');

for ii = 1:length(logs)
    
    fprintf('Error in iterations %d\t-  Error in solution %e\n', ...
        logs{ii}.acados_error_iter, logs{ii}.acados_error_sol)
    
end

plot_logs_with_acados(logs, [], false)