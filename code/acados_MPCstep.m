function output = acados_MPCstep( input )

% ACADOS_MPCSTEP Perform RTI step using the matlab interface of acados

%% perform RTI step

nlp = input.nlp;

nlp.set_field('lbx', 0, input.x0');
nlp.set_field('ubx', 0, input.x0');

% NOTE: provide initial guess at first call
if input.time == 0
    sol = nlp.solve(input.x0', zeros(3, 1));
else
    sol = nlp.solve();    
end

x_sol = sol.states();
u_sol = sol.controls();
info  = sol.info();

%% write to output struct

output.x   = horzcat(x_sol{:})';
output.u   = horzcat(u_sol{:})';

output.info.cpuTime = info.total_time;
output.info.nIterations = info.num_qp_iter;
output.info.objValue = nan;
output.info.status = info.status;

end

