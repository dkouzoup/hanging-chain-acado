function output = acados_MPCstep( input )

% ACADOS_MPCSTEP Perform RTI step using the matlab interface of acados

%% perform RTI step

nlp = input.nlp;

nlp.set_field('lbx', 0, input.x0');
nlp.set_field('ubx', 0, input.x0');

N = size(input.u, 1);

for ii = 1:N
    nlp.set_init_at_stage(input.x(ii,:)', input.u(ii,:)', ii-1); 
end
nlp.set_init_at_stage(input.x(N+1, :)', [], N); 

sol = nlp.solve();    

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

