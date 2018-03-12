function code_generate_acado_stats(sim_opts, qp_times, iters, fpath)

filename = sprintf('ocp_qp_data_nmasses_%d_nsteps_%d_solver_%s_warmstart_%d_ACADO_', ...
    sim_opts.NMASS, sim_opts.N, sim_opts.ACADOSOLVER, sim_opts.WARMSTART);

if fpath(end) ~= filesep
    fpath(end+1) = filesep;
end

cpu_times_files = fopen([fpath filename 'cpu_times.txt'], 'w');
print_double_vec(cpu_times_files, qp_times);
fclose(cpu_times_files);


iter_times_files = fopen([fpath filename 'iters.txt'], 'w');
print_int_vec(iter_times_files, iters);
fclose(iter_times_files);

end


function print_int_vec(fhandle, vec)

for ii = 1:length(vec)
    fprintf(fhandle, '%d\n', vec(ii));
end

end


function print_double_vec(fhandle, vec)

for ii = 1:length(vec)
    fprintf(fhandle, '%1.15e\n', vec(ii));
end

end

