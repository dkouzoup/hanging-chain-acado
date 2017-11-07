function delete_temp_files( )

% DELETE_TEMP_FILES delete all code generated files

warning_status = warning('query', 'MATLAB:DELETE:FileNotFound');
warning('off', 'MATLAB:DELETE:FileNotFound')
delete('temp_data.mat')
delete('mpc_data*.txt')
delete('sim_data*.txt')
delete('sim.cpp')
delete('mpc.cpp')
delete('mpc_RUN.m')
delete('mpc_RUN.mex*')
delete('sim_RUN.mex*')
delete('sim_RUN.m')
warning(warning_status.state, 'MATLAB:DELETE:FileNotFound');

if exist('export_MPC', 'dir')
    rmdir('export_MPC', 's')
end
if exist('export_SIM', 'dir')
    rmdir('export_SIM', 's')
end

end

