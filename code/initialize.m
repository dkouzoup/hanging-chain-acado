function initialize(USE_ACADO_DEV)

if nargin < 1 || isempty(USE_ACADO_DEV)
    USE_ACADO_DEV = 1;
end

% set up matlab path properly and compile libs

persistent first_matlab_run

if isempty(first_matlab_run)
    
    first_matlab_run = 0;
    
    disp('... initializing')
        
    addpath([pwd filesep 'utils'])
    
    % remove acado installations from path and add version of this repository
    remove_acado_from_path();
    
    curr_path = pwd;
    
    if USE_ACADO_DEV
        cd(['..' filesep 'external' filesep 'acado-dev' filesep 'interfaces' filesep 'matlab']);
    else
        cd(['..' filesep 'external' filesep 'acado' filesep 'interfaces' filesep 'matlab']);
    end
    make
    
    % compile blasfeo and hpmpc
    cd([curr_path filesep '..' filesep 'external' filesep 'blasfeo'])
    system('make static_library')
    cd([curr_path filesep '..' filesep 'external' filesep 'hpmpc'])
    system('make static_library BLASFEO_PATH=$(pwd)/../blasfeo/')
    
    cd(curr_path)
    
    % add helper functions missing from older matlab versions
    if verLessThan('matlab', 'R2016a')
        addpath([pwd filesep 'legacy'])
    end
    
    addpath(genpath([pwd filesep '../external/acados_install']))
    
else
    disp('... skipping initialization')
end

end

