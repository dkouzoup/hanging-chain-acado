function [ FOUND ] = remove_acado_from_path(  )

% REMOVE_ACADO_FROM_PATH Remove acado from matlab path

FOUND = 0;

% temporarily disable warning
warning_status = warning('query', 'MATLAB:rmpath:DirNotFound');
warning('off', 'MATLAB:rmpath:DirNotFound');

while 1
    
    % check if acado is in path
    curr_path = path;
    acado_folders = strfind(curr_path, 'interfaces/matlab/acado');
    
    if isempty(acado_folders)
        break
    else 
        FOUND = 1;
    end
    
    separator = strfind(curr_path, ':');
    separator(separator > acado_folders(2)) = [];
    separator = separator(end);
    
    % form absolute path of acado
    acado_path = curr_path(separator+1:acado_folders(2)-1);
    
    % remove all acado subfolders from path
    rmpath(genpath(acado_path))
    
end

% restore warning to its status
warning(warning_status.state, 'MATLAB:rmpath:DirNotFound');

end

