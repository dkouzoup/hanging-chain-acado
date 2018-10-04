
% Make function for Matlab


% seperator (MAC/WINDOWS)
sep    = pwd;
sep    = sep(1);

if strcmp(sep,'/') ~= 1
    warning('LAPACK and BLAS libraries in a windows PC must be added in a different way')
end

% Pre-processor commands
DFLAGS = '-DUSE_IN_MATLAB ';
LDFLAGS = [];

try
    if OPT.useExternalLibraries == 1
        DFLAGS = [DFLAGS ' -DUSE_EXTERNAL_LIBRARIES '];
        
        % Linker flags
        LDFLAGS = 'LDFLAGS="\$LDFLAGS -framework Accelerate"';
    end
catch
    %do nothing
end

try
    if OPT.terminationCondition == 0
        DFLAGS = [DFLAGS '-DSTOP_OPTIMAL_SOLUTION '];
    elseif OPT.terminationCondition == 1
        DFLAGS = [DFLAGS '-DSTOP_GRADIENT_MAP '];
    elseif OPT.terminationCondition == 2
        DFLAGS = [DFLAGS '-DSTOP_MAXIMUM_ITERATIONS '];
    end
catch
    % default termination condition
    DFLAGS = [DFLAGS '-DSTOP_OPTIMAL_SOLUTION ' ];
    warning('Termination condition not found in options, using gradient map.')
end

% C files
FILES   = [' gdfgm_mex.c ' 'gdfgm_solve.c ' 'gdfgm_utils.c ' 'gdfgm_algebra.c ' 'gdfgm_core.c ' 'timing.c '];

% mex command
eval(['mex -v ' DFLAGS LDFLAGS FILES]);



