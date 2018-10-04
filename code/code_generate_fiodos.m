function code_generate_fiodos(N, NX, NU, FIORDOS_EXPORT, FIORDOS_COMPILE)

%CODE_GENERATE_FIORDOS Code generate fiordos solver with time-varying
% dynamics and bounds (due to relative QP)


% 'dual' or 'primal-dual'
APPROACH = 'dual';
tol      = 1e-3;
maxit    = 100;

% remove any previously generated controller
if FIORDOS_EXPORT
    
    if isdir('fiordos_controller')
        rmdir('fiordos_controller', 's')
    end
    
    mkdir('fiordos_controller')
    
end

M = (NX/3 - 1)/2;

%% SET UP PROBLEM

% TODO: CAN WE DO IT WITHOUT INFS?
ZZ = SimpleSet(2*N);

% input bounds
ZZ.addSet(1:N, EssBox(NU, 'l', 'param', 'u', 'param'));

% state bounds (x0 in equality constraints)
ZZ.addSet(N+(1:N),EssBox(NX, 'l', 'param','u', 'param'));

op = OptProb('H', 'param.diag', 'g', 'param', 'X', ZZ, 'Ae', 'param', 'be', 'param', 'me', N*NX);

%% CODE GENERATE

cd('fiordos_controller')

if FIORDOS_EXPORT
    
    % TODO: CLIPPING NOT POSSIBLE??
    if strcmp(APPROACH, 'primal-dual')
        s = Solver(op, 'approach', APPROACH, 'algo', 'fgm');
        s.setSettings('approach', 'stopg', true, 'stopgEpsPrimal', tol, 'stopgEpsDual', tol, 'apprMaxit', maxit);
    else
        s = Solver(op, 'approach', APPROACH, 'algoOuter', 'fgm'); % 'algoInner', 'fgm',
        s.setSettings('algoOuter', 'stopg', true, 'stopgEps', tol, 'maxit', maxit);
    end
    
    s.generateCode('prefix','fiordos_mpc_', 'forceOverwrite', true);
    % s.setAutoPrecond('algo'); % can't precondition due to changing data
    % s.setSettings('approach', 'inlineA',true); % not much improvement
    rehash;
    
    % process generated code to add timings
    copyfile('../fiordos_utils/timing.c','./')
    copyfile('../fiordos_utils/timing.h','./')
    
    fid_in  = fopen('fiordos_mpc_mex.c', 'r');
    fid_tmp = fopen('tmp','wt');
    
    if fid_in <= 0
        error('Error reading file');
    end
    
    % scan and edit fiordos_mpc_mex.c file
    while ~feof(fid_in)
        
        % read line
        fline = fgetl(fid_in);
        
        % add header
        check_header = strfind(fline, '#include "fiordos_mpc_solver.h"');
        
        % tic toc
        check_solve = strfind(fline, 'fiordos_mpc_solve(&params, &settings, &result, &work);');
        
        % augment output
        check_output1a = strfind(fline, 'const char *fnames[] = {"la", "x", "iter", "exitflag"}');
        check_output1b = strfind(fline, 'const char *fnames[] = {"la", "x", "d", "iter", "exitflag"}');
        check_output2a = strfind(fline, 'MRES = mxCreateStructMatrix(1,1, 4, fnames)');
        check_output2b = strfind(fline, 'MRES = mxCreateStructMatrix(1,1, 5, fnames)');
        check_output3 = strfind(fline, 'mxSetField(MRES, 0, "exitflag", tmp1);');
        
        % replace escape characters for printf
        fline = strrep(fline, '\', '\\');
        fline = strrep(fline, '%', '%%');
        
        if ~isempty(check_header) %#ok<*STREMP>
            fline = [fline '\n#include "timing.h"'];
        end
        if ~isempty(check_output1a)
            fline = 'const char *fnames[] = {"la", "x", "iter", "exitflag", "cputime"};';
        end
        if ~isempty(check_output1b)
            fline = 'const char *fnames[] = {"la", "x", "d", "iter", "exitflag", "cputime"};';
        end
        if ~isempty(check_output2a)
            fline = 'MRES = mxCreateStructMatrix(1,1, 5, fnames);';
        end
        if ~isempty(check_output2b)
            fline = 'MRES = mxCreateStructMatrix(1,1, 6, fnames);';
        end
        if ~isempty(check_output3)
            fline = [fline '\n\n\t\ttmp1 = mxCreateDoubleMatrix(1,1,mxREAL);\n\t\t*mxGetPr(tmp1) = (double) cputime;\n\t\tmxSetField(MRES, 0, "cputime", tmp1);'];
        end
        if ~isempty(check_solve)
            fline = ['\tacados_timer tmr;\n\tacados_tic(&tmr);\n' fline '\n\tdouble cputime = acados_toc(&tmr);'];
        end
        
        % add end of line character
        fnew  = [fline '\n'];
        
        % print in temp file
        fprintf(fid_tmp, fnew);
    end
    
    % close files
    fclose(fid_in);
    fclose(fid_tmp);
    
    % copy and delete temp
    copyfile('tmp', 'fiordos_mpc_mex.c');
    delete('tmp');
    
    % process makefile
    fid_in  = fopen('fiordos_mpc_mex_make.m', 'r');
    fid_tmp = fopen('tmp','wt');
    
    if fid_in <= 0
        error('Error reading file');
    end
    
    
    % scan and edit fiordos_mpc_mex.c file
    while ~feof(fid_in)
        
        % read line
        fline = fgetl(fid_in);
        
        % add header
        check = strfind(fline, 'fiordos_mpc_solver.c');
        
        % replace escape characters for printf
        fline = strrep(fline, '\', '\\');
        fline = strrep(fline, '%', '%%');
        
        if ~isempty(check) %#ok<*STREMP>
            fline = [fline ' timing.c'];
        end
        
        % add end of line character
        fnew  = [fline '\n'];
        
        % print in temp file
        fprintf(fid_tmp, fnew);
    end
    
    % close files
    fclose(fid_in);
    fclose(fid_tmp);
    
    % copy and delete temp
    copyfile('tmp', 'fiordos_mpc_mex_make.m');
    delete('tmp');
    
end

if FIORDOS_COMPILE
    fiordos_mpc_mex_make;
end

cd('..');
addpath('fiordos_controller')

end

