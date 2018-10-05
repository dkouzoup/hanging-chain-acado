function fiordos_code_generate(N, NX, NU, W, WN, opts)

%CODE_GENERATE_FIORDOS Code generate fiordos solver with time-varying
%                      dynamics and bounds (due to relative QP)

% remove any previously generated controller
if opts.export
    
    if isdir('fiordos_controller')
        if strfind(path, 'fiordos_controller') %#ok<STRIFCND>
            rmpath([pwd filesep 'fiordos_controller'])
        end
        rmdir('fiordos_controller', 's')
    end
    
    mkdir('fiordos_controller')
    
end

%% SET UP PROBLEM

% TODO: can we do it without infs?
ZZ = SimpleSet(2*N);

% input bounds
ZZ.addSet(1:N, EssBox(NU, 'l', 'param', 'u', 'param'));

% state bounds (x0 in equality constraints)
ZZ.addSet(N+(1:N),EssBox(NX, 'l', 'param','u', 'param'));

Q  =  W(1:NX,1:NX);
R  =  W(NX+1:NX+NU,NX+1:NX+NU);
QN =  WN;
H  = blkdiag(kron(eye(N), R),  kron(eye(N-1), Q), QN);

op = OptProb('H', H, 'g', 'param', 'X', ZZ, 'Ae', 'param', 'be', 'param', 'me', N*NX);

%% CODE GENERATE

cd('fiordos_controller')

if opts.export
    
    if strcmp(opts.approach, 'primal-dual')
        s = Solver(op, 'approach', opts.approach, 'algo', 'fgm');
        s.setSettings('approach', 'stopg', true, 'stopgEpsPrimal', opts.tol, 'stopgEpsDual', opts.tol, 'apprMaxit', opts.maxit);
    else
        s = Solver(op, 'approach', opts.approach, 'algoOuter', 'fgm'); % 'algoInner', 'fgm',
        s.setSettings('algoOuter', 'stopg', true, 'stopgEps', opts.tol, 'maxit', opts.maxit);
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

if opts.compile
    fiordos_mpc_mex_make;
end

cd('..');
addpath('fiordos_controller')

end

