function dfgm_compile(N, NX, NU, opts)

% DFGM_COMPILE Compile DFGM for given problem dimensions

tmp = pwd;

cd('../external/mexedDGM/');
updateProblemDimensions(N, NX, NU);

OPT.useExternalLibraries = 1;

if opts.criterion == 1
    OPT.terminationCondition = 1; % 1 for gradient map, 2 for maxit, 3 for known optimal solution    
elseif opts.criterion == 2
    OPT.terminationCondition = 3; %#ok<*STRNU>
end

make
  
txt = [tmp '/mexedDGM.mexmaci64'];

copyfile('gdfgm_mex.mexmaci64', txt);

cd(tmp)

end

