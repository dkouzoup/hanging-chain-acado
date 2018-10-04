function dfgm_compile(N, NX, NU)

%DFGM_COMPILE Compile DFGM for given problem dimensions

tmp = pwd;

cd('../external/mexedDGM/');
updateProblemDimensions(N, NX, NU);

OPT.useExternalLibraries = 1;    
OPT.terminationCondition = 1; % 1 for gradient map, 2 for maxit, 3 for known optimal solution    

make
  
% txt = [tmp '/mexedDGM' num2str(N) '.mexmaci64'];
txt = [tmp '/mexedDGM.mexmaci64'];

copyfile('gdfgm_mex.mexmaci64', txt);

cd(tmp)

end

