clear all; clc;

% change HPMPC solver name to match the block size

M = 4;
B = 'N';

if isnumeric(B)
    Bstr = num2str(B);
else
    Bstr = B;
end

PATH = ['logs' filesep];

load([PATH 'loggings_M' num2str(M) '_HPMPC_B' Bstr '.mat'])

for ii = 1:length(loggings)
   
    if contains(loggings{ii}.solver, 'HPMPC')
        loggings{ii}.solver = ['HPMPC_B' Bstr];
    end
    
end

save([PATH 'loggings_M' num2str(M) '_HPMPC_B' Bstr '.mat'], 'loggings')
