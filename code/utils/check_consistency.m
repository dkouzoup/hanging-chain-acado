function answer = check_consistency(M)

% CHECK_CONSISTENCY Check whether all columns of matrix M, each
%                   representing a different RUN are the equal to each
%                   other


answer = true;

NRUNS = size(M, 2);

for ii = 2:NRUNS
    
    if ~isequal(M(:,ii), M(:,1))
        answer = false;
        break;
    end
end

end

