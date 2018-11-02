function QPCONDENSINGSTEPS = extract_block_size(SOLVER)

% EXTRACT_BLOCK_SIZE extract block size from solver name. Solver name
%                    follows the convention SOLVER_BXX, with XX the size of
%                    the blocks in partial condensing

if contains(SOLVER,'qpDUNES')
    bpos = strfind(SOLVER,'B');
    QPCONDENSINGSTEPS = str2double(SOLVER(bpos+1:end));
    if QPCONDENSINGSTEPS ~= 0 &&(QPCONDENSINGSTEPS < 1 || mod(N, QPCONDENSINGSTEPS)~= 0)
        error('Invalid block size for given horizon length N and QP solver.')
    end
elseif contains(SOLVER,'HPMPC')
    bpos = strfind(SOLVER,'B');
    QPCONDENSINGSTEPS = str2double(SOLVER(bpos+1:end));
    if ~(QPCONDENSINGSTEPS >=0)
        error('Invalid block size');
    end
else
    QPCONDENSINGSTEPS = [];
end

end

