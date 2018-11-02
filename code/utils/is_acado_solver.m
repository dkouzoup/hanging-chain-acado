function answer = is_acado_solver(SOLVER)

answer = false;

if ~contains(SOLVER, 'acados')
    
    if contains(SOLVER, 'HPMPC') || contains(SOLVER, 'qpOASES') || ...
            contains(SOLVER, 'qpDUNES') || contains(SOLVER, 'FORCES')
        
        answer = true;
    end
    
end

end

