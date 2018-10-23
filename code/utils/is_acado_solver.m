function answer = is_acado_solver(SOLVER)

if contains(SOLVER, 'HPMPC') || contains(SOLVER, 'qpOASES') || ...
        contains(SOLVER, 'qpDUNES') || contains(SOLVER, 'FORCES')
    
    answer = true;
    
else
    
    answer = false;
    
end

end

