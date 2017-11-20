function ode_rhs = chain_dynamics( x, v, xEnd, L, D, m, M, x0, f )

% takes an empty f of type acado.Expression in case an ACADO expression is 
% desired as output

g = [0; 0; -9.81];
f = f + repmat(g, M, 1);
Force = [];
for i = 1:M+1
    if i == 1
        dist = x((i-1)*3+1:i*3) - x0;
    elseif( i <= M )
        dist = x((i-1)*3+1:i*3) - x((i-2)*3+1:(i-1)*3);
    else
        dist = xEnd - x((M-1)*3+1:end);
    end

    scale = D/m*(1-L/norm(dist));
    F = scale*dist;

    Force = [Force; F];

    % mass on the right
    if i < M+1
        f((i-1)*3+1:i*3) = f((i-1)*3+1:i*3) - F;
    end
    % mass on the left
    if i > 1
        f((i-2)*3+1:(i-1)*3) = f((i-2)*3+1:(i-1)*3) + F;
    end
end

ode_rhs = [ v; ...
            f];

end

