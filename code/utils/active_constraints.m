function [nbx, nbu] = active_constraints(output, WALL, UMAX)

% calculate number of active constraints in trajectory

TOL = 1e-4;
N   = size(output.x, 1) - 1;
NX  = size(output.x, 2); 
NMF = (NX-3)/6;

lbx = [repmat(repmat([-inf WALL -inf], N+1,1), 1, NMF+1) -inf(N+1, 3*NMF)];
ubu = UMAX*ones(N, 3);
lbu = -ubu;

nbx = sum(sum(abs(output.x - lbx) <= TOL));
nbu = sum(sum(abs(output.u - lbu) <= TOL | abs(output.u - ubu) <= TOL));

end

