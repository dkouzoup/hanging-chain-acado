function [ Ai ] = build_inequality_constraints_sparsity(N, NX, NU)

% BUILD_INEQUALITY_CONSTRAINTS_SPARSITY construct sparsity of bounds 

M  = (NX/3-1)/2;
lx = -inf*ones(NX,1);
lx(2:3:3*(M+1)) = -1;
lx = repmat(lx, N, 1);
lu = -1*ones(N*NU,1);
lz = [lu; lx];
Ai = eye(N*(NX+NU));
Ai(isinf(lz), :) = [];

end

