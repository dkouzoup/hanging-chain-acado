function [ Ae ] = build_equality_constraints_sparsity(N, NX, NU)

% BUILD_EQUALITY_CONSTRAINTS_SPARSITY construct "wost-case" sparsity of 
%                                     equality constraints for dynamics

A = cell(N, 1);
B = cell(N, 1);
c = cell(N,1);

for ii = 1:N
    A{ii} = ones(NX, NX);
    B{ii} = ones(NX, NU);
    c{ii} = ones(NX,1);
end

x0 = ones(NX,1);

[Ae, ~] = build_equality_constraints(A, B, c, x0);

Ae(Ae~=0) = 1;

end

