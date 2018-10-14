function [Ae, be] = build_equality_constraints(A, B, c, x0)

% BUILD_EQUALITY_CONSTRAINTS transform dynamics to form Ae*Z = be with  
%                            Z = [u_0, ..., u_{N-1}, x_1, ..., x_N] 

N    = length(A);
NX   = size(A{1},1);
Aeu  = [];
Aex1 = [];
Aex2 = [];

for ii = 1:N
    Aeu  = blkdiag(Aeu, -B{ii});
    Aex1 = blkdiag(Aex1, eye(NX));
    if ii < N
        Aex2 = blkdiag(Aex2, -A{ii+1});
    end
end

Aex = Aex1 + [zeros(NX, size(Aex1,2)); Aex2 zeros(size(Aex2,1), NX)];
Ae  = [Aeu Aex];
be  = vertcat(c{:});

% x0 embedding
be(1:NX) = be(1:NX) + A{1}*x0;

end

