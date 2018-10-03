function code_generate_fiodos(N, NX, NU)

%CODE_GENERATE_FIORDOS Code generate fiordos solver with time-varying
% dynamics and bounds (due to relative QP)

% remove any previously generated controller
if isdir('fiordos_controller')
    rmdir('fiordos_controller', 's')
end

mkdir('fiordos_controller')

M = (NX/3 - 1)/2;

% weights
R  = 0.01*eye(3);
Q  = blkdiag(25*eye(3), 25*eye(3*M), 1*eye(3*M));
QN = blkdiag(25*eye(3), 25*eye(3*M), 1*eye(3*M));

%% SET UP PROBLEM

% TODO: CAN WE DO IT WITHOUT INFS?
ZZ = SimpleSet(2*N);

% input bounds
ZZ.addSet(1:N, EssBox(NU, 'l', 'param', 'u', 'param'));

% state bounds (x0 in equality constraints)
ZZ.addSet(N+(1:N),EssBox(NX, 'l', 'param','u', 'param'));

H = blkdiag(kron(eye(N), R),  kron(eye(N-1), Q), QN);

op = OptProb('H', H, 'g', 'param', 'X', ZZ, 'Ae', 'param', 'be', 'param', 'me', N*NX);


%% CODE GENERATE

tol = 1e-6;

cd('fiordos_controller')

s = Solver(op);
s.setSettings('approach', 'stopg', true, 'stopgEpsPrimal', tol, 'stopgEpsDual', tol);
s.generateCode('prefix','fiordos_mpc_', 'forceOverwrite', true);
rehash;
fiordos_mpc_mex_make;

cd('..');
addpath('fiordos_controller')

end

