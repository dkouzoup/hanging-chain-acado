function nlp = acados_setup(NMASS, ode_ca_fun, WALL, N, Ts, W, WN, Xref, SOLVER, WARMSTART)

% ACADOS_SETUP Create acados RTI solver
% TODO: write acados commit or add submodule

rng(0); % what for?

import casadi.*
import acados.*

M  = NMASS - 2;
NX = (2*M + 1)*3;
NU = 3;
BS = extract_block_size(SOLVER);

if isempty(BS)
    N2 = N;
else
    N2 = round(N/BS);
end

try
    nlp = ocp_nlp(N, NX, NU);
catch
    error('acados-matlab is probably not in your path');
end

% TODO: IRK_GL2
nlp.set_dynamics(ode_ca_fun, struct('integrator', 'rk4', 'step', Ts));

% TODO: do not hardcode those bounds
nlp.set_field('lbu', -1*ones(NU,1));
nlp.set_field('ubu', +1*ones(NU,1));

% TODO: why not doing this before init does not update constr. online?
nlp.set_field('lbx', 0, Xref);
nlp.set_field('ubx', 0, Xref);

lbx = -inf*ones(NX,1);
lbx(2:3:3*(M+1)) = WALL;
ubx = inf*ones(NX, 1);

ubx(isinf(ubx)) = 1e6;
lbx(isinf(lbx)) = -1e6;

for ii = 1:N
    nlp.set_field('lbx', ii, lbx);
    nlp.set_field('ubx', ii, ubx);
end

nlp.set_stage_cost(eye(NX + NU), [Xref; zeros(NU, 1)], W);

nlp.set_terminal_cost(eye(NX), Xref, WN);


if contains(SOLVER,'qpOASES_N3')

    error('acados does not support cubic condensing nor C++ qpOASES')

elseif contains(SOLVER,'qpOASES_N2')

    error('acados does not support C++ qpOASES')

elseif contains(SOLVER,'qpOASES_e_N3')

    error('acados does not support cubic condensing')

elseif contains(SOLVER,'qpOASES_e_N2')

    nlp.initialize_solver('rti', struct('qp_solver', 'qpoases'));

    if WARMSTART
        % TODO: warmstart
    end
    
elseif contains(SOLVER,'qpDUNES') && N2 <= 1

    nlp.initialize_solver('rti', struct('qp_solver', 'qpdunes', 'clipping', 1));

    if WARMSTART
        % TODO
    end

elseif contains(SOLVER,'qpDUNES') && N2 > 1

    nlp.initialize_solver('rti', struct('qp_solver', 'qpdunes', 'clipping', 0, 'N2', N2));

    if WARMSTART
        % TODO
    end
    
elseif contains(SOLVER,'HPIPM')

    nlp.initialize_solver('sqp', struct('qp_solver', 'hpipm', 'max_iter', 1, 'N2', N2));

else
    error('SPECIFIED SOLVER DOES NOT EXIST')
end

end

