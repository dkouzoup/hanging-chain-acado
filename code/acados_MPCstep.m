function output = acados_MPCstep( input )
    import acados.*
    
    persistent qp;
    
    solver = input.solver;
%     if strfind(input.solver, 'qpOASES')
%         
%         solver = 'qpoases';
%         
%     elseif strfind(input.solver, 'qpDUNES')
%         
%         solver = 'qpdunes';
%         
%     elseif strfind(input.solver, 'HPMPC')
%         
%         solver = 'hpmpc';
%         
%     end
         
    %% extract data from input
    
    NX           =  size(input.x,2);
    NU           =  size(input.u,2);
    N            =  size(input.u,1);
    M            =  (NX/3-1)/2;
    Q            =  input.W(1:NX,1:NX);
    R            =  input.W(NX+1:end,NX+1:end);
    QN           =  input.WN;
    state_traj   =  input.x';
    control_traj =  input.u';
    X_ref        =  [input.y(:,1:NX);input.yN]';
    U_ref        =  input.y(:,NX+1:end)';
    lbu          =  input.LBU;
    ubu          =  input.UBU;
    WALL         =  input.WALL;
   
    lbx          =  -inf*ones(NX,1);
    ubx          =  inf*ones(NX,1);
    lbx(2:3:3*(M+1)) = WALL;
    
    %% evaluate sensitivities
    
    SensX = cell(N,1);
    SensU = cell(N,1);
    SensC = cell(N,1);
    GradX = cell(N+1,1);
    GradU = cell(N,1);
 
    INT_NAME = sprintf('integrate_chain_M%d',M+2);
    for ii = 1:N
        sim_input.x = state_traj(:,ii);
        sim_input.u = control_traj(:,ii);
        [states,out] = eval([ INT_NAME,'(sim_input)']);
        if out.errorCode
            error('integrator is failed')
        end
        SensX{ii,1} = states.sensX; 
        SensU{ii,1} = states.sensU;
        SensC{ii,1} = states.value - state_traj(:,ii+1);
        GradX{ii,1} = Q*(sim_input.x-X_ref(:,ii));
        GradU{ii,1} = R*(sim_input.u-U_ref(:,ii));
    end
    GradX{N+1,1} = QN*(state_traj(:,N+1)-X_ref(:,N+1));


    %% solve with acados
    
    if isempty(qp)
        warning('first iteration')
        qp = ocp_qp(N, NX, NU);
   
        qp.set('Q', Q);
        qp.set('R', R);
        qp.set('Q', N, QN);
        
        opts = struct('max_iter', 1000, 'warm_start', input.warmstart);
        
        if strcmp(solver, 'qpdunes') || strcmp(solver, 'hpmpc')
            if ~isinf(input.N2)
                opts.N2 = input.N2;
            else
                opts.N2 = N;
                if strcmp(solver, 'qpdunes')
                    opts.clipping = 1;
                end
            end
        end
        
        % TODO: FIX SEGFAULTS IN HPMPC!!
        if strcmp(solver, 'hpmpc')
            opts.mu0 = 1e6;
            % opts.tol = 1e-12; % TODO: WHY SETTING THIS VIA OPTS BREAKS SOLVER?
        end
        
        qp.initialize_solver(solver, opts);
    end
    
    for ii = 0:N-1
        qp.set('A', ii, SensX{ii+1});
        qp.set('B', ii, SensU{ii+1});
        qp.set('b', ii, SensC{ii+1}); 
        qp.set('q', ii, GradX{ii+1});
        qp.set('r', ii, GradU{ii+1});        
    end
    qp.set('q', N, GradX{N+1});
    
    % specify initial condition
    for ii = 0:N-1
        qp.set('lbu', ii, lbu - control_traj(:,ii+1));
        qp.set('ubu', ii, ubu - control_traj(:,ii+1));
        qp.set('lbx', ii+1, lbx - state_traj(:,ii+2));
        qp.set('ubx', ii+1, ubx - state_traj(:,ii+2));
    end
    qp.set('lbx', 0, input.x0' - state_traj(:,1));
    qp.set('ubx', 0, input.x0' - state_traj(:,1));
    
    % solve QP
    acados_output = qp.solve();
    
    Delta_xopt = acados_output.states();
    Delta_uopt = acados_output.controls();

    %% write output
    
    output.x   = (state_traj   + [Delta_xopt{:}])';
    output.u   = (control_traj + [Delta_uopt{:}])';
    output.info.QP_time = acados_output.info.total_time;
    output.info.nIterations = acados_output.info.num_iter;

end

