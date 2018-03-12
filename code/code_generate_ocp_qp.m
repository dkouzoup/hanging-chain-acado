function code_generate_ocp_qp(qpdata, indx, fpath)

if fpath(end) ~= filesep
    fpath(end+1) = filesep;
end

ELIMINATE_X0 = 1;

%% pre-process data

N = qpdata.N;
 
idxb = cell(N+1, 1);
nbx  = zeros(1, N+1);
nbu  = zeros(1, N+1);

% concatenate dynamics

Av = [];
Bv = [];
for ii = 1:N
    Atmp = qpdata.A{ii};
    Btmp = qpdata.B{ii};
    Av   = [Av; Atmp(:)];
    Bv   = [Bv; Btmp(:)];
end
b = vertcat(qpdata.b{:});

% build idxb and nbx/nbu
for ii = 1:N+1
   
    not_inf_bounds_u = ~isinf(qpdata.lbu{ii})|~isinf(qpdata.ubu{ii});
    not_inf_bounds_x = ~isinf(qpdata.lbx{ii})|~isinf(qpdata.ubx{ii});
    
    nbu(ii)  = sum(not_inf_bounds_u);
    nbx(ii)  = sum(not_inf_bounds_x);
    idxb{ii} = find([not_inf_bounds_u; not_inf_bounds_x])-1;
end

nb = nbx+nbu;

% build lb/ub

lbvec = [];
ubvec = [];

for ii = 1:N+1
    lbtmp = [qpdata.lbu{ii}; qpdata.lbx{ii}];
    ubtmp = [qpdata.ubu{ii}; qpdata.ubx{ii}];
    infs  = isinf(lbtmp)&isinf(ubtmp);
    % remove double sided infs
    lbtmp(infs) = [];
    ubtmp(infs) = [];
    % make one sided infs a reasonably large number
    lbtmp(lbtmp < -1e9) = -1e9;
    ubtmp(ubtmp >  1e9) =  1e9;
    
    lbvec = [lbvec; lbtmp];
    ubvec = [ubvec; ubtmp];
end

% check that there are no bounds close to each other (for HPIPM etc)
diff = abs(ubvec-lbvec);
diff(qpdata.nu(1)+1:qpdata.nu(1)+qpdata.nx(1)) = [];
if any(diff(qpdata.nu(1)+1:end) < 1e-3)
    disp('detected upper/lower bounds almost equal to each other');
    keyboard
end

if ELIMINATE_X0    
    % TODO: change d[0] if ng[0] != 0
    
    nx_new = qpdata.nx;
    nx_new(1) = 0;
    
    x0 = qpdata.lbx{1};
    if abs(x0 - qpdata.ubx{1}) > 1e-10
        error('upper and lower bounds are not equal')
    end
    
    nbx_new = nbx;
    nbx_new(1) = nbx_new(1) - qpdata.nx(1);
    nb_new = nb;
    nb_new(1) = nbx_new(1) + nbu(1);
    
    b0 = qpdata.b{1}; 
    b0 = b0 + qpdata.A{1}*x0;
    
    b_new = b;
    b_new(1:length(b0)) = b0;
    
    % TODO make objective cell array
    % TODO test for non-zero S!
    r0 = qpdata.r(1:qpdata.nu(1));
    r0 = r0 + reshape(qpdata.Sv(1:qpdata.nx(1)*qpdata.nu(1)), qpdata.nx(1), qpdata.nu(1))'*x0;
    
    r_new = qpdata.r;
    r_new(1:length(r0)) = r0;
    
    idxb_new = idxb;
    
    idxb_new{1} = idxb_new{1}(1:nbu(1));
   
end

%% write data in .c file

filename = sprintf('ocp_qp_data_nmasses_%d_nsteps_%d_solver_%s_warmstart_%d.c', qpdata.nmasses, qpdata.N, qpdata.solver, qpdata.warmstart);

if indx == 0

    newdatafile = fopen([fpath filename], 'w');

    fprintf(newdatafile, '\n/* Dimensions */\n\n');

    fprintf(newdatafile,'int N = %d;\n\n', qpdata.N);
    code_generate_vec(newdatafile, 'int', 'nx', qpdata.nx);
    code_generate_vec(newdatafile, 'int', 'nu', qpdata.nu);
    code_generate_vec(newdatafile, 'int', 'nbx', nbx);
    code_generate_vec(newdatafile, 'int', 'nbu', nbu);
    code_generate_vec(newdatafile, 'int', 'nb', nb);
    code_generate_vec(newdatafile, 'int', 'ng', zeros(1, N+1));
    code_generate_vec(newdatafile, 'int', 'ns', zeros(1, N+1));

    % CANNOT CHANGE!
    code_generate_vec(newdatafile, 'int', 'idxb', vertcat(idxb{:})); 

    if ELIMINATE_X0
        % updated data
        code_generate_vec(newdatafile, 'int', 'nx_new', nx_new);
        code_generate_vec(newdatafile, 'int', 'nbx_new', nbx_new);
        code_generate_vec(newdatafile, 'int', 'nb_new', nb_new);
        code_generate_vec(newdatafile, 'int', 'idxb_new', vertcat(idxb_new{:})); 
    end
    
else
    newdatafile = fopen([fpath filename], 'a');
end

fprintf(newdatafile, '\n/* Objective */\n\n');

code_generate_vec(newdatafile, 'double', sprintf('Qv_%d',indx), qpdata.Qv); 
code_generate_vec(newdatafile, 'double', sprintf('Rv_%d',indx), qpdata.Rv); 
code_generate_vec(newdatafile, 'double', sprintf('Sv_%d',indx), qpdata.Sv); 
code_generate_vec(newdatafile, 'double', sprintf('q_%d',indx), qpdata.q); 
code_generate_vec(newdatafile, 'double', sprintf('r_%d',indx), qpdata.r);
if ELIMINATE_X0
    % truncated data
    code_generate_vec(newdatafile, 'double', sprintf('Qv_new_%d',indx), qpdata.Qv(qpdata.nx(1)*qpdata.nx(1)+1:end)); 
    code_generate_vec(newdatafile, 'double', sprintf('Sv_new_%d',indx), qpdata.Sv(qpdata.nx(1)*qpdata.nu(1)+1:end)); 
    code_generate_vec(newdatafile, 'double', sprintf('q_new_%d',indx), qpdata.q(qpdata.nx(1)+1:end)); 
    % updated data
    code_generate_vec(newdatafile, 'double', sprintf('r_new_%d',indx), qpdata.r);
end

fprintf(newdatafile, '\n/* Dynamics */\n\n');

code_generate_vec(newdatafile, 'double', sprintf('Av_%d',indx), Av); 
code_generate_vec(newdatafile, 'double', sprintf('Bv_%d',indx), Bv); 
code_generate_vec(newdatafile, 'double', sprintf('b_%d',indx), b); 
if ELIMINATE_X0
    % truncated data
    code_generate_vec(newdatafile, 'double', sprintf('Av_new_%d',indx), Av(qpdata.nx(2)*qpdata.nx(1)+1:end)); 
    % updated data
    code_generate_vec(newdatafile, 'double', sprintf('b_new_%d',indx), b_new);
end

fprintf(newdatafile, '\n/* Constraints */\n\n');

code_generate_vec(newdatafile, 'double', sprintf('lb_%d',indx), lbvec); 
code_generate_vec(newdatafile, 'double', sprintf('ub_%d',indx), ubvec);
if ELIMINATE_X0
    % truncated data
    lbvec_new = lbvec;
    ubvec_new = ubvec;
    lbvec_new(nbu(1)+1:nbu(1)+nbx(1)) = [];
    ubvec_new(nbu(1)+1:nbu(1)+nbx(1)) = [];
    code_generate_vec(newdatafile, 'double', sprintf('lb_new_%d',indx), lbvec_new);
    code_generate_vec(newdatafile, 'double', sprintf('ub_new_%d',indx), ubvec_new);
end

fprintf(newdatafile, '\n/* Solution */\n\n');

code_generate_vec(newdatafile, 'int', sprintf('acado_iter_%d',indx), qpdata.acado_iter); 
code_generate_vec(newdatafile, 'double', sprintf('acado_sol_%d',indx), qpdata.acado_sol);
if ELIMINATE_X0
    % truncated data
    acado_sol_new = qpdata.acado_sol;
    acado_sol_new(qpdata.nu(1)+1:qpdata.nu(1)+qpdata.nx(1)) = [];
    code_generate_vec(newdatafile, 'double', sprintf('acado_sol_new_%d',indx), acado_sol_new);
end

fclose(newdatafile);

end


function code_generate_vec(file_name, field_type, field_name, data)

len = length(data);

fprintf(file_name, '%s %s[%d] = { ', field_type, field_name, len);
for ii = 1:len
    if strcmp(field_type, 'int')
        fprintf(file_name,'%d, ', data(ii));
    elseif strcmp(field_type, 'double')
        fprintf(file_name,'%1.15e, ', data(ii));
    else 
        error('unknown field type, implement if statement.');
    end
end
fprintf(file_name,'};\n\n');

end
