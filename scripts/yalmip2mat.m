function yalmip2mat (filename,constraints,objective,ops,opt_variables,opt_output)
% DUMP MATLAB SPARSE MATRICES TO FILE
    % check file name ending
    assert(all(filename([end-3:end]) == '.mat'));
    %% EXTRACT MODEL
	% ops = sdpsettings('verbose',1,'solver','osqp', 'usex0',0); % warm start unsupported through yalmip interface
    yalmip_model = export(constraints,objective,ops,opt_variables,opt_output);
    ymod = yalmip2osqp(yalmip_model);

    % Save Problem in Maros Meszaros Format
    % min_x  1/2 x^T*P*x + q^T*x + r
    % s.t    l <= A*x <= u
    mazomeros = struct();
    [mazomeros.m,mazomeros.n] = size(ymod.A);
    mazomeros.P = ymod.P;
    mazomeros.q = ymod.q;
    mazomeros.r = 0;
    mazomeros.l = ymod.l;
    mazomeros.u = ymod.u;
    mazomeros.A = ymod.A; %TODO: yalmip inflates inequalities by a factor of 2
    
    %% DUMP MODEL
    % dump header
    save(filename, '-struct', 'mazomeros');
end



%options = ymod.options;
% alternative
%ops = sdpsettings('verbose',0,'solver','osqp','usex0',0, 'savesolverinput',1);
%diagnostics = optimize(constraints,objective,ops);
