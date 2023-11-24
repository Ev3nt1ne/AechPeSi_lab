classdef PLMPC < handle
    properties
        sys         % LinearSystem instance defining dynamics & constraints
        optimizer   % Yalmip optimizer
        N           % Prediction horizon
        Q           % state stage cost weight matrix x_i^T Q x
        R           % input stage cost weight matrix u_i^T R u
        P           % terminal cost weight matrix u_i^T P u
        alpha       % terminal set level for level set of terminal cost: 
                    % Xf = { x | x^T P x <= alpha} 
        F           % struct with polytopic reachable sets (MPT Polyhedra)
        K           % Tube control gain
        sm          % Set membership estimtor to update performance model
        x_prev      % Previous state measurement
        u_prev      % Previous input measurement
    end
    
    methods
        function obj = PLMPC(sys,N,Q,R,P,alpha,F,K,sm)
            obj.sys = sys;
            obj.N = N;
            obj.Q = Q;
            obj.R = R;
            obj.P = P;
            obj.alpha = alpha;
            obj.F = F;
            obj.K = K;
            obj.sm = sm;
            obj.x_prev = [];
            obj.u_prev = [];

            
            x_i=sdpvar(sys.n, N+1); %performance state
            z_i=sdpvar(sys.n, N+1); %nominal state
            u_i=sdpvar(sys.m, N); %input
            x_0=sdpvar(sys.n,1); %initial state
            Ap=sdpvar(sys.n,sys.n,'full'); %performance A Matrix
            Bp=sdpvar(sys.n,sys.m,'full'); %performance b Matrix

            % Define Cost Function & Constraints
            objective_MPC=0;
            
            % Initial State Constraints
            constraints_MPC=[x_i(:,1)==x_0];
            constraints_MPC=[constraints_MPC, z_i(:,1)==x_0];
            
            for i=1:N
                % Stage cost
                objective_MPC=objective_MPC+x_i(:,i)'*Q*x_i(:,i)+ u_i(:,i)'*R*u_i(:,i);
                %objective_MPC=objective_MPC+z_i(:,i)'*Q*z_i(:,i)+ u_i(:,i)'*R*u_i(:,i);
                
                % State Propagation Constraints
                constraints_MPC=[constraints_MPC, x_i(:,i+1) == Ap*x_i(:,i)+Bp*u_i(:,i)];
                constraints_MPC=[constraints_MPC, z_i(:,i+1) == sys.A*z_i(:,i)+sys.B*u_i(:,i)];


                % State & Input Constraints
                % Use MPT Toolbox to to compute tightened constraints
                X_t{i+1} = sys.Px - F{i+1};
                U_t{i} = sys.Pu - K*F{i};
                constraints_MPC=[constraints_MPC, X_t{i+1}.A*z_i(:,i+1)<=X_t{i+1}.b];
                constraints_MPC=[constraints_MPC, U_t{i}.A*u_i(:,i)<=U_t{i}.b];
            end

            % Terminal Cost
            objective_MPC=objective_MPC+x_i(:,N+1)'*P*x_i(:,N+1);
            %objective_MPC=objective_MPC+z_i(:,N+1)'*P*z_i(:,N+1);

            % Terminal constraint
            constraints_MPC=[constraints_MPC, z_i(:,N+1)'*P*z_i(:,N+1)<=alpha];

            % Optimizer allows to solve optimisation problem repeatedly in a fast
            % manner. Inputs are x_0 and outputs are x_0, Ap,Bp, x_i
            % To speed up the solve time, MOSEK can be used with an academic license.
            % Replace [] with sdpsettings('solver','mosek') if installed. Other solvers
            % can be used as well (OSQP,...)
            %ops = sdpsettings('verbose',1,'solver','mosek');
            ops = sdpsettings('verbose',1,'solver','sedumi');
            obj.optimizer=optimizer(constraints_MPC, objective_MPC, ops, {x_0, Ap, Bp}, {u_i,x_i,z_i});
        end
        
        function [u, U, X] = solve(obj, x)
            
            if ~isempty(obj.x_prev) %update set membership estimator
                obj.sm.update(x, obj.x_prev, obj.u_prev);
            end
          
            ABp = obj.sm.get_AB;
            Ap = ABp(:,1:obj.sys.n);
            Bp = ABp(:,obj.sys.n+1:end);

            % Call optimizer and check solve status
            [sol, flag] = obj.optimizer(x, Ap, Bp);
            
            U = sol{1};
            X = sol{2};
            u = U(:,1);
            if flag~=0
                warning(yalmiperror(flag))
                if any(isnan(u))
                    u = zeros(obj.sys.m,1);
                    warning('MPC returned NaN, overwriting with 0')
                end
            end
            obj.x_prev = x;
            obj.u_prev = u; 
        end
        
        function reset(obj, Omega)
            obj.sm.reset(Omega)
            obj.x_prev = [];
            obj.u_prev = [];
        end
    end
end

