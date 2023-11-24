classdef RobustMPC
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
    end
    
    methods
        function obj = RobustMPC(sys,N,Q,R,P,alpha,F,K)
            obj.sys = sys;
            obj.N = N;
            obj.Q = Q;
            obj.R = R;
            obj.P = P;
            obj.alpha = alpha;
            obj.F = F;
            obj.K = K;
            
            x_i=sdpvar(sys.n, N+1); %state
            u_i=sdpvar(sys.m, N); %input
            x_0=sdpvar(sys.n,1); %initial state


            % Define Cost Function & Constraints
            objective_MPC=0;
            
            % Initial State Constraints
            constraints_MPC=[x_i(:,1)==x_0];
            
            for i=1:N
                % Stage cost
                objective_MPC=objective_MPC+x_i(:,i)'*Q*x_i(:,i)+ u_i(:,i)'*R*u_i(:,i);
                
                % State Propagation Constraints
                constraints_MPC=[constraints_MPC, x_i(:,i+1) == sys.A*x_i(:,i)+sys.B*u_i(:,i)];

                % State Constraints
                % Use MPT Toolbox to to compute tightened constraints
                X_t{i+1} = sys.Px - F{i+1};
                U_t{i} = sys.Pu - K*F{i};
                constraints_MPC=[constraints_MPC, X_t{i+1}.A*x_i(:,i+1)<=X_t{i+1}.b];
                constraints_MPC=[constraints_MPC, U_t{i}.A*u_i(:,i)<=U_t{i}.b];
            end

            % Terminal Cost
            objective_MPC=objective_MPC+x_i(:,N+1)'*P*x_i(:,N+1);

            % Terminal constraint
            constraints_MPC=[constraints_MPC, x_i(:,N+1)'*P*x_i(:,N+1)<=alpha];

            % Optimizer allows to solve optimisation problem repeatedly in a fast
            % manner. Inputs are x_0 and outputs are u_i, x_i
            % To speed up the solve time, MOSEK can be used with an academic license.
            % Replace [] with sdpsettings('solver','mosek') if installed. Other solvers
            % can be used as well (OSQP,...)
            ops = sdpsettings('verbose',1,'solver','sedumi');
            obj.optimizer=optimizer(constraints_MPC, objective_MPC, ops, {x_0}, {u_i,x_i});
        end
        
        function [u, U, X] = solve(obj, x)
            % Call optimizer and check solve status
            [sol, flag] = obj.optimizer(x);
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
        end
    end
end

