classdef NominalMPC
    properties
        sys         % LinearSystem instance defining dynamics & constraints
        optimizer   % Yalmip optimizer
        N           % Prediction horizon
        Q           % state stage cost weight matrix x_i^T Q x
        R           % input stage cost weight matrix u_i^T R u
        P           % terminal cost weight matrix u_i^T P u
        alpha       % terminal set level for level set of terminal cost: 
                    % Xf = { x | x^T P x <= alpha} 
    end
    
    methods
        function obj = NominalMPC(sys,N,Q,R,P,alpha)
            % class constructor setting up optimization problem of MPC
            
            obj.sys = sys;
            obj.N = N;
            obj.Q = Q;
            obj.R = R;
            obj.P = P;
            obj.alpha = alpha;
            
            x_i=sdpvar(sys.n, N+1); %state
            u_i=sdpvar(sys.m, N); %input
            x_0=sdpvar(sys.n,1); %initial state
            
            % Define Cost Function & Constraints
            objective_MPC=0;
            
            % Initial State Constraint
            constraints_MPC=[x_i(:,1)==x_0];
                        
            for i=1:N
                % Stage cost
                objective_MPC=objective_MPC+x_i(:,i)'*Q*x_i(:,i)+ u_i(:,i)'*R*u_i(:,i);
                
                % State Propagation Constraints
                constraints_MPC=[constraints_MPC, x_i(:,i+1) == sys.A*x_i(:,i)+sys.B*u_i(:,i)];

                % State & Input Constraints
                constraints_MPC=[constraints_MPC, sys.Px.A*x_i(:,i+1)<=sys.Px.b];
                constraints_MPC=[constraints_MPC, sys.Pu.A*u_i(:,i)<=sys.Pu.b];
            end

            % Terminal Cost
            objective_MPC=objective_MPC+x_i(:,N+1)'*P*x_i(:,N+1);


            % Terminal constraint
            constraints_MPC=[constraints_MPC, x_i(:,N+1)'*P*x_i(:,N+1)<=alpha];         

            % Optimizer allows to solve optimisation problem repeatedly in a fast
            % manner. Inputs are x_0 and outputs are predicted sequences u_i, x_i
            % To speed up the solve time, MOSEK can be used with an academic license.
            % Replace [] with sdpsettings('solver','mosek') if installed. Other solvers
            % can be used as well (OSQP,...)
            ops = sdpsettings('verbose',1,'solver','sedumi');
            %ops = sdpsettings('verbose',1,'solver','mosek');
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

