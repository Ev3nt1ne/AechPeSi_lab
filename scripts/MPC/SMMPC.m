classdef SMMPC < handle
    properties
        sys         % LinearSystem instance defining dynamics & constraints
        optimizer   % Yalmip optimizer
        N           % Prediction horizon
        Q           % state stage cost weight matrix x_i^T Q x
        R           % input stage cost weight matrix u_i^T R u
        P           % terminal cost weight matrix u_i^T P u
        alpha       % terminal set level for level set of terminal cost: 
                    % Xf = { x | x^T P x <= alpha} 
        K           % Tube control gain
        sm          % Set membership estimtor to update performance model
        x_prev      % Previous state measurement
        u_prev      % Previous input measurement
    end
    
    methods
        function obj = SMMPC(sys,N,Q,R,P,alpha,K,sm)
            obj.sys = sys;
            obj.N = N;
            obj.Q = Q;
            obj.R = R;
            obj.P = P;
            obj.alpha = alpha;
            obj.K = K;
            obj.sm = sm;
            obj.x_prev = [];
            obj.u_prev = [];

            
            x_i=sdpvar(sys.n, N+1); %performance state
            z_i=sdpvar(sys.n, N+1); %nominal state
            u_i=sdpvar(sys.m, N); %input
            x_0=sdpvar(sys.n,1); %initial state
            PxA=sdpvar(sys.nx, sys.n, N+1,'full'); %A matrices of state constraints
            PuA=sdpvar(sys.nu,sys.m, N,'full'); %A Matrices of input constraints
            Pxb=sdpvar(sys.nx,N+1,'full');  %b vectors of state constraints
            Pub=sdpvar(sys.nu,N,'full');    %b vectors of input constraints
            Ap=sdpvar(sys.n,sys.n,'full'); %A Matrix
            Bp=sdpvar(sys.n,sys.m,'full'); %b Matrix

            % Define Cost Function & Constraints
            objective_MPC=0;
            
            % Initial State Constraints
            constraints_MPC=[x_i(:,1)==x_0];
            
            for i=1:N
                % Stage cost
                objective_MPC=objective_MPC+x_i(:,i)'*Q*x_i(:,i)+ u_i(:,i)'*R*u_i(:,i);
                
                % State Propagation Constraints
                constraints_MPC=[constraints_MPC, x_i(:,i+1) == Ap*x_i(:,i)+Bp*u_i(:,i)];

                % State Constraints
                constraints_MPC=[constraints_MPC, PxA(:,:,i+1)*x_i(:,i+1)<=Pxb(:,i+1)];
                constraints_MPC=[constraints_MPC, PuA(:,:,i)*u_i(:,i)<=Pub(:,i)];
            end

            % Terminal Cost
            objective_MPC=objective_MPC+x_i(:,N+1)'*P*x_i(:,N+1);

            % Terminal constraint
            constraints_MPC=[constraints_MPC, x_i(:,N+1)'*P*x_i(:,N+1)<=alpha];

            % Optimizer allows to solve optimisation problem repeatedly in a fast
            % manner. Inputs are x_0, Ap,Bp, PxA, Pxb, PuA, Pub and outputs are u_i, x_i
            % To speed up the solve time, MOSEK can be used with an academic license.
            % Replace [] with sdpsettings('solver','mosek') if installed. Other solvers
            % can be used as well (OSQP,...)
            ops = sdpsettings('verbose',1,'solver','sedumi');
            %ops = sdpsettings('verbose',1,'solver','mosek');
            obj.optimizer=optimizer(constraints_MPC, objective_MPC, ops, {x_0, Ap,Bp, PxA, Pxb, PuA, Pub}, {u_i,x_i});
        end
        
        function [u, U, X] = solve(obj, x)
            if ~isempty(obj.x_prev) %update Estimator
                obj.sm.update(x, obj.x_prev, obj.u_prev);
            end
                     
            % Define current disturbance set
            W = obj.sm.W + obj.sm.estimateW_theta(obj.sys.Px,obj.sys.Pu);
            
            % Compute forward reachable sets
            F{1} = Polyhedron(zeros(0,2),[]);
            for i=1:obj.N
                F{i+1} = (obj.sys.A+obj.sys.B*obj.K)*F{i} + W;
                aux = obj.sys.Px - F{i+1};
                aux.minHRep;
                PxA(:,:,i+1) = aux.A;
                Pxb(:,i+1) = aux.b;
                aux = obj.sys.Pu - obj.K*F{i};
                aux.minHRep;
                PuA(:,:,i) = aux.A;
                Pub(:,i) = aux.b;
			end
            
            ABp = obj.sm.get_AB;
            Ap = ABp(:,1:obj.sys.n);
            Bp = ABp(:,obj.sys.n+1:end);
            
            % Call optimizer and check solve status
            [sol, flag] = obj.optimizer(x, Ap, Bp, PxA, Pxb, PuA, Pub);
            
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

