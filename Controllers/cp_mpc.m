classdef cp_mpc < mpc_hpc & CP
	%CP_MPC Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
        ylmp_constraints;                % save yalmip constraints 
        ylmp_objective;                  % save yalmip objective
        ylmp_opt_variables;              % save the yalmip optimization variables used
        ylmp_opt_output;                 % save the yalmip variable that is extracted after optimization

        ops;                             % options for calling osqp

        mz;                              % QP problem formulation matrices and vectors
        osqps;                           % osqp solver object - used for warm-starting        
        dimx;                            % MPC state dimension
        dimu;                            % MPC output dimension

        cutoffA;                         % sparsify cutoff for discretized state-to-state evolution matrix
        cutoffB;                         % "        "      for discretized input-to-state "         "

        Fdiscretization_step;
	end
	
	methods
		function obj = cp_mpc(hpc)
			%CP_MPC Construct an instance of this class
			%   Detailed explanation goes here

            %%% Object Initialization %%
            % Call superclass constructor before accessing object
            % You cannot conditionalize this statement
            obj = obj@CP(hpc);
            obj = obj@mpc_hpc();
		
            obj.ops = osqp.default_settings();
            obj.ops.warm_start = 0;
            obj.cutoffA = 0;
            obj.cutoffB = 0;
		end
		
	end

	methods
		function obj = setup_mpc(obj)
			
			% It's good practice to start by clearing YALMIPs internal database 
			% Every time you call sdpvar etc, an internal database grows larger
			yalmip('clear')

            obj.dimx = obj.lNs;
            obj.dimu = obj.lNc;
			x = sdpvar(repmat(obj.lNs,1,obj.Nhzn+1),repmat(1,1,obj.Nhzn+1));
            if obj.Nhzn == 1
                u = sdpvar(obj.lNc,1);
            else
			    u = sdpvar(repmat(obj.lNc,1,obj.Nhzn),repmat(1,1,obj.Nhzn));
            end
			
			ot = sdpvar(1,1); %assume ambient temperatre horizon independent

			ly_uref = sdpvar(obj.lNc,1);
            ly_usum = sdpvar(1,1); % total power budget

			Bu = obj.Bd_ctrl(:,1:obj.lNc);
			Bd = obj.Bd_ctrl(:,obj.lNc+1:end);
					
			constraints = [];
			objective = 0;

            if obj.ops.warm_start
                % dummy input constraints, requried for OSQP as inputs must be parametric.
                %% This is required to avoid refactorization of KKT linear system matrix.
                %% set to weird constant to verify index placement in l-vector by Yalmip
                constraints = [constraints, x{1} == 273+80];
                constraints = [constraints, ot == 273+25];
                constraints = [constraints, ly_uref == 18];
                constraints = [constraints, ly_usum == 17*obj.lNc];
            end
            
            if obj.Nhzn == 1
                constraints = [constraints, x{2} == obj.Ad_ctrl*x{1}+Bu*u+Bd*ot];
			    constraints = [constraints, obj.Cc*x{2} <= obj.T_target - obj.Cty(1,:)'];
			    constraints = [constraints,sum(u) <= ly_usum - sum(obj.Ctu(1,:),2)];
				objective = objective + (u-ly_uref)'*obj.Rt*(u-ly_uref) + u'*obj.Rs*(u) + x{2}'*obj.Q*x{2};
            else
			    for k = 1 : obj.Nhzn
			        constraints = [constraints, x{k+1} == obj.Ad_ctrl*x{k}+Bu*u{k}+Bd*ot];
			        constraints = [constraints, obj.Cc*x{k+1} <= obj.T_target - obj.Cty(k,:)'];
			        constraints = [constraints,sum(u{k}) <= ly_usum - sum(obj.Ctu(k,:),2)];

                    %TODO: Ct I did only for y, and only for max
					if (~isinf(obj.umin))
						constraints = [constraints, u{k} >= obj.umin];
					end
					if (~isinf(obj.umax))
						constraints = [constraints, u{k} <= obj.umax - obj.Ctu(k,:)'];
					end

				    objective = objective + (u{k}-ly_uref)'*obj.Rt*(u{k}-ly_uref) + u{k}'*obj.Rs*(u{k}) + ...
                        ((x{k+1}-obj.T_amb)/(obj.T_target(1)-obj.T_amb))'*obj.Q*((x{k+1}-obj.T_amb)/(obj.T_target(1)-obj.T_amb)); %x{k+1}'*obj.Q*x{k+1};
                end
            end
			

			%ops = sdpsettings('verbose',1,'solver','quadprog', 'usex0',1);
			lops = sdpsettings('verbose',1,'solver','osqp', 'usex0',0); %You have specified an initial point, but the selected solver (OSQP) does not support warm-starts through YALMIP
			lops.quadprog.TolCon = 1e-2;
			lops.quadprog.TolX = 1e-5;
			lops.quadprog.TolFun = 1e-3;
			lops.convertconvexquad = 0;
			%lquadprog.MaxPCGIter = max(1, ops.quadprog.MaxPCGIter * 3);
			lops.quadprog.MaxIter = 50;

			%ops.warmstart = 1;

			%ops.usex0 = 1;
			%{
			for i=1:(obj.Nhzn+1)
				assign(x{i}, (273.15+25)*ones(obj.lNs,1));
			end
			for i=1:(obj.Nhzn)
				assign(u{i}, (2)*ones(obj.lNc,1));
			end
			%}

			%just for Andrino:s
			lops.savesolveroutput = 1;
			lops.osqp.rho = 0.1; %0.1 %0.04
			lops.osqp.eps_abs = 0.01; %0.01
			lops.osqp.eps_rel = 0.01; %0.01
			lops.osqp.check_termination = 1;
			%ops.osqp.max_iter = 17;
			lops.osqp.warm_start = 1;

            % save yalmip model for potential extraction to external solver or code generation

			
            % extract model from yalmip to QP formulation
            obj.ylmp_opt_variables = {x{1},ot,ly_uref,ly_usum};
            obj.ylmp_opt_output = {u{1},x{2}};
            obj.ylmp_constraints = constraints;
            obj.ylmp_objective = objective;
            if obj.ops.warm_start
                %lops.osqp = obj.ops;
                ymod = export(constraints,objective,sdpsettings('solver','osqp'));
                % save model .mat file
                lmz = struct(); % save in mazomeros file format
                % cut away the -inf to inf constraints yalmip adds at the end of constraint matrix
                numconstr = length(sdpvar(constraints));
                ymod.A = ymod.A(1:numconstr,:);
                ymod.l = ymod.l(1:numconstr,:);
                ymod.u = ymod.u(1:numconstr,:);



                % setup in maros meszaros problem format
                [lmz.m,lmz.n] = size(ymod.A);
                lmz.P = ymod.P;
                lmz.q = ymod.q;
                lmz.r = 0;
                lmz.l = ymod.l;
                lmz.u = ymod.u;
                lmz.A = ymod.A;
                obj.mz = lmz;

                % call python for optimizing (potentially obmitable)
                %filename = '/tmp/model.mat';
                %fileopt = '/tmp/model_opt.mat';
                %save(filename, '-struct', 'mz');
                %!../optimize_mat2mat.py /tmp/model.mat /tmp/model_opt.mat
                %obj.mz = load(fileopt);

                % setup osqp solver
                obj.osqps = osqp();
                obj.osqps.setup(obj.mz.P, obj.mz.q, obj.mz.A, obj.mz.l, obj.mz.u ,obj.ops);
            else
                % generate optimizer
                % OSQP) does not support warm-starts through YALMIP
				%{
			    yops = sdpsettings('verbose',1,'solver','osqp', 'usex0',0);
			    yops.quadprog.TolCon = 1e-2;
			    yops.quadprog.TolX = 1e-5;
			    yops.quadprog.TolFun = 1e-3;
			    yops.convertconvexquad = 0;
			    yops.quadprog.MaxIter = 50;
                obj.mpc_ctrl = optimizer(constraints,objective,yops,obj.ylmp_opt_variables,obj.ylmp_opt_output);
				%}
				obj.mpc_ctrl = optimizer(constraints,objective,lops,obj.ylmp_opt_variables,obj.ylmp_opt_output);
				%obj.mpc_ctrl = optimizer(constraints,objective,ops,{x{1},ot,ly_uref,ly_usum},{u{1}, x{2}});
			    obj.mpc_ctrl
            end

			
		end %lin_mpc_setup
		function uout = call_mpc(obj, x, ot, uref, usum)

            if ~obj.ops.warm_start
			    [uout, err,~,~,optimizer_object, sol] = obj.mpc_ctrl({x, ot, uref, usum});
            else
                % setup input to osqp by setting equality constraints
                % constraints = [x0,ot,uref,usum, state evolution,....]
                linput = obj.mz.l;
                uinput = obj.mz.u;
                linput(1:obj.dimx) = x;
                uinput(1:obj.dimx) = x;
                linput(obj.dimx+1) = ot;
                uinput(obj.dimx+1) = ot;
                linput(obj.dimx+2:obj.dimx+1+obj.dimu) = uref;
                uinput(obj.dimx+2:obj.dimx+1+obj.dimu) = uref;
                linput(obj.dimx+obj.dimu+2) = usum;
                uinput(obj.dimx+obj.dimu+2) = usum;
                obj.osqps.update('l',linput,'u',uinput);

                % solve with osqp
                solwarm = obj.osqps.solve();

                % extract control input u0 and future temperature state x1 from primal
                %% solwarm.x = xprim = [x0,x1,x2,u0,u1,ot,u_T,P_budget]
                x1 = solwarm.x(obj.dimx+1:2*obj.dimx);
                u0 = solwarm.x((1+obj.Nhzn)*obj.dimx+1:(1+obj.Nhzn)*obj.dimx+obj.dimu);
                uout = cell(2,1);
                uout{1} = u0;
                uout{2} = x1;
            end

			% Save solver stats
			if obj.save_solver_stats
                stats = {};
                stats.x0 = x;
                stats.uref = uref;
                stats.ot = ot;
                stats.usum = usum - sum(obj.Ctu(1,:),2);
                if obj.ops.warm_start
                    stats.lin = linput;
                    stats.uin = uinput;
                    stats.prim = solwarm.x;
                    stats.dua = solwarm.y;
                    stats.info = solwarm.info;
                    stats.x1 = x1;
                    stats.u0 = u0;
                else
                    stats.sol = sol;
                    stats.u0 = uout{1};
                    stats.x1 = uout{2};
                end
				obj.solver_stats = [obj.solver_stats stats];
			end
		end %lin_mpc
	end

	%TODO don't know if these go here!!! or inside MPC!
	methods
		function [obj, comms] = init_fnc(obj, hc, chip, ctrl_id, Nsim)

			obj.initialize(chip, Nsim);

            % sparsify
            if obj.cutoffA ~= 0
                nnzS = nnz(obj.Ad_ctrl);
                obj.Ad_ctrl(obj.Ad_ctrl < obj.cutoffA) = 0;
                nnzE = nnz(obj.Ad_ctrl);
                disp(sprintf('sparsifying A by x%f.  nnz(A) goes from %d to %d.',nnzS/nnzE, nnzS, nnzE));
            end
            if obj.cutoffB ~= 0
                nnzS = nnz(obj.Bd_ctrl);
                obj.Bd_ctrl(obj.Bd_ctrl < obj.cutoffB) = 0;
                nnzE = nnz(obj.Bd_ctrl);
                disp(sprintf('sparsifying B by x%f.  nnz(B) goes from %d to %d.',nnzS/nnzE, nnzS, nnzE));
            end

            obj.umin = obj.lPmin;
            obj.umax = obj.lPmax;

			obj = obj.setup_mpc();

			obj.xlplot = zeros(Nsim+1, obj.lNs);
			obj.xlplot(1,:) = obj.T_amb*ones(obj.lNs,1);

			obj.tmpc = zeros(Nsim+1, obj.lNs);

			obj.failed = 0;
			obj.wl = [ones(obj.lNc,1) zeros(obj.lNc, obj.lipl -1)];

			obj.pw_storage = 0;
			obj.pw_adapt = 0;
			obj.pw_old = [];
			obj.pw_old{1} = 1*obj.lNc;
			obj.pw_old{2} = 1*obj.lNc;

			obj.pbold = 0;
			obj.pbc = 0;

			obj.wl = [ones(obj.lNc,1) zeros(obj.lNc, obj.lipl -1)];

			obj.T_target = ones(obj.lNc, 1)*obj.core_Tcrit;

			obj.f_ma = zeros(obj.lNc,1);

			%TODO
			obj.output_mpc = 1*ones(obj.lNc,1);

			obj.solver_stats = [];

            comms{1} = 0;
            comms{2} = 0;
		end
		function [F,V, comm, obj] = ctrl_fnc(obj, f_ref, pwbdg, pvt, i_pwm, i_wl, ctrl_id, ctrl_comm)

			obj.ex_count = obj.ex_count + 1;

			T = pvt{obj.PVT_T};
			process = pvt{obj.PVT_P};

            chippwbdg = pwbdg(2);
            totpwbdg = pwbdg(1);

			% Process Workload
			obj.wl = obj.wl*(1-obj.alpha_wl) + i_wl*obj.alpha_wl;
			Ceff = obj.wl * (obj.pw_ceff)';
						
			% Process Power Budget
			% (?)
	
			% Adapt Measured&Computed Power
			obj.pw_adapt = obj.cp_pw_adapt(obj.pw_adapt, i_pwm, obj.pw_old{1}, obj.pbc);
			obj.pw_old{1} = obj.pw_old{2};

            if obj.pbc==1
                obj.f_ma = zeros(size(obj.f_ma));
            end
			
			% Choose Voltage
			fctrl = f_ref - obj.f_ma;
            fctrl(fctrl<=0) = obj.lFmin;
			FD = diag(fctrl)*obj.lVDom;	
			V = obj.find_dom_sharedV(obj.lFVT, FD, obj.voltage_rule);
			F = f_ref;

			% Compute Power
			if isempty(obj.pw_stat_exp)
                pex = 1;
            else
				maxT = 125+273.15;
				pex = exp(obj.lVDom*V*obj.pw_stat_exp(1) + (min(obj.Cc*T,ones(obj.lNc,1)*maxT))*obj.pw_stat_exp(2) + ones(obj.lNc,1)*obj.pw_stat_exp(3));		
			end
			pu = Ceff.*F.*(obj.lVDom*(V.*V)) + (obj.pw_stat_lin(1).*(obj.lVDom*V) + process*obj.pw_stat_lin(1)).*pex;

			% TODO observer
			obj.Tobs = T;			

			% MPC
			obj.xlplot(obj.ex_count+1,:) = 0;
			%mpc_pw_target = repmat(chippwbdg, obj.Nhzn,1+obj.lvd);
			%res = obj.call_mpc(obj.Tobs, obj.T_amb*1000, pu, mpc_pw_target);
			res = obj.call_mpc(obj.Tobs, obj.T_amb*1000, pu, chippwbdg);
			tt = isnan(res{1});
			% too complex to make it vectorial
			%dp = res{1}(~tt);
			%obj.output_mpc = repmat(tt, length(tt), length(dp))*dp + ...tt.*obj.output_mpc;
			if ~tt
				obj.output_mpc = res{1};
				%TODO
				%else
			end
			obj.tmpc(obj.ex_count+1, :) = res{2}';

			%TODO outmpc_min
			outmpc_min = 0.1;
			obj.output_mpc(obj.output_mpc<outmpc_min) = outmpc_min;
			obj.failed = obj.failed + ((sum(tt)>0) | (sum(obj.output_mpc<outmpc_min)>0));

			% Compute Freq
			if isempty(obj.pw_stat_exp)
                pex = 1;
            else
				maxT = 125+273.15;
				pex = exp(obj.lVDom*V*obj.pw_stat_exp(1) + (min(obj.Cc*T,ones(obj.lNc,1)*maxT))*obj.pw_stat_exp(2) + ones(obj.lNc,1)*obj.pw_stat_exp(3));		
            end
			F_og = (obj.output_mpc - (obj.pw_stat_lin(1).*(obj.lVDom*V) + process*obj.pw_stat_lin(3)).*pex) ./ (obj.lVDom*(V.*V)) ./ Ceff;

			% Process Freq
			%	Check vs maxF, Temp hysteresis, etc.					
			fmaxi = f_ref+0.001;
			F = F_og + (F_og>fmaxi).*(fmaxi - F_og);
			F = F + (F<obj.lFmin).*(obj.lFmin*ones(obj.lNc,1) - F);
			% discretize
			if obj.Fdiscretization_step > 0
				F = fix(F/obj.Fdiscretization_step) * obj.Fdiscretization_step;
			end
			%	Check vs maxF, Temp hysteresis, etc.
			fmaxi = obj.V2FM(obj.lFVT, obj.lVDom * V);
			F = F + (F>fmaxi).*(fmaxi - F);
			F = F + (F<obj.lFmin).*(obj.lFmin*ones(obj.lNc,1) - F);
	
			%if powerbudget has changed || freq changed
			% interpolate a parabola
			%	parabola depends on the changes (powe >> freq, and Delta of changes)
			% f_MA = f_ref - f_app --> sat 0
			%	Sat 0 is a choice, without we have TurboBoost!
			%TODO Turbo boosting is not working!
			turbo_boost = zeros(obj.lNc,1);
			tt = (f_ref - F_og)>-turbo_boost;
			obj.f_ma = obj.f_ma*(1-obj.alpha_ma) + obj.alpha_ma*(tt.*(f_ref - F_og) + ~tt.*0); %This is 0 and not turbo_boost!
			
			%TEST: TODO REMOVE
			pu = Ceff.*F.*(obj.lVDom*(V.*V)) + (obj.pw_stat_lin(1).*(obj.lVDom*V) + process*obj.pw_stat_lin(3)).*pex;
			obj.pw_old{2} = sum(pu);

            % Distributed Algorithm:
            comm{1} = 0;
            comm{2} = 0;
            aa = 0;
            bb= 0;
            if ctrl_id==1
                adab = 1;
            else
                adab = 2;
            end
            for i=1:length(ctrl_comm)
                if ~isempty(ctrl_comm{i})
                    aa = aa + ctrl_comm{i}{1} + adab;
                    bb = bb + ctrl_comm{i}{2} + 1;
                end
            end
            comm{1} = aa;
            comm{2} = bb;

		end

		function [obj] = cleanup_fnc(obj,hc)
		end
		function [obj] = plot_fnc(obj, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop)
			disp(strcat('[CTRL][MPC] number of times the optimization algorithm failed: ',int2str(obj.failed), '/', int2str(obj.lNsim)));
				
			figure();
			%plot(obj.Ts*sim_mul*[1:Nsim]', pceff(2:end,:)*20+293, 'g'); hold on;
			%plot(t1, cpxplot(2:end,:) - 273, 'b'); hold on; grid on;
			%plot(t2, [NaN*ones(1,size(obj.tmpc,2)); (obj.tmpc(2:end-1,:) - 273)], 'm'); hold on;
			plot(t1, [(cpxplot(2+2:end,:) - 273); NaN*ones(2,size(obj.tmpc,2))], 'b'); hold on;
			plot(t2, obj.tmpc(2:end,:) - 273.15, 'm'); hold on;
			xlabel("Time [s]");
			ylabel("Temperature [T]");
		end
	end

	methods
		function res = intl_bis_fnc(obj, in, args)
			res = obj.pw_function(in, args, 1, 0, 0);
		end
	end
end

