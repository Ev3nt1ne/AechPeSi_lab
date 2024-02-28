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
	end
	
	methods
		function obj = cp_mpc()
			%CP_MPC Construct an instance of this class
			%   Detailed explanation goes here
			
            obj.ops = osqp.default_settings();
            obj.ops.warm_start = false;
		end
		
	end

	methods
		function obj = setup_mpc(obj, hc)
			
			% It's good practice to start by clearing YALMIPs internal database 
			% Every time you call sdpvar etc, an internal database grows larger
			yalmip('clear')

            obj.dimx = hc.Ns;
            obj.dimu = hc.Nc;
			x = sdpvar(repmat(hc.Ns,1,obj.Nhzn+1),repmat(1,1,obj.Nhzn+1));
			u = sdpvar(repmat(hc.Nc,1,obj.Nhzn),repmat(1,1,obj.Nhzn));
			
			ot = sdpvar(1,1); %assume ambient temperatre horizon independent

			ly_uref = sdpvar(hc.Nc,1);
            ly_usum = sdpvar(1,1); % total power budget

			Bu = obj.Bd_ctrl(:,1:hc.Nc);
			Bd = obj.Bd_ctrl(:,hc.Nc+1:end);
					
			constraints = [];

            if obj.ops.warm_start
                % dummy input constraints, requried for OSQP as inputs must be parametric.
                %% This is required to avoid refactorization of KKT linear system matrix.
                %% set to weird constant to verify index placement in l-vector by Yalmip
                constraints = [constraints, x{1} == 0];
                constraints = [constraints, ot == 1];
                constraints = [constraints, ly_uref == 2];
                constraints = [constraints, ly_usum == 3];
            end

			for k = 1 : obj.Nhzn
			    constraints = [constraints, x{k+1} == obj.Ad_ctrl*x{k}+Bu*u{k}+Bd*ot];
			    constraints = [constraints, hc.Cc*x{k+1} <= obj.T_target - obj.Cty(k,:)'];
			    constraints = [constraints,sum(u{k}) <= ly_usum - sum(obj.Ctu(k,:),2)];
		    end
			
			objective = 0;
			for k = 1:obj.Nhzn
				objective = objective + (u{k}-ly_uref)'*obj.R*(u{k}-ly_uref) + u{k}'*obj.R2*(u{k}) + x{k+1}'*obj.Q*x{k+1};
			end

			%ops = sdpsettings('verbose',1,'solver','quadprog', 'usex0',1);
			ops = sdpsettings('verbose',1,'solver','osqp', 'usex0',0); %You have specified an initial point, but the selected solver (OSQP) does not support warm-starts through YALMIP
			ops.quadprog.TolCon = 1e-2;
			ops.quadprog.TolX = 1e-5;
			ops.quadprog.TolFun = 1e-3;
			ops.convertconvexquad = 0;
			%ops.quadprog.MaxPCGIter = max(1, ops.quadprog.MaxPCGIter * 3);
			ops.quadprog.MaxIter = 50;
			
            % extract model from yalmip to QP formulation
            obj.ylmp_opt_variables = {x{1},ot,ly_uref,ly_usum};
            obj.ylmp_opt_output = {u{1},x{2}};
            obj.ylmp_constraints = constraints;
            obj.ylmp_objective = objective;
            if obj.ops.warm_start
                ops.osqp = obj.ops;
                ymod = export(constraints,objective,ops);
                % save model .mat file
                mz = struct(); % save in mazomeros file format
                % cut away the -inf to inf constraints yalmip adds at the end of constraint matrix
                numconstr = length(sdpvar(constraints));
                ymod.A = ymod.A(1:numconstr,:);
                ymod.l = ymod.l(1:numconstr,:);
                ymod.u = ymod.u(1:numconstr,:);

                % setup in maros meszaros problem format
                [mz.m,mz.n] = size(ymod.A);
                mz.P = ymod.P;
                mz.q = ymod.q;
                mz.r = 0;
                mz.l = ymod.l;
                mz.u = ymod.u;
                mz.A = ymod.A;
                obj.mz = mz;

                % call python for optimizing (potentially obmitable)
                %filename = '/tmp/model.mat';
                %fileopt = '/tmp/model_opt.mat';
                %save(filename, '-struct', 'mz');
                %!../optimize_mat2mat.py /tmp/model.mat /tmp/model_opt.mat
                %obj.mz = load(fileopt);

                % setup osqp solver
                obj.osqps = osqp();
                obj.osqps.setup(obj.mz.P, obj.mz.q, obj.mz.A, obj.mz.l, obj.mz.u ,ymod.options());
            else
                % generate optimizer
                obj.mpc_ctrl = optimizer(constraints,objective,ops,obj.ylmp_opt_variables,obj.ylmp_opt_output);
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
		function [obj] = init_fnc(obj, hpc_class, Nsim)

			disc = obj.discreate_system(hpc_class, obj.Ts_ctrl);
			obj.Ad_ctrl = disc.A;
			obj.Bd_ctrl = disc.B;

			obj = obj.setup_mpc(hpc_class);

			obj.xlplot = zeros(Nsim+1, hpc_class.Ns);
			obj.xlplot(1,:) = hpc_class.t_init;

			obj.tmpc = zeros(Nsim+1, hpc_class.Ns);

			obj.ex_count = 0;
			obj.failed = 0;
			obj.wl = [ones(hpc_class.Nc,1) zeros(hpc_class.Nc, hpc_class.ipl -1)];
			obj.lNsim = Nsim;	

			obj.pw_storage = 0;
			obj.pw_adapt = 0;
			obj.pw_old = [];
			obj.pw_old{1} = 1*hpc_class.Nc;
			obj.pw_old{2} = 1*hpc_class.Nc;

			obj.pbold = 0;
			obj.pbc = 0;

			obj.wl = [ones(hpc_class.Nc,1) zeros(hpc_class.Nc, hpc_class.ipl -1)];

			obj.T_target = ones(hpc_class.Nc, 1)*hpc_class.core_limit_temp;

			obj.f_ma = zeros(hpc_class.Nc,1);

			%TODO
			obj.output_mpc = 1*ones(hpc_class.Nc,1);

			obj.solver_stats = [];

		end
		function [F,V,obj] = ctrl_fnc(obj, hc, target_index, pvt, i_pwm, i_wl)

			obj.ex_count = obj.ex_count + 1;

			f_ref = hc.frtrc(min(target_index, size(hc.frtrc,1)),:)';
			p_budget = hc.tot_pw_budget(min(target_index, length(hc.tot_pw_budget)));
			T = pvt{hc.PVT_T};
			process = pvt{hc.PVT_P};

			% Process Workload
			obj.wl = obj.wl*(1-obj.alpha_wl) + i_wl*obj.alpha_wl;
			Ceff = obj.wl * (hc.dyn_ceff_k)';
						
			% Process Power Budget
			% (?)
	
			% Adapt Measured&Computed Power
			obj.pw_adapt = obj.cp_pw_adapt(obj.pw_adapt, i_pwm, obj.pw_old{1}, obj.pbc);
			obj.pw_old{1} = obj.pw_old{2};
			
			% Choose Voltage
			FD = diag(f_ref-obj.f_ma)*hc.VDom;
			V = obj.cp_voltage_choice(hc, FD);
			F = f_ref;

			% Compute Power
			if hc.leak_exp
				maxT = 125+273.15;
				pex = exp(hc.VDom*V*hc.leak_exp_vdd_k + (min(hc.Cc*T,ones(hc.Nc,1)*maxT)-273.15)*hc.leak_exp_t_k + ones(hc.Nc,1)*hc.leak_exp_k);
			else
				pex = 1;
			end
			pu = Ceff.*F.*(hc.VDom*(V.*V)) + (hc.leak_vdd_k.*(hc.VDom*V) + process*hc.leak_process_k).*pex;

			% TODO observer
			Tobs = T;			

			% MPC
			obj.xlplot(obj.ex_count+1,:) = 0;
			%mpc_pw_target = repmat(p_budget, obj.Nhzn,1+hc.vd);
			%res = obj.call_mpc(Tobs, hc.temp_amb*1000, pu, mpc_pw_target);
			res = obj.call_mpc(Tobs, hc.temp_amb*1000, pu, p_budget);
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
			if hc.leak_exp
				maxT = 125+273.15;
				pex = exp(hc.VDom*V*hc.leak_exp_vdd_k + (min(hc.Cc*res{2},ones(hc.Nc,1)*maxT)-273.15)*hc.leak_exp_t_k + ones(hc.Nc,1)*hc.leak_exp_k);
			else
				pex = 1;
			end
			F_og = (obj.output_mpc - (hc.leak_vdd_k.*(hc.VDom*V) + process*hc.leak_process_k).*pex) ./ (hc.VDom*(V.*V)) ./ Ceff;

			% Process Freq
			%	Check vs maxF, Temp hysteresis, etc.					
			fmaxi = f_ref+0.001;
			F = F_og + (F_og>fmaxi).*(fmaxi - F_og);
			F = F + (F<hc.F_min).*(hc.F_min*ones(hc.Nc,1) - F);
			% discretize
			if hc.F_discretization_step > 0
				F = fix(F/hc.F_discretization_step) * hc.F_discretization_step;
			end
	
			%if powerbudget has changed || freq changed
			% interpolate a parabola
			%	parabola depends on the changes (powe >> freq, and Delta of changes)
			% f_MA = f_ref - f_app --> sat 0
			%	Sat 0 is a choice, without we have TurboBoost!
			%TODO Turbo boosting is not working!
			turbo_boost = zeros(hc.Nc,1);
			tt = (f_ref - F_og)>-turbo_boost;
			obj.f_ma = obj.f_ma*(1-obj.alpha_ma) + obj.alpha_ma*(tt.*(f_ref - F_og) + ~tt.*0); %This is 0 and not turbo_boost!
			
			%TEST: TODO REMOVE
			pu = Ceff.*F.*(hc.VDom*(V.*V)) + (hc.leak_vdd_k.*(hc.VDom*V) + process*hc.leak_process_k).*pex;
			obj.pw_old{2} = sum(pu);

		end

		function [obj] = cleanup_fnc(obj, hpc_class)
		end
		function [obj] = plot_fnc(obj, hc, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop)
			disp(strcat('[CTRL][MPC] number of times the optimization algorithm failed: ',int2str(obj.failed), '/', int2str(obj.lNsim)));
				
			figure();
			%plot(obj.Ts*sim_mul*[1:Nsim]', pceff(2:end,:)*20+293, 'g'); hold on;
			plot(t1, cpxplot(2:end,:) - 273, 'b'); hold on; grid on;
			plot(t2, [NaN*ones(1,size(obj.tmpc,2)); (obj.tmpc(2:end-1,:) - 273)], 'm'); hold on;
			%plot(t2, [(obj.tmpc(2+2:end,:) - 273); NaN*ones(2,size(obj.tmpc,2))], 'm'); hold on;
			%plot(t2, obj.tmpc(2:end,:) - 273.15, 'm'); hold on;
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

