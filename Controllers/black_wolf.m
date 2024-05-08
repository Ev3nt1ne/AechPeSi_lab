classdef black_wolf < mpc_hpc & CP
	%UNTITLED Summary of this class goes here
	%   Detailed explanation goes here
	properties
		%alpha_wl = 0.4;				% Moving Average filter parameter for Workload
	end
	properties(SetAccess=protected, GetAccess=public)	
		psac;
		psoff_lut;
		F0v;
		T0v;

		prevF;
		prevV;
		C2;

	end
	
	methods
		function obj = black_wolf()
			%UNTITLED Construct an instance of this class
			%   Detailed explanation goes here
		end		
	end

	methods
		function obj = setup_mpc(obj, hc)
			
			% It's good practice to start by clearing YALMIPs internal database 
			% Every time you call sdpvar etc, an internal database grows larger
			yalmip('clear')

			x = sdpvar(repmat(hc.Ns,1,obj.Nhzn+1),repmat(1,1,obj.Nhzn+1));
			u = sdpvar(repmat(hc.Nc,1,obj.Nhzn),repmat(1,1,obj.Nhzn));
			
			ot = sdpvar(1,1); %TODO: decide if single or with horizon
			w = sdpvar(hc.Nc,1); %TODO: decide if single or with horizon
			h0v = sdpvar(hc.Nc,1);

			%sdpvar ly_xref;
			ly_uref = sdpvar(hc.Nc,1);
			%sdpvar ly_yref;
			ly_usum = sdpvar(obj.Nhzn,hc.vd+1);
			%ly_Cty = dpvar(repmat(obj.Nout,1,obj.Nhzn),repmat(1,1,obj.Nhzn));
			
			h0 = obj.psac(1);
			h1 = obj.psac(2);
			h2 = obj.psac(3);

			if length(h0) < hc.Nc
				h0 = ones(hc.Nc,1)*h0;
			end

			Bu = obj.Bd_ctrl(:,1:hc.Nc);
			Bd = obj.Bd_ctrl(:,hc.Nc+1:end);
					
			constraints = [];
			%constraints = [x{2} == (Adl_mpc + C2*h1)*x{1}-(C2*h1*273.15*ones(hc.Ns,1))+(Bu+(hc.Cc'*h2))*(u0*w)+Bd*ot + hc.Cc*(h0 + h0v)];

			for k = 1 : obj.Nhzn
					%constraints = [constraints, x{k+1} == (obj.Ad_ctrl + obj.C2*h1)*x{k}-(obj.C2*h1*273.15*ones(hc.Ns,1))+(Bu+(hc.Cc'*h2))*(u{k}.*w)+Bd*ot + hc.Cc'*(h0 + h0v)];
					constraints = [constraints, x{k+1} == (obj.Ad_ctrl + Bu*diag(h1)*hc.Cc)*x{k}-(Bu*diag(h1*273.15)*hc.Cc*ones(hc.Ns,1)) ...
									+ Bu*(u{k}.*(w+h2)) + Bd*ot + Bu*h0 + hc.Cc'*h0v];

					%TODO: Ct I did only for y, and only for max
					if (~isinf(obj.umin))
						%constraints = [constraints, u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k}-273.15)) + h0 >= obj.umin];
					end
					if (~isinf(obj.uMax))
						%constraints = [constraints, u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k}-273.15)) + h0 <= obj.uMax - obj.Ctu(k,:)'];
					end
					if (~isinf(obj.xmin))
						%constraints = [constraints, x{k+1} >= obj.xmin];
					end
					if (~isinf(obj.xMax))
						%constraints = [constraints, x{k+1} <= obj.xMax];
					end
					if (~isinf(obj.ymin))
						%constraints = [constraints, obj.C)*x{k+} >= obj.ymin];
					end
					%if (~isinf(obj.yMax))
						constraints = [constraints, hc.Cc*x{k+1} <= obj.T_target - obj.Cty(k,:)'];
					%end
					%if (~isinf(obj.usum(1)))
						constraints = [constraints,sum(u{k}.*w*(h2+1) + h1*(hc.Cc*(x{k}-273.15)) + h0v+h0) <= ly_usum(k,1) - sum(obj.Ctu(k,:),2)];
						%TODO should I also add the "resulting" thing. i.e.
						%	the above formula with x{k+1}?
						%	maybe no, because I use horizon for this
					%end
					%if (~isinf(obj.usum(2:end)))
						%constraints = [constraints, obj.VDom'*(u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k}-273.15)) + h0) <= ly_usum(k,2:end)' - (obj.Ctu(k,:)*obj.VDom)'];
					%end
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
			obj.mpc_ctrl = optimizer(constraints,objective,ops,{x{1},w,ot,h0v,ly_uref,ly_usum},{u{1}, x{2}});
			obj.mpc_ctrl
			
		end %lin_mpc_setup
		function uout = call_mpc(obj, x, w, ot, h0v, uref, usum)
			
			%{
			if nargin < 5 || ...
				error("Needs current x, u0, w, ot");
			end
			if nargin < 7 || ...
				usum = obj.usum;
			end
			if nargin < 8 || ...
				uref = obj.uref;
			end
			
			if size(usum,1) < obj.Nhzn
				usum(size(usum,1)+1:obj.Nhzn,:) = repmat(usum(size(usum,1),:), obj.Nhzn - size(usum,1),1);
			end
			%}
			
			[uout, problem,~,~,optimizer_object] = obj.mpc_ctrl({x, w, ot, h0v, uref, usum});

			% Analyze error flags
			%if problem
			%	warning(yalmiperror(problem))
			%end			
		end %lin_mpc
	end

	%TODO don't know if these go here!!! or inside MPC!
	methods
		function [obj] = init_fnc(obj, hpc_class, Nsim)

			disc = obj.discreate_system(hpc_class, obj.Ts_ctrl);
			obj.Ad_ctrl = disc.A;
			obj.Bd_ctrl = disc.B;
			
			obj.C2 = eye(hpc_class.Ns);
			for k=2:hpc_class.full_model_layers
				obj.C2([k:hpc_class.full_model_layers:end]) = 0;
			end
			obj.C2([end-hpc_class.add_states+1:end], :) = 0;

			% Voltage
			%[obj.psac(1),obj.psac(2),obj.psac(3)] = hpc_class.pws_ls_approx([0.5 1.2], [20 90], 0.9, 1/3.497, 1.93, 1);
			[obj.psac(1),obj.psac(2),obj.psac(3)] = hpc_class.pws_ls_approx([0.5 1.2], [20 90], 0.9, [6.659 -1.979], -1.48, 1);
			% Freq
			%[obj.psac(1),obj.psac(2),obj.psac(3)] = hpc_class.pws_ls_approx([hpc_class.F_min hpc_class.F_max], [20 90], 0.9, 0.2995, 0, 0);
			%TODO parametrize
			[obj.psoff_lut, Fv, Tv] = hpc_class.pws_ls_offset(8, 8, 10);
			%obj.psoff_lut=0; Fv=1; Tv=1;
			%[obj.psoff_lut, Fv, Tv] = hpc_class.pws_ls_offset(8, 16, 10);
			obj.F0v = ones(hpc_class.Nc, length(Fv))*diag(Fv);
			obj.T0v = ones(hpc_class.Nc, length(Tv))*diag(Tv);

			obj = obj.setup_mpc(hpc_class);

			obj.xlplot = zeros(Nsim+1, hpc_class.Ns);
			obj.xlplot(1,:) = hpc_class.t_init;

			obj.tmpc = zeros(Nsim+1, hpc_class.Ns);

			obj.ex_count = 0;
			obj.failed = 0;
			obj.wl = [ones(hpc_class.Nc,1) zeros(hpc_class.Nc, hpc_class.ipl -1)];
			obj.lNsim = Nsim;
			obj.prevF = hpc_class.F_min*ones(hpc_class.Nc,1);
			obj.prevV = hpc_class.V_min*ones(hpc_class.vd,1);
			%TODO
			obj.output_mpc = 1*ones(hpc_class.Nc,1);
			

		end
		function [F,V,obj] = ctrl_fnc(obj, hc, target_index, pvt, i_pwm, i_wl)

			obj.ex_count = obj.ex_count + 1;

			%pp = obj.output_mpc;

			f_ref = hc.frtrc(min(target_index, size(hc.frtrc,1)),:)';
			p_budget = hc.tot_pw_budget(min(target_index, length(hc.tot_pw_budget)));
			T = pvt{hc.PVT_T};
			process = pvt{hc.PVT_P};

			% Process Workload
			obj.wl = obj.wl*(1-obj.alpha_wl) + i_wl*obj.alpha_wl;
			Ceff = obj.wl * (hc.dyn_ceff_k)';

			%TODO here hc.Cc
			[~, Tidx] = min(abs(obj.T0v - hc.Cc*T),[],2);
			% before overwriting F
			[~, Fidx] = min(abs(obj.F0v - obj.prevF),[],2);
			
			h0v = zeros(hc.Nc, 1);
			for i=1:hc.Nc
				atidx = min(Tidx(i)+0, 9);
				h0v(i,1) = obj.psoff_lut(atidx, Fidx(i));
			end
			
			% Choose Voltage
			domain_paired = 0;
			if (domain_paired)
				FD = diag(f_ref)*hc.VDom;
				V = obj.compute_sharedV(hc, FD, obj.voltage_rule);
			else
				fd = diag(f_ref)*ones(length(f_ref), 15);
				V = hc.FV_table(sum(fd' > hc.FV_table(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue
			end
			F = f_ref;

			if (domain_paired)
				input_mpc = F .* (hc.VDom*(V .* V));
			else
				input_mpc = F .* (V .* V);
			end

			% TODO observer
			Tobs = T;

			% TODO: here I could use the "real" formula A*x + B*u, where u
			%	is given by the identified model (so Pd+Ps*exp() with 
			%	identified parameters). But since I don't have an
			%	identification atm, I will use this to not use the real
			%	model which seems to me as too simplicistic
			power_mpc = obj.output_mpc.*Ceff + ...
				(hc.VDom*obj.prevV*hc.leak_vdd_k + hc.leak_process_k*process).* ...
				exp(hc.leak_exp_vdd_k*hc.VDom*obj.prevV + (hc.Cc*Tobs-273.15)*hc.leak_exp_t_k + hc.leak_exp_k);
			
			%state_MPC = (obj.Ad_ctrl + obj.C2*obj.psac(2))*Tobs-(obj.C2*obj.psac(2)*273.15*ones(hc.Ns,1)) + ...
			%	(obj.Bd_ctrl(:,1:hc.Nc)+(hc.Cc'*obj.psac(3)))*(obj.output_mpc.*Ceff)+obj.Bd_ctrl(:,hc.Nc+1:end)*(1000*hc.temp_amb) + hc.Cc'*(obj.psac(1) + h0v);

			state_MPC = obj.Ad_ctrl*Tobs + obj.Bd_ctrl*[power_mpc; 1000*hc.temp_amb];

			%TBD is it needed?
			%Update T0v
			%TODO here hc.Cc
			[~, Tidx] = min(abs(obj.T0v - hc.Cc*state_MPC),[],2);
			for i=1:hc.Nc
				%Here it does not change a lot if I add stuff. Fixeing it to 9 does not
				%solve the overhead problem.
				atidx = min(Tidx(i)+0, 9);
				h0v(i,1) = obj.psoff_lut(atidx, Fidx(i));
			end
			
			%state_MPC = Tobs;

			% MPC
			obj.xlplot(obj.ex_count+1,:) = 0;
			mpc_pw_target = repmat(p_budget, obj.Nhzn,1+hc.vd);
			res = obj.call_mpc(state_MPC, Ceff, hc.temp_amb*1000, h0v, input_mpc, mpc_pw_target);
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

			% Bisection
			V =[];
			for vi=1:hc.vd
				fp = f_ref.*hc.VDom(:,vi);
				cidx = [1:1:hc.Nc]' .* (fp>0);
				cidx = cidx(cidx>0);
				lim_sup = fp(fp>0);		
				lim_inf = hc.F_min*ones(sum(hc.VDom(:,vi)),1);
				in_bis = obj.output_mpc.*hc.VDom(:,vi);	
				in_bis = in_bis(in_bis>0);

				Fc = obj.bisection(lim_inf, lim_sup, @obj.intl_bis_fnc, in_bis, ...
							16, 1e-3*length(in_bis), hc.F_discretization_step*length(in_bis));

				%discretize
				if hc.F_discretization_step > 0
					Fc = round(Fc/hc.F_discretization_step) * hc.F_discretization_step;
				end

				% output
				F(cidx) = Fc;

				% TODO: Reverse Voltage choice, here again I use max.
				V(vi,1) = hc.FV_table(sum(max(Fc) > hc.FV_table(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue
			end

			%disp(obj.ex_count)
			%toto = [hc.Cc*Tobs - 273.15, hc.Cc*state_MPC- 273.15, hc.Cc*res{2} - 273.15, pp, obj.output_mpc, obj.prevF, F, hc.VDom*V,  i_wl*(hc.dyn_ceff_k)'];
			%toto = ["T0", "T1", "T2-mpc", "prevMPC T0-T1", "mpc T1-T2", "prev F T0-T1", "F T1-T2", "V T1-T2", "wl T-1 - T0"; toto]

			obj.prevF = F;
			obj.prevV = V;

		end

		function [obj] = cleanup_fnc(obj, hpc_class)
		end
		function [obj] = plot_fnc(obj, hc, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop)
			disp(strcat('[CTRL][MPC] number of times the optimization algorithm failed: ',int2str(obj.failed), '/', int2str(obj.lNsim)));
				
			figure();
			%plot(obj.Ts*sim_mul*[1:Nsim]', pceff(2:end,:)*20+293, 'g'); hold on;
			plot(t1, cpxplot(2:end,:) - 273.15, 'b'); hold on; grid on;
			plot(t2, [NaN*ones(1,size(obj.tmpc,2)); (obj.tmpc(2:end-1,:) - 273.15)], 'm'); hold on;
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

