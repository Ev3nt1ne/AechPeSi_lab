classdef black_wolf < mpc_hpc & CP
	%UNTITLED Summary of this class goes here
	%   Detailed explanation goes here
	properties
		%alpha_wl = 0.4;				% Moving Average filter parameter for Workload

        %RE;
        Fdiscretization_step;
        lwl_mem;
        Q2 = 0;
        not_update_lin = 0;

	end
	properties(SetAccess=protected, GetAccess=public)	
		psac;
		psoff_lut;
		V0v;
		T0v;
        Pm_max;
        Pm_min;

		prevF;
		prevV;
		C2;

        Pw_diff;
	end
	
	methods
		function obj = black_wolf(hpc)
			%UNTITLED Construct an instance of this class
			%   Detailed explanation goes here
            obj = obj@CP(hpc);
		end		
	end

	methods
		function obj = setup_mpc(obj)
			
			% It's good practice to start by clearing YALMIPs internal database 
			% Every time you call sdpvar etc, an internal database grows larger
			yalmip('clear')

			x = sdpvar(repmat(obj.lNs,1,obj.Nhzn+1),repmat(1,1,obj.Nhzn+1));
			u = sdpvar(repmat(obj.lNc,1,obj.Nhzn),repmat(1,1,obj.Nhzn));
			
			ot = sdpvar(1,1); %TODO: decide if single or with horizon
			w = sdpvar(obj.lNc,1); %TODO: decide if single or with horizon
			h0v = sdpvar(obj.lNc,1);
            wSm = sdpvar(obj.lNc,1);

			%sdpvar ly_xref;
			ly_uref = sdpvar(obj.lNc,1);
			%sdpvar ly_yref;
			ly_usum = sdpvar(obj.Nhzn,obj.lvd+1);
			%ly_Cty = dpvar(repmat(obj.Nout,1,obj.Nhzn),repmat(1,1,obj.Nhzn));
			
			h0 = obj.psac(1);
			h1 = obj.psac(2);
			h2 = obj.psac(3);

			if length(h0) < obj.lNc
				h0 = ones(obj.lNc,1)*h0;
			end

			Bu = obj.Bd_ctrl(:,1:obj.lNc);
			Bd = obj.Bd_ctrl(:,obj.lNc+1:end);

            Toff = obj.Tmpc_off*ones(obj.lNs,1);
					
			constraints = [];
			%constraints = [x{2} == (Adl_mpc + C2*h1)*x{1}+(Bu+(obj.Cc'*h2))*(u0*w)+Bd*ot + obj.Cc*(h0 + h0v)];

			for k = 1 : obj.Nhzn
					%constraints = [constraints, x{k+1} == (obj.Ad_ctrl + obj.C2*h1)*x{k}-(obj.C2*h1*273.15*ones(obj.lNs,1))+(Bu+(obj.Cc'*h2))*(u{k}.*w)+Bd*ot + obj.Cc'*(h0 + h0v)];
					constraints = [constraints, x{k+1} == (obj.Ad_ctrl + Bu*diag(h1)*obj.Cc)*(x{k}) ... %-(Bu*diag(h1*273.15)*obj.Cc*ones(obj.lNs,1)) ...
									+ Bu*(u{k}.*(w+h2)) + Bd*ot + Bu*h0 + obj.Cc'*h0v ...
                                    + obj.Ad_ctrl*(-Toff) + Toff];

					%TODO: Ct I did only for y, and only for max
					if (~isinf(obj.umin))
						constraints = [constraints, u{k}.*(w+h2) + h1*(obj.Cc(1:obj.lNc,:)*(x{k})) + h0 + h0v >= obj.umin];
					end
					if (~isinf(obj.umax))
						constraints = [constraints, u{k}.*(w+h2) + h1*(obj.Cc(1:obj.lNc,:)*(x{k})) + h0 + h0v <= (obj.umax - obj.Ctu(k,:)')];
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
						constraints = [constraints, obj.Cc*x{k+1} <= (obj.T_target + obj.Tmpc_off - obj.Cty(k,:)')];
					%end
					%if (~isinf(obj.usum(1)))
                        Pw = u{k}.*(w+h2) + h1*(obj.Cc*(x{k})) + h0v+h0;
						constraints = [constraints,sum(Pw) <= (ly_usum(k,1) - sum(obj.Ctu(k,:),2))];
						%TODO should I also add the "resulting" thing. i.e.
						%	the above formula with x{k+1}?
						%	maybe no, because I use horizon for this
					%end
					%if (~isinf(obj.usum(2:end)))
						%constraints = [constraints, obj.VDom'*(u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k})) + h0) <= ly_usum(k,2:end)' - (obj.Ctu(k,:)*obj.VDom)'];
					%end
            end
			
			objective = 0;

			for k = 1:obj.Nhzn %TODO T_target(1) not good!!!!!
				objective = objective + (((u{k}-ly_uref)-obj.Pm_min)/(obj.Pm_max-obj.Pm_min))'*obj.Rt*(((u{k}-ly_uref)-obj.Pm_min)/(obj.Pm_max-obj.Pm_min)) + ...
                        ((u{k}-obj.Pm_min)/(obj.Pm_max-obj.Pm_min))'*obj.Rs*((u{k}-obj.Pm_min)/(obj.Pm_max-obj.Pm_min)) + ...
                        ((x{k+1}-obj.T_amb-Toff)/(obj.T_target(1)-obj.T_amb))'*obj.Q*((x{k+1}-obj.T_amb-Toff)/(obj.T_target(1)-obj.T_amb)) + ...
                        ((u{k}-obj.Pm_min)/(obj.Pm_max-obj.Pm_min).*wSm)'*obj.R*((u{k}-obj.Pm_min)/(obj.Pm_max-obj.Pm_min).*wSm) + ...
                        (x{k+1}-obj.T_amb)'*(-obj.Q2*1)*(x{k+1}-obj.T_amb) + (-obj.Cc'*(h0+h0v)/h1)'*(obj.Q2*1)*(-obj.Cc'*(h0+h0v)/h1);
                %energy
                %at = (ly_uref./u{k}).*wSm;
                %Ev = [Pw; at; h1*(obj.Cc*(x{k+1}-x{k})); (at-1)];
                %objective = objective + Ev'*obj.RE*Ev;
			end
			
			%ops = sdpsettings('verbose',1,'solver','quadprog', 'usex0',1);
			ops = sdpsettings('verbose',1,'solver','osqp', 'usex0',0); %You have specified an initial point, but the selected solver (OSQP) does not support warm-starts through YALMIP
			ops.quadprog.TolCon = 1e-2;
			ops.quadprog.TolX = 1e-5;
			ops.quadprog.TolFun = 1e-3;
			ops.convertconvexquad = 0;
			%ops.quadprog.MaxPCGIter = max(1, ops.quadprog.MaxPCGIter * 3);
			ops.quadprog.MaxIter = 200;%50;
			obj.mpc_ctrl = optimizer(constraints,objective,ops,{x{1},w,ot,h0v,wSm,ly_uref,ly_usum},{u{1}, x{2}});
			obj.mpc_ctrl
			
		end %lin_mpc_setup
		function uout = call_mpc(obj, x, w, ot, h0v, wSm, uref, usum)
			
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
			
			[uout, problem,~,~,~, sol] = obj.mpc_ctrl({x, w, ot, h0v, wSm, uref, usum});

            obj.ctrl_info = [obj.ctrl_info; sol.solvertime];

			% Analyze error flags
			%if problem
			%	warning(yalmiperror(problem))
			%end			
		end %lin_mpc
	end

	%TODO don't know if these go here!!! or inside MPC!
	methods
		function [obj, comms] = init_fnc(obj, hc, chip, ctrl_id, Nsim)

			obj.initialize(chip, Nsim);

            obj.ctrl_info = [];

			% Voltage
			%[obj.psac(1),obj.psac(2),obj.psac(3)] = hc.pws_ls_approx([0.5 1.2], [20 90], 0.9, 1/3.497, 1.93, 1);
            if obj.not_update_lin 
               obj.not_update_lin = 0;
            else
			    [obj.psac(1),obj.psac(2),obj.psac(3)] = hc.pws_ls_approx(chip,[0.5 1.2], [20+273.15 90+273.15], obj.Tmpc_off, 0.9, [6.659 -1.979], -1.48, 1);
			    % Freq
			    %[obj.psac(1),obj.psac(2),obj.psac(3)] = hc.pws_ls_approx([obj.lFmin obj.lFmax], [20 90], 0.9, 0.2995, 0, 0);
			    %TODO parametrize
                show = 0;
			    [obj.psoff_lut, Fv, Tv] = hc.pws_ls_offset(chip, obj, 5, 4, hc.temp_amb, chip.core_crit_temp, show);

                obj.V0v = ones(obj.lNc, length(Fv))*diag(Fv);
			    obj.T0v = ones(obj.lNc, length(Tv))*diag(Tv);
            end


            %todo?
            obj.Pm_max = obj.lFmax*obj.lVmax^2;		
            obj.Pm_min = obj.lFmin*obj.lVmin^2;

            obj.umin = chip.core_min_power; %obj.Pm_min;
            obj.umax = chip.core_max_power; %obj.Pm_max;

            %TODO: check if all the matrix are ok and defined, otherwise
            %       fix them

			obj = obj.setup_mpc();

			obj.xlplot = zeros(Nsim+1, obj.lNs);
			obj.xlplot(1,:) = obj.T_amb*ones(obj.lNs,1);

			obj.tmpc = zeros(Nsim+1, obj.lNs);
            obj.Pw_diff = zeros(Nsim+1, 1);

			obj.ex_count = 0;
            obj.hyst_idx = -1;
			obj.failed = 0;
			obj.wl = [ones(obj.lNc,1) zeros(obj.lNc, obj.lipl -1)];
			obj.lNsim = Nsim;
			obj.prevF = obj.lFmin*ones(obj.lNc,1);
			obj.prevV = obj.lVmin*ones(obj.lvd,1);
            %rounding initial guess
            rig = 5;
            obj.Tobs = fix(hc.t_init{ctrl_id}/rig)*rig; %hc.T_amb*ones(obj.lNs,1);
			%TODO
			obj.output_mpc = 1*ones(obj.lNc,1);

            obj.hyst_duration = 4;

            obj.pw_old = 1*obj.lNc;

            % Observer
            % create the dicrete system
            %sys = ss(obj.Ad_ctrl, obj.Bd_ctrl, obj.C, zeros(size(obj.C,1), ...
            %    size(obj.Bd_ctrl,2)), obj.Ts_ctrl);
            %To use kalman, you must provide a model sys that has an input 
            %   for the noise w. Thus, sys is not the same as Plant, because 
            %   Plant takes the input un = u + w. You can construct sys by 
            %   creating a summing junction for the noise input.
            %sys = sys*[1 1];
            %For this example, assume both noise sources have unit covariance 
            %   and are not correlated (N = 0).
            %[kalmf,L,P] = kalman(sys,Q,R,N);
            [obj.LK, P, Z] = dlqe(obj.Ad_ctrl, obj.Gw, obj.C, obj.Qcov, obj.Rcov);

            par = zeros(3,1);
            par(1) = 0.2;
            par(2) = obj.pw_old;
            par(3) = hc.toto_pw_budget(1);

            [obj, comms] = obj.init_grad_track_comm(hc, ctrl_id, par);
		end
		function [F,V, comm, obj] = ctrl_fnc(obj, f_ref, pwbdg, pvt, i_pwm, i_wl, ctrl_id, ctrl_comm)

			obj.ex_count = obj.ex_count + 1;

            %if (obj.ex_count == 300)
            %    aa = 1;
            %end

			%pp = obj.output_mpc;
			T = pvt{obj.PVT_T};
            Tc = T(1:obj.lNc);
			process = pvt{obj.PVT_P};

            chippwbdg = pwbdg(2);
            totpwbdg = pwbdg(1);

            %%%%%%%
            % Distributed Algorithm:
            par = zeros(3,1);
            par(1) = sum(i_wl*(1:obj.lipl)')/obj.lNc;
            %TODO improve this
            par(2) = (obj.pw_old+chippwbdg)/2;
            par(3) = totpwbdg;
            [dist_pw, dist_grad, obj] = obj.grad_track_alg(ctrl_comm, ctrl_id,par);

            comm{1} = dist_pw;
            obj.gta_pl_x(obj.ex_count+1, :) = dist_pw';
            comm{2} = dist_grad;
            obj.gta_pl_grad(obj.ex_count+1, :) = dist_grad';

            chippwbdg = min(chippwbdg, dist_pw(ctrl_id));
            %%%%%%%

            % Process Workload
			obj.wl = obj.wl*(1-obj.alpha_wl) + i_wl*obj.alpha_wl;
			Ceff = obj.wl * (obj.pw_ceff)';

			%TODO here obj.Cc
			%[~, Tidx] = min(abs(obj.T0v - obj.Cc*T),[],2);
            Tidx = sum(obj.T0v < Tc,2)+1;
			% before overwriting F
			%[~, Vidx] = min(abs(obj.V0v - obj.prevV),[],2);
            Vidx = sum(obj.V0v < obj.lVDom*obj.prevV,2)+1;
			
			h0v = zeros(obj.lNc, 1);
			for i=1:obj.lNc
				atidx = min(Tidx(i)+0, 9);
				h0v(i,1) = obj.psoff_lut(atidx, Vidx(i));
			end
			
			% Choose Voltage
			domain_paired = 0;
			if (domain_paired)
				FD = diag(f_ref)*obj.lVDom;
				V = obj.find_dom_sharedV(obj.lFVT, FD, obj.voltage_rule);
            else
                %TODO: 15
				fd = diag(f_ref)*ones(length(f_ref), size(obj.lFVT,1));
				V = obj.F2VM(obj.lFVT, fd);
                %V = obj.lFVT(sum(fd' > obj.lFVT(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue
			end
			F = f_ref;

			if (domain_paired)
				input_mpc = F .* (obj.lVDom*(V .* V));
			else
				input_mpc = F .* (V .* V);
			end

			% TODO observer
			%below?

			% TODO: here I could use the "real" formula A*x + B*u, where u
			%	is given by the identified model (so Pd+Ps*exp() with 
			%	identified parameters). But since I don't have an
			%	identification atm, I will use this to not use the real
			%	model which seems to me as too simplicistic
			power_mpc = obj.output_mpc.*Ceff + ...
				(obj.lVDom*obj.prevV*obj.pw_stat_lin(1) + obj.pw_stat_lin(3)*process).* ...
				exp(obj.pw_stat_exp(1)*obj.lVDom*obj.prevV + (Tc)*obj.pw_stat_exp(2) + obj.pw_stat_exp(3));
			
			%state_MPC = (obj.Ad_ctrl + obj.C2*obj.psac(2))*obj.Tobs + ...
			%	(obj.Bd_ctrl(:,1:obj.lNc)+(obj.Cc'*obj.psac(3)))*(obj.output_mpc.*Ceff)+obj.Bd_ctrl(:,obj.lNc+1:end)*(1000*obj.T_amb) + obj.Cc'*(obj.psac(1) + h0v);

			%state_MPC = obj.Ad_ctrl*obj.Tobs + obj.Bd_ctrl*[power_mpc; 1000*obj.T_amb];
            obj.Tobs = (obj.Ad_ctrl - obj.LK*obj.C)*obj.Tobs + obj.Bd_ctrl*[power_mpc; 1000*obj.T_amb] + obj.LK*T;
            state_MPC = obj.Tobs;

			%TBD is it needed?
			%Update T0v
			%TODO here obj.Cc
			%[~, Tidx] = min(abs(obj.T0v - obj.Cc*state_MPC),[],2);
            Tidx = sum(obj.T0v < obj.Cc*state_MPC,2)+1;
			for i=1:obj.lNc
				%Here it does not change a lot if I add stuff. Fixeing it to 9 does not
				%solve the overhead problem.
				atidx = min(Tidx(i)+0, 9);
				h0v(i,1) = obj.psoff_lut(atidx, Vidx(i));
			end
			
			%state_MPC = obj.Tobs;

			% MPC
			obj.xlplot(obj.ex_count+1,:) = 0;
			mpc_pw_target = repmat(chippwbdg, obj.Nhzn,1+obj.lvd);
            wSm = obj.wl*(1-obj.lwl_mem)';
			res = obj.call_mpc(state_MPC+obj.Tmpc_off, Ceff, (obj.T_amb)*1000, h0v, wSm, input_mpc, mpc_pw_target);
			tt = isnan(res{1});
			% too complex to make it vectorial
			%dp = res{1}(~tt);
			%obj.output_mpc = repmat(tt, length(tt), length(dp))*dp + ...tt.*obj.output_mpc;
			if (~tt)
                if (obj.hyst_idx < obj.ex_count)
				    obj.output_mpc = res{1};
                end
			%TODO
            else
                obj.output_mpc = obj.Pm_min; %obj.output_mpc/3*2;
                obj.hyst_idx = obj.ex_count + obj.hyst_duration;
			end
			obj.tmpc(obj.ex_count+1, :) = res{2}';

			%TODO outmpc_min
			outmpc_min = 0.1;
			obj.output_mpc(obj.output_mpc<outmpc_min) = outmpc_min;
			obj.failed = obj.failed + ((sum(tt)>0) | (sum(obj.output_mpc<outmpc_min)>0));

			% Bisection
			V =[];
			for vi=1:obj.lvd
				fp = f_ref.*obj.lVDom(:,vi);
				cidx = [1:1:obj.lNc]' .* (fp>0);
				cidx = cidx(cidx>0);
				lim_sup = fp(fp>0);		
				lim_inf = obj.lFmin*ones(sum(obj.lVDom(:,vi)),1);
				in_bis = obj.output_mpc.*obj.lVDom(:,vi);	
				in_bis = in_bis(in_bis>0);

				Fc = obj.bisection(lim_inf, lim_sup, @obj.intl_bis_fnc, in_bis, ...
							16, 1e-3*length(in_bis), obj.Fdiscretization_step*length(in_bis));

				%discretize
				if obj.Fdiscretization_step > 0
					Fc = round(Fc/obj.Fdiscretization_step) * obj.Fdiscretization_step;
				end

				% output
				F(cidx) = Fc;

				% TODO: Reverse Voltage choice, here again I use max.
				V(vi,1) = obj.F2VM(obj.lFVT, max(Fc)); %+1e-6 to fix matlab issue
			end

			%disp(obj.ex_count)
			%toto = [obj.Cc*obj.Tobs - 273.15, obj.Cc*state_MPC- 273.15, obj.Cc*res{2} - 273.15, pp, obj.output_mpc, obj.prevF, F, obj.lVDom*V,  i_wl*(hc.dyn_ceff_k)'];
			%toto = ["T0", "T1", "T2-mpc", "prevMPC T0-T1", "mpc T1-T2", "prev F T0-T1", "F T1-T2", "V T1-T2", "wl T-1 - T0"; toto]

            obj.pw_old = sum((F.*(obj.lVDom*V).*(obj.lVDom*V)).*Ceff + ...
				(obj.lVDom*V*obj.pw_stat_lin(1) + obj.pw_stat_lin(3)*process).* ...
				exp(obj.pw_stat_exp(1)*obj.lVDom*V + (Tc)*obj.pw_stat_exp(2) + obj.pw_stat_exp(3)));

			obj.prevF = F;
			obj.prevV = V;

		end

		function [obj] = cleanup_fnc(obj,hc)
		end
		function [obj] = plot_fnc(obj, t1, t2, xop, uop, fop, vop, wlop)
			disp(strcat('[CTRL][MPC] number of times the optimization algorithm failed: ',int2str(obj.failed), '/', int2str(obj.lNsim)));
				
			figure();
			%plot(obj.Ts*sim_mul*[1:Nsim]', pceff(2:end,:)*20+293, 'g'); hold on;
			%plot(t1, cpxplot(2:end,:) - 273.15, 'b'); hold on; grid on;
			%plot(t2, [NaN*ones(1,size(obj.tmpc,2)); (obj.tmpc(2:end-1,:) - 273.15)], 'm'); hold on;
			plot(t1, [(xop(2+2:end,:) - 273.15); NaN*ones(2,size(obj.tmpc,2))], 'b'); hold on;
			plot(t2, obj.tmpc(2:end,:) - 273.15 -obj.Tmpc_off, 'm');
			xlabel("Time [s]");
			ylabel("Temperature [T]");

            figure();
			plot(t2, obj.gta_pl_grad(2:end,:), 'b');
			xlabel("Time [s]");
			ylabel("Gradient of GTA");
		end
	end

	methods
		function res = intl_bis_fnc(obj, in, args)
			res = obj.pw_function(in, args, 1, 0, 0);
		end
	end



end

