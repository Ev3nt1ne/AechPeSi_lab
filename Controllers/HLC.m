classdef HLC < controller
	%CP Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
        hp = 100e-3;
	
	end

	properties(SetAccess=protected, GetAccess=public)	
		pw_storage = 0;
		pw_adapt = 0;
		pw_old;
		wl;
		pbold = 0;
		pbc = 0;

		f_ma;
		pid_integ;

		Adctrl;
		Bdctrl;

        Fp;
        Vp;

	end
	
	methods
		function obj = HLC()
			%CP Construct an instance of this class
			%   Detailed explanation goes here
			
		end
		function [obj] = init_fnc(obj, hc, Nsim)
			%TODO understand which one I actually really need
			%TODO understand which needs to go out and which in
			obj.pw_storage = 0;
			obj.pw_adapt = 0;
			obj.pw_old = [];
			obj.pw_old{1} = 1*hc.Nc;
			obj.pw_old{2} = 1*hc.Nc;
			obj.pbold = 0;
			obj.pbc = 0;
			obj.ex_count = 0;
			obj.wl = [ones(hc.Nc,1) zeros(hc.Nc, hc.ipl -1)];
			obj.T_target = ones(hc.Nc, 1)*hc.core_limit_temp;

			disc = obj.discreate_system(hc, obj.Ts_ctrl);
			obj.Adctrl = disc.A;
			obj.Bdctrl= disc.B;

            obj.Fp = hc.F_max*ones(hc.Nc,1);
            obj.Vp = 1.1*ones(hc.vd,1);
		end
		function [F,V, obj] = ctrl_fnc(obj, hc, target_index, pvt, i_pwm, i_wl)

			obj.ex_count = obj.ex_count + 1;

			f_ref = hc.frtrc(min(target_index, size(hc.frtrc,1)),:)';
			p_budget = hc.tot_pw_budget(min(target_index, length(hc.tot_pw_budget)));

			if p_budget~=obj.pbold
				obj.pbold = p_budget;
				obj.pbc = 1;
			else
				obj.pbc = 0;
			end

			T = pvt{hc.PVT_T};
			process = pvt{hc.PVT_P};

			obj.pw_storage = 0;

            fhp = round(obj.hp / obj.Ts_ctrl);
            fms = round(1e-3/obj.Ts_ctrl);

            if (mod(obj.ex_count, fhp) == 0) || ...
                    (sum(hc.Cc*T >= obj.T_target) && (obj.ex_count > fms))

            obj.ex_count = 0;

			% Process Workload
			%obj.wl = obj.wl*(1-obj.alpha_wl) + i_wl*obj.alpha_wl;
			%Ceff = iwl * (hc.dyn_ceff_k)';

			% Process Power Budget
			% (?)
	
			% Adapt Measured&Computed Power
			%obj.pw_adapt = obj.cp_pw_adapt(obj.pw_adapt, i_pwm, obj.pw_old{1}, obj.pbc);
			%obj.pw_old{1} = obj.pw_old{2};

			% Voltage
			FD = diag(f_ref)*hc.VDom;			
			V = obj.compute_sharedV(hc, FD, 100);
			F = f_ref;

			% Compute Power
			pu = hc.power_compute(F,hc.VDom*V,hc.Cc*T,i_wl,process);

			% initial test
			xi = T;
			nj = 5;
			for i=1:nj
				xi = obj.Adctrl*xi + obj.Bdctrl*[pu; hc.temp_amb*1000];
			end

			dfj = hc.F_discretization_step;%*4;
			utsens = i_wl * (1 - hc.wl_mem_weigth');

			%starting the real control
			if sum(hc.Cc*xi >= obj.T_target)
				xi = obj.Adctrl*T;
				cc = hc.Cc*(obj.Bdctrl(:,1:hc.Nc)*pu) > (obj.T_target - hc.Cc*xi);
				pun = pu;
				it = 1;
				fl = F <= hc.F_min + 1e-4;
				while sum(cc)

					deltaO = - utsens * (dfj) ./ (F .* (F - dfj));
					deltaO(~cc) = 1e10;
					%deltaO(fl) = 1e10;
					[~ , pos] = min(abs(deltaO));
					F(pos) = F(pos) - dfj;
					FD = diag(F)*hc.VDom;			
					V = obj.compute_sharedV(hc, FD, 100);
					Vn = hc.VDom*V;
					Tn = hc.Cc*xi;

					pun(pos) = hc.power_compute(F(pos),Vn(pos),Tn(pos),i_wl(pos,:),process(pos));
					%xi = T;
					%for i=1:nj
					%	xi = obj.Adctrl*xi + obj.Bdctrl*[pu; hc.temp_amb*1000];
					%end
					xi = obj.Adctrl*T;
					cc = hc.Cc*(obj.Bdctrl(:,1:hc.Nc)*pun) > (obj.T_target - hc.Cc*xi);

					fl = F <= hc.F_min + 1e-4;
					F(fl) = hc.F_min;
					% check this 
					cc = cc & ~fl;
					it = it + 1;
					if mod(it, 100) == 1
						aa = 1;
					end
				end
				%it;
            end
            obj.Fp = F;
            obj.Vp = V;
            else
                F = obj.Fp;
                V = obj.Vp;
            end
			
		end
		function [obj] = cleanup_fnc(obj, hc)
		end
		function [obj] = plot_fnc(obj, hc, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop)
		end
		
	end
		
end

