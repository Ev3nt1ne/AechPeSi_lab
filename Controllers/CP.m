classdef CP < controller
	%CP Summary of this class goes here
	%   Detailed explanation goes here
	
	properties

		% PID Controller
		kp = 1.9518;				% PID Kp
		ki = 73.931;				% PID Ki
		aw_up = 0.05;				% PID Anti-Windup Upper Saturation
		aw_down_c = -0.75;			% PID Anti-Windup Down Sat. Coefficient
		T_margin = 5;				% PID Temperature Margin
		voltage_rule = 100;			% Percentile in the Choice of the Voltage in the domain

		sat_up = 0;					% PID Upper Saturation
		
		pid_e_down = 2.5;			% PID banding, down constraint (absolute value)
		pid_e_up = 1;				% PID banding, up constraint (absolute value)
		pid_e_band_coeff = 0.7;		% PID banding coefficient
		
		ctrl_MA = 1;				% If CP should use Moving Average Voltage Selection
		ctrl_fixedv = 0;			% if CP should use Fixed Voltage = V_Max
		ctrl_fixedv_value = 1.2;	% the value of the fixed Voltage
		iterative_fv = 0;			% iterative fv solution. Generally: = ~ctrl_MA
		pw_inverse = 0;				% use Inverse pw reduction wrt temp
		
		alpha_wl = 0.4;				% Moving Average filter parameter for Workload
		alpha_ma = 0.1;
		
		dummy_pw = 0;

		o_inc_steps; %Tper+3; %TODO:
		o_dec_steps = 20;
		o_hys_steps = 5; %15		
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

		nOut;
		o_inc_st = 0;
		o_dec_st = 0;
		o_inc_step_count = 0; %v_inc_steps;
		o_dec_step_count = 0; %v_dec_steps;
		o_hys_act = 0;
		o_hys_step_count = 0; %v_hys_steps;
		hys_p_count = 0;
	end
	
	methods
		function obj = CP()
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

			obj.f_ma = hc.F_min*ones(hc.Nc,1);
			obj.pid_integ = zeros(hc.Nc,1);

			obj.o_inc_st = 0;
			obj.o_dec_st = 0;
			obj.o_inc_step_count = 0; %v_inc_steps;
			obj.o_dec_step_count = 0; %v_dec_steps;
			obj.o_inc_steps = 1+3;
			obj.o_hys_act = 0;
			obj.o_hys_step_count = 0; %v_hys_steps;
			obj.hys_p_count = 0;
			obj.nOut = hc.F_min*ones(hc.Nc,1);
		end
		function [F,V, obj] = ctrl_fnc(obj, hc, target_index, pvt, i_pwm, i_wl)

			obj.ex_count = obj.ex_count + 1;

			f_ref = hc.frtrc(min(target_index, size(hc.frtrc,1)),:)';
			p_budget = hc.tot_pw_budget(min(target_index, length(hc.tot_pw_budget)));

			if p_budget~=obj.pbold
				obj.pbold = p_budget;
				obj.pbc = 1;
				obj.hys_p_count = 0;
				obj.o_hys_act = 0;
				obj.o_inc_st = 0;
				obj.o_dec_st = 0;
				obj.o_inc_step_count = 0; %v_inc_steps;
				obj.o_dec_step_count = 0; %v_dec_steps;
				obj.o_hys_step_count = 0;
			else
				obj.pbc = 0;
			end

			T = pvt{hc.PVT_T};
			process = pvt{hc.PVT_P};

			obj.pw_storage = 0;

			% Process Workload
			obj.wl = obj.wl*(1-obj.alpha_wl) + i_wl*obj.alpha_wl;
			Ceff = obj.wl * (hc.dyn_ceff_k)';

			% Process Power Budget
			% (?)
	
			% Adapt Measured&Computed Power
			obj.pw_adapt = obj.cp_pw_adapt(obj.pw_adapt, i_pwm, obj.pw_old{1}, obj.pbc);
			obj.pw_old{1} = obj.pw_old{2};

			% Choose Voltage
			obj.f_ma(obj.f_ma<0) = 0; %saturate F_MA
			FD = diag(f_ref - obj.f_ma)*hc.VDom;			
			V = obj.compute_sharedV(hc, FD, obj.voltage_rule);
			F = f_ref;

			% Compute Power
			pu = (Ceff.*F.*(hc.VDom*V) + hc.leak_vdd_k).*(hc.VDom*V) + process*hc.leak_process_k;

			% DispatchPower
			%TODO, quad power budget dispatching
			delta_p = sum(pu) - p_budget + obj.pw_adapt;
			if (delta_p > 0)
				[pu, ~] = obj.cp_pw_dispatcher(T, hc.core_limit_temp, obj.T_target, delta_p, pu, hc.min_pw_red); %todo: core_limit_temp(1:obj.vd)
			end

			if sum(isempty(pu)) || sum(pu<=0)
				disp("error, something wrong");
				F = hc.F_min * ones(hc.Nc,1);
				V = hc.V_min * ones(hc.vd,1);
				%TODO should I trigger some other cleanup?
				return;
			end

			% PID
			[upid, obj.pid_integ] = obj.cp_pid(hc, T, obj.T_target, obj.pid_integ, pu);
			% Controller Output
			pu = pu + upid;
			%cpdebug(s+1,:) = upid;

			% Compute Freq
			F = (((pu - process*hc.leak_process_k) ./ (hc.VDom*V)) - hc.leak_vdd_k) ./ (hc.VDom*V) ./ Ceff;
		
			% Save OG freq
			F_og = F;


			%% HYSTERESIS F-FILTER ACTIVATION
			%{
			F_inc = (F > (obj.nOut+1e-3));
			F_dec = (F < (obj.nOut-1e-3));
			
			obj.o_inc_st = obj.o_inc_st | F_inc;
			obj.o_inc_step_count = (obj.o_inc_steps.*F_inc) + ((~F_inc).*(obj.o_inc_step_count - obj.o_inc_st));
			
			obj.o_inc_st = (obj.o_inc_step_count>0); % | o_dec_st;
			
			obj.o_dec_st = obj.o_inc_st & ( obj.o_dec_st | F_dec);
			%if o_dec_st
			%	s=s;
			%end
			%o_dec_step_count = (o_dec_steps.*V_dec) + ((~V_dec).*(o_dec_step_count - o_dec_st));
			obj.o_dec_step_count = ((obj.o_dec_step_count>0).*(obj.o_dec_step_count-1)) + (obj.o_dec_step_count<=0).*(obj.o_dec_st.*obj.o_dec_steps);
			
			obj.o_dec_st = (obj.o_dec_step_count>0);
			
			o_hys_act_hold = obj.o_hys_act;
			
			obj.o_hys_act = obj.o_hys_act | (obj.o_dec_st & F_inc);
			
			%Initialize
			obj.o_hys_step_count = ((obj.o_hys_step_count>0).*obj.o_hys_step_count) + (obj.o_hys_step_count<=0).*(obj.o_hys_act.*obj.o_hys_steps);
			%Reduce
			%pull = (Ceff.*F.*(obj.VDom*V) + obj.leak_vdd/1000).*(obj.VDom*V) + d_p*obj.leak_process/1000;
			f_hys_red = (i_pwm - p_budget) < -(p_budget*(0.125+obj.hys_p_count/10));
			%(sum(pull) - p_budget + pw_adapt) < -6 ;
			%(delta_p < (p_budget - 6)); %TODO: 6 depends on the #of cores, #of domains, and the ratio between the two.
			%TODO: also depends on the power budget itself, and the PMAX, PMIN
			obj.o_hys_step_count = (obj.o_hys_act & obj.o_hys_step_count) .* ...
				(obj.o_hys_step_count - f_hys_red);
			
			obj.o_hys_act = (obj.o_hys_step_count>0);
			
			obj.hys_p_count = obj.hys_p_count + (o_hys_act_hold ~= obj.o_hys_act) .* obj.o_hys_act;

			F = F.*((~obj.o_hys_act)|(F<=obj.nOut)) + obj.nOut.*(obj.o_hys_act & (F>obj.nOut));
			
			obj.nOut = F;


			% Save OG freq
			F_og = F;
			%}

			%%
	
			% Process Freq
			%	Check vs maxF, Temp hysteresis, etc.
			fmaxi = hc.FV_table(sum(hc.VDom * V > ones(hc.Nc,1)*hc.FV_table(:,1)'+1e-6, 2) + 1, 3);
			F = F + (F>fmaxi).*(fmaxi - F);
			F = F + (F<hc.F_min).*(hc.F_min*ones(hc.Nc,1) - F);
			
			%TEST: TODO REMOVE
			pu = (Ceff.*F.*(hc.VDom*V) + hc.leak_vdd_k).*(hc.VDom*V) + process*hc.leak_process_k;
			obj.pw_old{2} = sum(pu);

			%TODO
			alpha_MA = 0.037;
			%if powerbudget has changed || freq changed
			% interpolate a parabola
			%	parabola depends on the changes (powe >> freq, and Delta of changes)
			% f_MA = f_ref - f_app --> sat 0
			%	Sat 0 is a choice, without we have TurboBoost!
			%TODO Turbo boosting is not working!
			turbo_boost = zeros(hc.Nc,1);
			tt = (f_ref - F_og)>-turbo_boost;
			obj.f_ma = obj.f_ma*(1-alpha_MA) + alpha_MA*(tt.*(f_ref - F_og) + (~tt).*0); %This is 0 and not turbo_boost!

		end
		function [obj] = cleanup_fnc(obj, hc)
		end
		function [obj] = plot_fnc(obj, hc, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop)
		end
		
	end
		

	methods
		function [upid, ointeg] = cp_pid(obj, hc, T, pid_target, pid_integ, pu)
			
			kp_l = obj.kp;
			ki_l = obj.ki*obj.Ts_ctrl;

			% PID error
			e = pid_target - T;

			% PID error banding
			eA = ~(obj.pid_e_down>0).*(e>0) + ~(obj.pid_e_up>0).*(e<=0);
			if (obj.pid_e_down>0)
				eA = (e*(1-obj.pid_e_band_coeff)/obj.pid_e_down + obj.pid_e_band_coeff) .* (e>0);
			end
			if (obj.pid_e_up>0)
				eA = eA + ((-e*(1-obj.pid_e_band_coeff)/obj.pid_e_up + obj.pid_e_band_coeff) .* (e<=0));
			end
			%eA = eA*obj.pid_e_band_coeff + ~eA;	
			eA = eA + (eA>1).*(1-eA);
			e = e .* eA;

			% PID Integral
			ointeg = pid_integ + ki_l*e;
			eI = ointeg>obj.aw_up;
			ointeg = eI*obj.aw_up + (~eI).*ointeg;
			eI = ointeg<(pu*obj.aw_down_c);
			ointeg = eI.*(pu*obj.aw_down_c) + (~eI).*ointeg;

			% PID Proporitonal (& Output)
			upid = kp_l*e + ointeg;

			% PID Output Saturation
			eU = upid>obj.sat_up;
			upid = eU*obj.sat_up + (~eU).*upid;
			%TODO
			eU = upid < (-(pu-ones(hc.Nc,1).*hc.min_pw_red*0.7));
			upid = eU.*(-(pu-ones(hc.Nc,1).*hc.min_pw_red*0.7)) + (~eU).*upid;
		end
		function [pu, pw_storage] = cp_pw_dispatcher(obj, T, core_crit_temp, core_limit_temp, delta_p, ipu, min_pw_red)
			
			safety_pw = -1e-3;
			
			if obj.dummy_pw
				perc = delta_p / sum(ipu);
				red = ipu*perc;
				pu = ipu - red;
				
				%check on pw
				%since the difference is negligible (0.7369 vs 0.8503), but also when
				%	leakage, I will just use MAX_CEFF.
				%min_pw_red = computed above
				tt = (pu <= min_pw_red*0.7);
				pu = pu-(pu.*tt) + min_pw_red*0.7.*tt;
				pw_storage = delta_p - sum(pu);
			elseif obj.pw_inverse
				%compute alpha
				%check number
				safety = 1;
				tt = (core_limit_temp - T - safety)>=0;
				tt_value = core_crit_temp - core_limit_temp;
				if tt_value < safety
					tt_value = safety1; %0.1; %TODO (considering that a margin of 10°C is a max -otherwise too much performance loss-)
				end
				pu = ipu;
				pw_storage = -delta_p;
				while(pw_storage < safety_pw) %-5/obj.vd)
					%1
					Alpha = (core_crit_temp - T).*tt + tt_value.*(~tt);
					
					%2 
					% this seems worse
					%Alpha = 1./(core_limit_temp - T).*tt + 0.99.*((tt - 1) * (-1));

					%test on numbers .....

					%normalize
					aa = sum(Alpha);
					Alpha = Alpha/aa;	

					%test on numbers .....

					%check on power
					%moved below application
					%delta_power*Alpha > pu-obj.core_min_power
					%alpha = (pu-obj.core_min_power)/delta_p	

					%apply
					pu = pu + Alpha*pw_storage;
					%check on pw
					%since the difference is negligible (0.7369 vs 0.8503), but also when
					%	leakage, I will just use MAX_CEFF.
					%min_pw_red = computed above
					tt2 = (pu <= min_pw_red*0.7);
					pu = pu-(pu.*tt2) + min_pw_red*0.7.*tt2;
					pw_storage = (sum(ipu) - sum(pu)) - delta_p;
					if sum(tt2)==length(ipu)
						break;
					end
				end
			else
				%compute alpha
				%check number
				safety = 5;
				tt = (core_limit_temp - T - safety)>=0;
				tt_value = core_crit_temp - core_limit_temp;
				tt2 = zeros(length(ipu),1);
				if tt_value < safety
					tt_value = safety; %0.1; %TODO (considering that a margin of 10°C is a max -otherwise too much performance loss-)
				end
				pu = ipu;
				pw_storage = -delta_p;
				while(pw_storage < safety_pw) %-5/obj.vd)
					%1
					Alpha = 1./(core_crit_temp - T).*tt + 1./tt_value.*(~tt);

					%2 
					% this seems worse
					%Alpha = 1./(core_limit_temp - T).*tt + 0.99.*((tt - 1) * (-1));

					%test on numbers .....

					%normalize
					Alpha = Alpha.*(~tt2);
					aa = sum(Alpha);
					Alpha = Alpha/aa;	

					%test on numbers .....

					%check on power
					%moved below application
					%delta_power*Alpha > pu-obj.core_min_power
					%alpha = (pu-obj.core_min_power)/delta_p	

					%apply
					pu = pu + Alpha*pw_storage;
					%check on pw
					%since the difference is negligible (0.7369 vs 0.8503), but also when
					%	leakage, I will just use MAX_CEFF.
					%min_pw_red = computed above
					tt2 = (pu <= min_pw_red*0.7);
					pu = pu-(pu.*tt2) + min_pw_red*0.7.*tt2;
					pw_storage = (sum(ipu) - sum(pu)) - delta_p;
					if sum(tt2)==length(ipu)
						break;
					end
				end
			end		
		
		end

	end

	methods(Static)
		function [pw_adapt] = cp_pw_adapt(prev_pw_adapt, pw_m, pw_est, pbc, a1, a2)
			if nargin < 6
				alpha_up = 0.04;
			else
				alpha_up = a1;
			end
			if nargin < 7
				alpha_down = 0.04; %0.08
			else
				alpha_down = a2;
			end

			delta = pw_m-pw_est;
			
			alpha_pw = (delta>0)*alpha_down + (delta<=0)*alpha_up;
			
			%TODO: Here I need to do the thing that connect the alpha with
			%a poly function when power budget change. atm power budget is
			%constant
			
			pw_adapt = prev_pw_adapt*(1-alpha_pw) + delta*alpha_pw;	
			
			%here? above?
			if pbc
				pw_adapt = delta;
			end
		end
	end
end

