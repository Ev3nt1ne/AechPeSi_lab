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
		
		dummy_pw = 0;
		
	end
	
	methods
		function obj = CP()
			%CP Construct an instance of this class
			%   Detailed explanation goes here
			
		end
		
	end
		

	methods
		function [voltage_choice] = cp_voltage_choice(obj, hpc_class, Ft)
			for v=1:hpc_class.vd
				extrV = sum( Ft(:,v) > (hpc_class.FV_table(:,3) + [zeros(hpc_class.FV_levels-1,1); inf])', 2);
				vote_cast(:,v) = extrV(nonzeros(hpc_class.VDom(:,v).*[1:hpc_class.Nc]')) + 1;
			end
			if size(vote_cast,1) == 1
				%problem: matrix become array and prctile does not work anymore
				vote_cast(2,:) = vote_cast;
			end
			voltage_choice = hpc_class.FV_table(round(prctile(vote_cast,obj.voltage_rule)),1);
		end
		function [pu, pw_storage] = cp_pw_dispatcher(obj, T, core_crit_temp, pid_target, delta_p, ipu, min_pw_red)
			
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
				tt = (pid_target - T - safety)>=0;
				tt_value = core_crit_temp - pid_target;
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
					%Alpha = 1./(pid_target - T).*tt + 0.99.*((tt - 1) * (-1));

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
				tt = (pid_target - T - safety)>=0;
				tt_value = core_crit_temp - pid_target;
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
					%Alpha = 1./(pid_target - T).*tt + 0.99.*((tt - 1) * (-1));

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

