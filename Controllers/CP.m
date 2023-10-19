classdef CP < controllers
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
		function obj = CP(inputArg1,inputArg2)
			%CP Construct an instance of this class
			%   Detailed explanation goes here
			obj.Property1 = inputArg1 + inputArg2;
		end
		
		function outputArg = method1(obj,inputArg)
			%METHOD1 Summary of this method goes here
			%   Detailed explanation goes here
			outputArg = obj.Property1 + inputArg;
		end
	end
end

