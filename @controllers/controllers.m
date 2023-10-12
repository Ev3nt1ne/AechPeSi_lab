classdef controllers
	%CONTROLLERS Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		Ts_ctrl (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 5e-4;			% Controller Ts
		Ts_obs (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 1e-3;			% Observer Ts
		Ts_mpc (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 5e-3;			% MPC Ts
		Ts_input (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 1e-3;			% Commands min Ts
		Obs_poles = [0.8 0.1];		% Poles of the Luemberg Observer

		tsim = 1;					% Simulation Time in [s]
		tesim = 2;					% Free Simulation Time in [s]

		x_init;						% Initial Conditions 
		urplot;						% Input Reference Plot
		frplot;						% Freq Reference Plot
		zrplot;						% Power Noise Plot
		wrplot;						% Wokrload Plot

		freq_fact = 8;				% Times per seconds expected Target Freq Changes

		%% All Controllers
		tot_pw_budget;
		quad_pw_budget;
		
		%% MPC Controller
		xref;
		uref;
		yref;
		usum;
		Q;
		R;
		R2;
		Nhzn = 3;
		Ctx;
		Ctu;
		Cty;
		
		mpc_robustness = 1;
		
		umin;
		uMax;
		xmin;
		xMax;
		ymin;
		yMax;
		
		%% SoA/PID Controller
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
		
		core_crit_temp = 358.15;
	end

	properties(Dependent)
		Ni_nl;						% Number of inputs for the non-linear case
		Ni_c;						% Number of controllable inputs
		Ni_nc;						% Number of non-controllable inputs

		Ad_mpc;
		Bd_mpc;
		Ad_ctrl;
		Bd_ctrl;
		Ad_obs;
		Bd_obs;
	end
	
	methods
		function obj = controllers(inputArg1,inputArg2)
			%CONTROLLERS Construct an instance of this class
			%   Detailed explanation goes here
			%obj.Property1 = inputArg1 + inputArg2;
		end
		
		function outputArg = method1(obj,inputArg)
			%METHOD1 Summary of this method goes here
			%   Detailed explanation goes here
			%outputArg = obj.Property1 + inputArg;
		end
	end
end

