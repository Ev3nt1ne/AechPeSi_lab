classdef controller
	%CONTROLLERS Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		Ts_ctrl (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 5e-4;			% Controller Ts
		Ts_obs (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 1e-3;			% Observer Ts
		Ts_input (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 1e-3;			% Commands min Ts
		
		Obs_poles = [0.8 0.1];		% Poles of the Luemberg Observer

		core_crit_temp = 358.15;

		

		%% All Controllers
		tot_pw_budget;
		quad_pw_budget;
		min_pw_red;
		
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
		function obj = controller(inputArg1,inputArg2)
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

