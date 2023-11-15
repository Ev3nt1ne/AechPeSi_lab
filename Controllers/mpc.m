classdef mpc < controller
	%MPC Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
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


		Ts_obs (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 1e-3;			% Observer Ts

		Obs_poles = [0.8 0.1];		% Poles of the Luemberg Observer


	end

	properties(Dependent)
		Ad_obs;
		Bd_obs;
	end

	properties(SetAccess=protected, GetAccess=public)		
		mpc_ctrl;					% Persistent variable to optimize controller (prev: Controller)
	end
	
	methods
		function obj = mpc()
			%MPC Construct an instance of this class
			%   Detailed explanation goes here

			%% TODO
			% check Yalmip
			% check OSQP
			
		end
	end

	methods(Abstract=true)
		[obj] = setup_mpc(obj)
		[uout] = call_mpc(obj)
	end
end

