classdef controller
	%CONTROLLERS Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		Ts_ctrl (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 5e-4;			% Controller Ts

		
	end

	properties(Dependent)

		Ad_ctrl;
		Bd_ctrl;

	end
	
	methods
		function obj = controller()
			%CONTROLLERS Construct an instance of this class
			%   Detailed explanation goes here
			%obj.Property1 = inputArg1 + inputArg2;
		end
		
	end

	methods(Abstract=true)
		[obj] = init_fnc(obj, hpc_class)
		[F,V,obj] = ctrl_fnc(obj, hpc_class, target_index, pvt, i_pwm, i_wl)
		[] = cleanup_fnc(obj, hpc_class)
		[] = plot_fnc(obj, hpc_class)
	end
end

