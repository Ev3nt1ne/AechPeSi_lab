classdef controller
	%CONTROLLERS Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		Ts_ctrl (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 5e-4;			% Controller Ts
		T_target;
		
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
		[obj] = cleanup_fnc(obj, hpc_class)
		[obj] = plot_fnc(obj, hpc_class)
	end

	methods(Static)
		function [r] = pw_function(f, pu, Ci, ks1, ks2)
			vdd_alpha = 0.3095; %0.2995
			vdd_offset = 0.07;
			%softmax_alpha = 10;
			%maxfv = sum(f.*exp(softmax_alpha*f)) / sum(exp(softmax_alpha*f));
			
			r = ( (max(f)*vdd_alpha + vdd_offset)^2*Ci.*f + ks1*(max(f)*vdd_alpha+vdd_offset)*ones(length(f),1) + ks2) - pu;
			
		end
	end
end

