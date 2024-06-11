classdef dumb_test < controller
	%CP Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
        freq = 0.4;
	
	end

	properties(SetAccess=protected, GetAccess=public)	

	end
	
	methods
		function obj = dumb_test()
			%CP Construct an instance of this class
			%   Detailed explanation goes here
			
		end
		function [obj] = init_fnc(obj, hc, Nsim)
			%TODO understand which one I actually really need
			%TODO understand which needs to go out and which in
		end
		function [F,V, obj] = ctrl_fnc(obj, hc, target_index, pvt, i_pwm, i_wl)

            F = obj.freq;
            V = hc.FV_table(sum(F > hc.FV_table(:,3))+1,1);
			
		end
		function [obj] = cleanup_fnc(obj, hc)
		end
		function [obj] = plot_fnc(obj, hc, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop)
		end
		
	end
		
end

