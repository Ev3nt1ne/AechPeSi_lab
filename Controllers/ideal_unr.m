classdef ideal_unr < controller
	%IDEAL_UNR Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		%
	end

	properties(SetAccess=protected, GetAccess=public)
		widx;
	end
	
	methods
		function obj = ideal_unr()
			%IDEAL_UNR Construct an instance of this class
			%   Detailed explanation goes here
			
		end
		
		function [obj] = init_fnc(obj, hc, Nsim)
			obj.widx = 0;
		end
		function [F,V, obj] = ctrl_fnc(obj, hc, target_index, pvt, i_pwm, i_wl)

			obj.ex_count = obj.ex_count + 1;

			F = hc.frtrc(min(target_index, size(hc.frtrc,1)),:)';
			FD = diag(F)*hc.VDom;			
			V = obj.compute_sharedV(hc, FD, 100);
		end
		function [obj] = cleanup_fnc(obj, hc)
			obj.widx = hc.wl_index;
		end
		function [obj] = plot_fnc(obj, hc, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop)
		end

	end
end



