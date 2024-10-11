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
		function obj = ideal_unr(chip)
			%IDEAL_UNR Construct an instance of this class
			%   Detailed explanation goes here
            obj = obj@controller(chip);
			
		end
		
		function [obj, comms] = init_fnc(obj, hc, chip, ctrl_id, Nsim)
			obj.widx = 0;
            comms = 0;
		end
		function [F,V,comm,obj] = ctrl_fnc(obj, f_ref, pwbdg, pvt, i_pwm, i_wl,ctrl_id, ctrl_comm)

			obj.ex_count = obj.ex_count + 1;

			F = f_ref;
			FD = diag(F)*obj.lVDom;			
			V = obj.find_dom_sharedV(obj.lFVT, FD, 100);

            comm = 0;
		end
		function [obj] = cleanup_fnc(obj, hc)
			obj.widx = hc.wl_index;
		end
		function [obj] = plot_fnc(obj, hc, t1, t2, xop, uop, fop, vop, wlop)
		end

	end
end



