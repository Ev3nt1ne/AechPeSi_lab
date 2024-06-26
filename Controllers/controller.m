classdef controller < handle
	%CONTROLLERS Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		Ts_ctrl (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 5e-4;			% Controller Ts
		T_target = 85 + 273.15;

		C {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1];

	end

	properties(Dependent)
		%Ad_ctrl;
		%Bd_ctrl;
	end

	properties(SetAccess=protected, GetAccess=public)	
		osunix;
		Ad_ctrl;
		Bd_ctrl;
		ex_count = 0;
		lNsim;
	end
	
	methods
		function obj = controller()
			%CONTROLLERS Construct an instance of this class
			%   Detailed explanation goes here
			%obj.Property1 = inputArg1 + inputArg2;
			if isunix
				obj.osunix = 1;
			else
				obj.osunix = 0;
			end

		end
		
	end

	methods(Abstract=true)
		[obj] = init_fnc(obj, hpc, Nsim)
		[F,V,obj] = ctrl_fnc(obj, target_index, pvt, i_pwm, i_wl)
		[obj] = cleanup_fnc(obj)
		[obj] = plot_fnc(obj, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop)
	end

	methods(Static)
		function [r] = pw_function(f, pu, Ci, ks1, ks2)
			%TODO parametrize these values
			vdd_alpha = 0.3095; %0.2995
			vdd_offset = 0.07;
			%softmax_alpha = 10;
			%maxfv = sum(f.*exp(softmax_alpha*f)) / sum(exp(softmax_alpha*f));
			
			r = ( (max(f)*vdd_alpha + vdd_offset)^2*Ci.*f + ks1*(max(f)*vdd_alpha+vdd_offset)*ones(length(f),1) + ks2) - pu;
        end
		function res = bisection(lim_inf, lim_sup, fnc, args, it, tolx, tolf)

			N = length(lim_inf);

			xi = fnc(lim_inf, args);	
			xs = fnc(lim_sup, args);

			res = zeros(N,1);

			%boundaries:
			sn = xs>0;
			so = xi>0;
			c_lim_inf = lim_inf + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);
			c_lim_sup = lim_sup + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_sup);
			%fo = fo + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);

			for nri=1:it
				res = (c_lim_inf+c_lim_sup)/2;
				xs = fnc(res, args);
				if (sum(abs(xs)) <= tolx) || (sum(abs(c_lim_sup-res))<= tolf)
					break;
				end
				sn = xs>0;
				so = xi>0;

				xi = xi + (sn==so).*(xs-xi);
				c_lim_inf = c_lim_inf + (sn==so).*(res-c_lim_inf);
				c_lim_sup = c_lim_sup + (sn~=so).*(res-c_lim_sup);				
			end
        end
        function [FMax] = V2FM(FVT, V)
            dimV = length(V);
            FMax = FVT(sum(V > ones(dimV,1)*FVT(:,1)'+1e-6, 2) + 1, 3);
        end
        function [VMax] = F2VM(FVT, F)
            dimV = length(F);
            VMax = FVT(sum(F > ones(dimV,1)*FVT(:,3)'+1e-6, 2) + 1, 1);
        end
        function [voltage_choice] = find_dom_sharedV(FVT, Ft, vrule)
            lvd = size(Ft,2);
            lNc = size(Ft,1);
            lFVlev = size(FVT,1);
            voltage_choice = -1*ones(lvd,1);
            comp_mat = (FVT(:,3) + [1e-6*ones(lFVlev-1,1); inf])';
            lvdom = Ft>0;
			for v=1:lvd
				extrV = sum( Ft(:,v) > comp_mat, 2);
				%vote_cast(:,v) = extrV(nonzeros(obj.VDom(:,v).*[1:obj.Nc]')) + 1;
				vote_cast = extrV(nonzeros(lvdom(:,v).*[1:lNc]')) + 1;
				voltage_choice(v,1) = FVT(round(prctile(vote_cast,vrule)),1);
			end
			%if size(vote_cast,1) == 1
				%problem: matrix become array and prctile does not work anymore
			%	vote_cast(2,:) = vote_cast;
			%end
			%voltage_choice = obj.FV_table(round(prctile(vote_cast,obj.voltage_rule)),1);
        end   
    end

	%% Dependent Variables
	methods
		%% SET
		%{
		function obj = set.Ts_ctrl(obj, val)
			warning("[CTRL]: Rememer to run the function ' ' to compute the discreate matrixes.")
			obj.Ts_ctrl = val;
		end
		%}
	end
end

