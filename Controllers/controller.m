classdef controller < handle
	%CONTROLLERS Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		Ts_ctrl (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 5e-4;			% Controller Ts
		T_target = 85 + 273.15;

		C {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1];

        T_amb = 25 + 273.15;

        adj_mat;
        comm_mat;

        %Maybe do a class alone for this???
        multi_der_fcn;
        multi_prev_grad;
        multi_prev_x;
        multi_prev_y;
        multi_learn_rate;
	end

	properties(Dependent)
		%Ad_ctrl;
		%Bd_ctrl;
	end

	properties(SetAccess=protected, GetAccess=public)	
		osunix;

        ex_count = 0;
		lNsim;

		Ad_ctrl;
		Bd_ctrl;
        Cc;
        lNc;
        lNh;
        lNv
        lNs;
        core_Tcrit;

        lvd;
        lVDom;
        lipl;
        pw_stat_lin;
        pw_stat_exp;
        pw_dyn;
        pw_ceff;
        lPmin;
        lPmax;

        lFVT;
        lFmin;
        lFmax;
        lVmin;
        lVmax;

    end

    properties(SetAccess=immutable, GetAccess=public)
        PVT_T = 3;
        PVT_V = 2;
        PVT_P = 1;
	end
	
	methods
		function obj = controller(chip)
			%CONTROLLERS Construct an instance of this class
			%   Detailed explanation goes here
			%obj.Property1 = inputArg1 + inputArg2;
			if isunix
				obj.osunix = 1;
			else
				obj.osunix = 0;
            end

            obj = obj.initialize(chip, []);

        end

        function obj = initialize(obj, chip, Nsim)

            obj.ex_count = 0;
		    obj.lNsim = Nsim;

            [obj.Ad_ctrl, obj.Bd_ctrl, obj.Cc, ...
                obj.lNc, obj.lNh, obj.lNv, ...
                obj.core_Tcrit] = chip.thermal_give_model(obj.Ts_ctrl);
            obj.lNs = length(obj.Ad_ctrl);

            [obj.pw_stat_lin, obj.pw_stat_exp, obj.pw_dyn, obj.pw_ceff, ...
                obj.lPmin, obj.lPmax] = chip.power_give_model();
            [obj.lvd, obj.lVDom] = chip.domain_give_model();
            obj.lipl = length(obj.pw_ceff);
    
            obj.lFVT = chip.give_fvtable();
            obj.lFmin = obj.lFVT(1,2);
            obj.lFmax = obj.lFVT(end,3);
            obj.lVmin = obj.lFVT(1,1);
            obj.lVmax = obj.lFVT(end,1);
        end
        function [obj, comm] = init_grad_track_comm(obj, hc, ctrl_id)
            if isempty(obj.comm_mat)
                error("[CTRL] Communication matrix empty");
            end
            obj.adj_mat = obj.create_adj_mat(obj.comm_mat, [2,3]);
            ll = length(hc.chip_pw_budget);
            c1 = ones(ll,1);
            for i=1:ll
                c1(i) = hc.chip_pw_budget{i}(1);
            end
            comm{1} = c1;
            grad = obj.multi_der_fcn(c1,ctrl_id);
            comm{2} = grad;
            obj.multi_prev_x = c1;
            obj.multi_prev_y = grad;
            obj.multi_prev_grad = grad;
        end
        function [x, y, obj] = grad_track_alg(obj, ctrl_comm, ctrl_id)
            GAMMA = obj.multi_learn_rate;
            % Gather informations from neighboors
            u = 0;
            v = 0;
            for idx=1:length(ctrl_comm)
                if (~isempty(ctrl_comm{idx}) && (idx ~= ctrl_id))
                    val = ctrl_comm{idx}{1};
                    u = u + obj.adj_mat(ctrl_id,idx)*val;
                    val = ctrl_comm{idx}{2};
                    v = v + obj.adj_mat(ctrl_id,idx)*val;
                end
            end
            x = obj.adj_mat(ctrl_id,ctrl_id)*obj.multi_prev_x + u;            
            x = x - GAMMA*obj.multi_prev_y;
            obj.multi_prev_x = x;
            
            y = obj.adj_mat(ctrl_id,ctrl_id)*obj.multi_prev_y + v;
            grad = obj.multi_der_fcn(x,ctrl_id);
            y = y + (grad - obj.multi_prev_grad);
            obj.multi_prev_y = y;
            obj.multi_prev_grad = grad;
        end
	end

	methods(Abstract=true)
		[obj, comms] = init_fnc(obj, hc, chip, ctrl_id, Nsim)
		[F,V,comm,obj] = ctrl_fnc(obj, f_ref, pwbdg, pvt, i_pwm, i_wl,ctrl_id, ctrl_comm)
		[obj] = cleanup_fnc(obj, hc)
		[obj] = plot_fnc(obj, t1, t2, xop, uop, fop, vop, wlop)
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
            FMax = FVT(sum(V > ones(dimV,1)*FVT(:,1)'+1e-6, 2) + 1, 3); %+1e-6 to fix matlab issue
        end
        function [VMax] = F2VM(FVT, F)
            dimV = size(F,1);
            VMax = FVT(sum(F > ones(dimV,1)*FVT(:,3)'+1e-6, 2) + 1, 1); %+1e-6 to fix matlab issue
        end
        function [voltage_choice] = find_dom_sharedV(FVT, Ft, vrule)
            lvd = size(Ft,2);
            lNc = size(Ft,1);
            lFVlev = size(FVT,1);
            voltage_choice = -1*ones(lvd,1);
            comp_mat = (FVT(:,3) + [1e-6*ones(lFVlev-1,1); inf])'; %+1e-6 to fix matlab issue
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
        function [Adj] = create_adj_mat(mat, mul, max_iter, tol)
            if (nargin < 4) || isempty(tol)
                tol = 1e-4;
            end
            if (nargin < 3) || isempty(max_iter)
                max_iter = 1000;
            end
            if (nargin < 2) || isempty(mul)
                mul = [1,1];
            end
            if length(mul) ~= 2
                warning("[CTRL] Wrong mul input, setting default [1,1]");
                mul = [1,1];
            end

            diagmat = diag(mat);
            am = mat - eye(size(mat)).*diagmat;
            am = am*mul(2);
            % convert zeros in 1
            if (sum(diagmat==0)>0)
                md = min(diagmat(diagmat>0));
                if isempty(md) || (md > 1)
                    md = 1;
                end
                %diagmat(diagmat==0) = md;
                diagmat = diagmat + md;
            end
            diagmat = diagmat*mul(1);
            Adj = am + eye(size(mat)).*diagmat;

            iter = 0;
            convcheck = 0;
            while ( (iter<max_iter) && (~convcheck))
                iter = iter + 1;
                % Normalize rows
                rows_sum = sum(Adj,2);
                Adj = Adj ./ rows_sum;
        
                % Normalize columns
                cols_sum = sum(Adj,1);
                Adj = (Adj' ./ cols_sum')';
        
                % Check convergence
                convcheck = (sum(sum(Adj,1)) <= tol) && (sum(sum(Adj,2)) <= tol);
            end
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

