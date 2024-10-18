classdef Fuzzy < CP
	%FUZZY Summary of this class goes here
	%   Detailed explanation goes here
		
	properties
		%parameters
		Tper=8;
		pw_ad_aup = 0.015; %0.05; %0.025;
		pw_ad_adown = 0.02; %0.075; %0.025; %0.1;
		pw_ad_achange = 0.5; %0.5;
		pw_ad_steps = 20;
		%{
		o_inc_steps; %Tper+3; %TODO:
		o_dec_steps = 20;
		o_hys_steps = 5; %15
        %}
        Tsensor_amplitude;
        Fdiscretization_step;
	end

	properties(SetAccess=protected, GetAccess=public)	
		%global var
		pw_ad_step_count = 0;
		VCred;
		%{
		nOut;
		o_inc_st = 0;
		o_dec_st = 0;
		o_inc_step_count = 0; %o_inc_steps;
		o_dec_step_count = 0; %o_dec_steps;
		o_hys_act = 0;
		o_hys_step_count = 0; %o_hys_steps;
		hys_p_count = 0;
		%}
		derT = 0;
		T_old = 0;
        T_count = 0;
	end
	
	methods
		function obj = Fuzzy(hpc)
			%FUZZY Construct an instance of this class
			%   Detailed explanation goes here
			
         %%% Object Initialization %%
         % Call superclass constructor before accessing object
         % You cannot conditionalize this statement
         obj = obj@CP(hpc);
			% 
		end

		function [obj, comms] = init_fnc(obj, hc, chip, ctrl_id, Nsim)
			%TODO understand which one I actually really need
			%TODO understand which needs to go out and which in
            obj.initialize(chip, Nsim);

			obj.o_inc_st = 0;
			obj.o_dec_st = 0;
			obj.pw_storage = 0;
			obj.pw_adapt = 0;
			obj.pw_old = [];
			obj.pw_old{1} = 1*obj.lNc;
			obj.pw_old{2} = 1*obj.lNc;
			obj.pw_ad_step_count = 0;
			obj.o_inc_step_count = 0; %o_inc_steps;
			obj.o_dec_step_count = 0; %o_dec_steps;
			obj.o_inc_steps = obj.Tper+3;
			obj.o_hys_act = 0;
			obj.o_hys_step_count = 0; %o_hys_steps;
			obj.hys_p_count = 0;
			obj.pbold = 0;
			obj.pbc = 0;
			obj.ex_count = 0;
            obj.T_count = 0;
			obj.derT = zeros(obj.lNc,1);
			obj.T_old = obj.T_amb;
			obj.wl = [ones(obj.lNc,1) zeros(obj.lNc, obj.lipl -1)];
			obj.VCred = zeros(obj.lNc,1);
			obj.nOut = obj.lVmin*ones(obj.lvd,1);

            par = zeros(3,1);
            par(1) = 0.2;
            par(2) = obj.pw_old{2};
            par(3) = hc.toto_pw_budget(1);

            [obj, comms] = obj.init_grad_track_comm(hc, ctrl_id, par);
		end		
		function [F,V, comm, obj] = ctrl_fnc(obj, f_ref, pwbdg, pvt, i_pwm, i_wl, ctrl_id, ctrl_comm)

			if obj.T_count >= 10*obj.Tper
				obj.T_count = 2*obj.Tper;
			end
			obj.T_count = obj.ex_count + 1;
            obj.ex_count = obj.ex_count + 1;
            
            chippwbdg = pwbdg(2);
            totpwbdg = pwbdg(1);

            %%%%%%%
            % Distributed Algorithm:
            par = zeros(3,1);
            par(1) = sum(i_wl*(1:obj.lipl)')/obj.lNc;
            %TODO here in fuzzy this can be moved below the Thermal
            %reduction so I can correctly put the target power without
            %considering the thermal control
            par(2) = (obj.pw_old{2}+ obj.pw_adapt+chippwbdg)/2;
            par(3) = totpwbdg;
            [dist_pw, dist_grad, obj] = obj.grad_track_alg(ctrl_comm, ctrl_id,par);

            comm{1} = dist_pw;
            obj.gta_pl_x(obj.ex_count+1, :) = dist_pw';
            comm{2} = dist_grad;
            obj.gta_pl_grad(obj.ex_count+1, :) = dist_grad';

            chippwbdg = min(chippwbdg, dist_pw(ctrl_id));
            %%%%%%%

            %TODO: here what I should do is:
            %   A) Filter chippwbdg
            %   B) compare: abs(obj.pbold-Filtered_chippwbdg) > K, (e.g. K=5)
            %   
			if pwbdg(2)~=obj.pbold
				obj.pbold = pwbdg(2);
				obj.pbc = 1;
				obj.hys_p_count = 0;
				obj.o_hys_act = 0;
				obj.o_inc_st = 0;
				obj.o_dec_st = 0;
				obj.o_inc_step_count = 0; %o_inc_steps;
				obj.o_dec_step_count = 0; %o_dec_steps;
				obj.o_hys_step_count = 0;
				obj.pw_ad_step_count = obj.pw_ad_steps;
			else
				obj.pbc = 0;
				obj.pw_ad_step_count = obj.pw_ad_step_count - 1;
				obj.pw_ad_step_count = (obj.pw_ad_step_count<0)*0 + (obj.pw_ad_step_count>=0)*obj.pw_ad_step_count;
			end

			T = pvt{obj.PVT_T};
			process = pvt{obj.PVT_P};
			
			if (mod(obj.T_count,obj.Tper)==2) % && (obj.T_count>obj.Tper+1)
				obj.derT = T - obj.T_old;
				obj.T_old = T;
			end

			obj.pw_storage = 0;

			% Process Workload
			obj.wl = obj.wl*(1-obj.alpha_wl) + i_wl*obj.alpha_wl;
			Ceff = obj.wl * (obj.pw_ceff)';

			% Process Power Budget
			% (?)
	
			% Adapt Measured&Computed Power
			aup = obj.pw_ad_aup + obj.pw_ad_step_count * (obj.pw_ad_achange - obj.pw_ad_aup) / obj.pw_ad_steps;
			adown = obj.pw_ad_adown + obj.pw_ad_step_count * (obj.pw_ad_achange - obj.pw_ad_adown) / obj.pw_ad_steps;
			
			obj.pw_adapt = obj.cp_pw_adapt(obj.pw_adapt, i_pwm, obj.pw_old{1}, (obj.pbc | (obj.pw_ad_step_count>(obj.pw_ad_steps-5))), aup, adown);
			obj.pw_old{1} = obj.pw_old{2};

			% Choose Voltage
			FD = diag(f_ref)*obj.lVDom;
			V = obj.find_dom_sharedV(obj.lFVT, FD, obj.voltage_rule);
			F = f_ref;

			% Control Temperature:
			if (mod(obj.T_count,obj.Tper)==2)
				Tband = [45.0 65.0 80.0 85 ];
				Tband = Tband + 273.15;
				derBand = [0 0.5 1.0 2.0];
				derBand = derBand + obj.Tsensor_amplitude/2;
				redmat =	[2	2	1	0	-1;
						 	2	1	1	0	-1;
						 	1	1	0	-1	-2;
						 	1	0	-1	-2	-3;
						 	0	0	-2	-3	-4];
	
				tr(:,2) = sum(T > Tband,2);
				tr(:,1) = sum(obj.derT > derBand,2);
				tr = tr + 1; %indexing
	
				obj.VCred = obj.VCred + redmat((tr(:,2)-1)*size(redmat,1) + tr(:,1));
				%saturate: (no turbo atm)
				obj.VCred = obj.VCred + (obj.VCred > 0).*(0-obj.VCred);
			end
			
			VDred = min(diag(obj.VCred)*obj.lVDom);
			%VDred = VDred + (VDred > 0).*(0-VDred); %Saturation here is
			%	way worse. Dunno Why.
			V = V + fix(VDred/2)'*0.05; %fix is better than floor!!!
			%V = V + floor(VDred/2)'*0.05;
			%saturate V
			V = V + (V < obj.lVmin).*(obj.lVmin-V);

			%cap F:
			% maybe remove?
            fmaxi = obj.V2FM(obj.lFVT, obj.lVDom * V);
			F = F + (F>fmaxi).*(fmaxi - F);
			F = F + obj.lVDom*(mod(VDred,2).*sign(VDred))'*0.05;

			% Compute Power
			pu = (Ceff.*F.*(obj.lVDom*V) + obj.pw_stat_lin(1)).*(obj.lVDom*V) + ...
                process*obj.pw_stat_lin(3);
			
			% DispatchPower
			%TODO, quad power budget dispatching
			delta_p = sum(pu) - chippwbdg + obj.pw_adapt;
			if (delta_p > 0)
				dompw = obj.lVDom'*pu;
				domT = ((obj.lVDom'*T ./ sum(obj.lVDom)') + max(T.*obj.lVDom)' ) / 2; %here mean*2 / 3??
				%TODO this T_target conversion is wrong!
				[resdpw, ~] = obj.cp_pw_dispatcher(domT, obj.core_Tcrit, obj.T_target(1:obj.lvd), delta_p, dompw, obj.lPmin); %todo: pid_target(1:obj.vd)
				deltapd = dompw-resdpw;
				for vdi=1:obj.lvd
					pd = pu.*obj.lVDom(:,vdi);
					pd = pd(pd>0);
					cd = Ceff.*obj.lVDom(:,vdi);
					cd = cd(cd>0);
					cidx = [1:1:obj.lNc]' .* obj.lVDom(:,vdi);
					cidx = cidx(cidx>0);
					[pd, ~] = obj.cp_pw_dispatcher_c(cd,deltapd(vdi), pd, obj.lPmin);
					pu(cidx) = pd;
				end
				%pw_storage = pw_storage + pws;	
			end

			if sum(isempty(pu)) || sum(pu<=0)
				disp("error, something wrong");
				F = obj.lFmin * ones(obj.lNc,1);
				V = obj.lVmin * ones(obj.lvd,1);
				%TODO should I trigger some other cleanup?
				return;
			end

			% Compute Freq
			%% Newton-Rapson		
			for vi=1:obj.lvd
				fp = f_ref.*obj.lVDom(:,vi);	
				cidx = [1:1:obj.lNc]' .* (fp>0);
				cidx = cidx(cidx>0);
				lim_sup = fp(fp>0);		
				lim_inf = obj.lFmin*ones(sum(obj.lVDom(:,vi)),1);
				pup = pu.*obj.lVDom(:,vi);
				pup = pup(pup>0);
				Cip = Ceff.*obj.lVDom(:,vi);
				Cip = Cip(Cip>0);
				d_pk = (process*obj.pw_stat_lin(3)).*obj.lVDom(:,vi);
				d_pk = d_pk(d_pk>0);
				xi = obj.pw_function(lim_inf, pup, Cip, obj.pw_stat_lin(1),d_pk);			
				xs = obj.pw_function(lim_sup, pup, Cip, obj.pw_stat_lin(1),d_pk);
				
				fo = zeros(length(cidx),1);
				
				%boundaries:
				sn = xs>0;
				so = xi>0;
				c_lim_inf = lim_inf + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);
				c_lim_sup = lim_sup + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_sup);
				%fo = fo + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);
				
				for nri=1:16
					fo = (c_lim_inf+c_lim_sup)/2;
					xs = obj.pw_function(fo, pup, Cip, obj.pw_stat_lin(1),d_pk);
					if (sum(abs(xs)) <= 1e-3*length(cidx)) || (sum(abs(c_lim_sup-fo))<= obj.Fdiscretization_step*length(cidx))
						break;
					end
					sn = xs>0;
					so = xi>0;
					
					xi = xi + (sn==so).*(xs-xi);
					c_lim_inf = c_lim_inf + (sn==so).*(fo-c_lim_inf);
					c_lim_sup = c_lim_sup + (sn~=so).*(fo-c_lim_sup);				
				end
				
				%discretize
				if obj.Fdiscretization_step > 0
					fo = round(fo/obj.Fdiscretization_step) * obj.Fdiscretization_step;
				end
				F(cidx) = fo;
				V(vi) = obj.F2VM(obj.lFVT,max(fo));
			end	

			%% HYSTERESIS V-FILTER ACTIVATION
			V_inc = (V > (obj.nOut+1e-3));
			V_dec = (V < (obj.nOut-1e-3));
			
			obj.o_inc_st = obj.o_inc_st | V_inc;
			obj.o_inc_step_count = (obj.o_inc_steps.*V_inc) + ((~V_inc).*(obj.o_inc_step_count - obj.o_inc_st));
			
			obj.o_inc_st = (obj.o_inc_step_count>0); % | o_dec_st;
			
			obj.o_dec_st = obj.o_inc_st & ( obj.o_dec_st | V_dec);
			%if o_dec_st
			%	s=s;
			%end
			%o_dec_step_count = (o_dec_steps.*V_dec) + ((~V_dec).*(o_dec_step_count - o_dec_st));
			obj.o_dec_step_count = ((obj.o_dec_step_count>0).*(obj.o_dec_step_count-1)) + (obj.o_dec_step_count<=0).*(obj.o_dec_st.*obj.o_dec_steps);
			
			obj.o_dec_st = (obj.o_dec_step_count>0);
			
			o_hys_act_hold = obj.o_hys_act;
			
			obj.o_hys_act = obj.o_hys_act | (obj.o_dec_st & V_inc);
			
			%Initialize
			obj.o_hys_step_count = ((obj.o_hys_step_count>0).*obj.o_hys_step_count) + (obj.o_hys_step_count<=0).*(obj.o_hys_act.*obj.o_hys_steps);
			%Reduce
			%pull = (Ceff.*F.*(obj.VDom*V) + obj.leak_vdd/1000).*(obj.VDom*V) + d_p*obj.leak_process/1000;
			v_hys_red = (i_pwm - chippwbdg) < -(chippwbdg*(0.125+obj.hys_p_count/10));
			%(sum(pull) - chippwbdg + pw_adapt) < -6 ;
			%(delta_p < (chippwbdg - 6)); %TODO: 6 depends on the #of cores, #of domains, and the ratio between the two.
			%TODO: also depends on the power budget itself, and the PMAX, PMIN
			obj.o_hys_step_count = (obj.o_hys_act & obj.o_hys_step_count) .* ...
				(obj.o_hys_step_count - v_hys_red);
			
			obj.o_hys_act = (obj.o_hys_step_count>0);
			
			obj.hys_p_count = obj.hys_p_count + (o_hys_act_hold ~= obj.o_hys_act) .* obj.o_hys_act;
			
			%{
			ppap(s,[1:obj.vd]) = o_inc_st;
			ppap(s,[obj.vd+1:obj.vd*2]) = o_inc_step_count;
			ppap(s,[obj.vd*2+1:obj.vd*3]) = o_dec_st;
			ppap(s,[obj.vd*3+1:obj.vd*4]) = o_dec_step_count;
			ppap(s,[obj.vd*4+1:obj.vd*5]) = o_hys_act;
			ppap(s,[obj.vd*5+1:obj.vd*6]) = o_hys_step_count;
			ppap(s,[obj.vd*6+1:obj.vd*7]) = V_inc;
			ppap(s,[obj.vd*7+1:obj.vd*8]) = V_dec;
			%}		
	
			V = V.*((~obj.o_hys_act)|(V<=obj.nOut)) + obj.nOut.*(obj.o_hys_act & (V>obj.nOut));
			
			obj.nOut = V;

			%pwap(s) = pw_adapt;

			% Save OG freq
			%F_og = F;
	
			% Process Freq
			%	Check vs maxF, Temp hysteresis, etc.
			fmaxi = obj.V2FM(obj.lFVT, obj.lVDom * V);
			F = F + (F>fmaxi).*(fmaxi - F);
			F = F + (F<obj.lFmin).*(obj.lFmin*ones(obj.lNc,1) - F);
			
			%TEST: TODO REMOVE
			pu = (Ceff.*F.*(obj.lVDom*V) + obj.pw_stat_lin(1)).*(obj.lVDom*V) + process*obj.pw_stat_lin(3);
			obj.pw_old{2} = sum(pu);

		end
		function [obj] = cleanup_fnc(obj, hc)
		end
		function [obj] = plot_fnc(obj, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop)
            figure();
			plot(t2, obj.gta_pl_grad(2:end,:), 'b');
			xlabel("Time [s]");
			ylabel("Gradient of GTA");
		end

		function [pu, pw_storage] = cp_pw_dispatcher_c(obj, Ceff, delta_p, ipu, min_pw_red)
			
			safety_pw = -1e-3;
			
			%perc = delta_p / sum(ipu);
			%red = ipu*perc;
			%pu = ipu - red;
			
			pu = ipu;
			tt = zeros(length(ipu),1);
			%pw_storage = (sum(ipu) - sum(pu)) - delta_p;
			pw_storage = -delta_p;
			%here it should be pw_storage = - delta_p; CHECK!
			while(pw_storage < safety_pw) %-5/obj.vd)
				Alpha = 1./(Ceff+(sum(Ceff)/length(Ceff)));
				Alpha = Alpha.*(~tt);
				%normalize
				aa = sum(Alpha);
				Alpha = Alpha/aa;
				%apply
				red = Alpha*pw_storage;
				red = red + (-red>Alpha.*pu).*(-Alpha.*pu - red);
				%pu = pu + Alpha*pw_storage;
				pu = pu + red;

				%check on pw
				%since the difference is negligible (0.7369 vs 0.8503), but also when
				%	leakage, I will just use MAX_CEFF.
				%min_pw_red = computed above
				tt = (pu <= min_pw_red);
				pu = pu-(pu.*tt) + min_pw_red.*tt;
				pw_storage = (sum(ipu) - sum(pu)) - delta_p;
				if sum(tt)==length(ipu)
					break;
				end
			end				
		end
	end
end

