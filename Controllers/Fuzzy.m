classdef Fuzzy < CP
	%FUZZY Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		%global var
		pw_storage = 0;
		pw_adapt = 0;
		pw_old;
		pw_ad_step_count = 0;
		wl;
		VCred;
		Vn;
		v_inc_st = 0;
		v_dec_st = 0;
		v_inc_step_count = 0; %v_inc_steps;
		v_dec_step_count = 0; %v_dec_steps;
		v_hys_act = 0;
		v_hys_step_count = 0; %v_hys_steps;
		hys_p_count = 0;
		pbold = 0;
		pbc = 0;
		ex_count = 0;
		derT = 0;
		T_old = 0;

		%parameters
		Tper=8;
		pw_ad_aup = 0.015; %0.05; %0.025;
		pw_ad_adown = 0.02; %0.075; %0.025; %0.1;
		pw_ad_achange = 0.5; %0.5;
		pw_ad_steps = 20;
		v_inc_steps; %Tper+3; %TODO:
		v_dec_steps = 20;
		v_hys_steps = 5; %15
	end
	
	methods
		function obj = Fuzzy()
			%FUZZY Construct an instance of this class
			%   Detailed explanation goes here
			

			% 
		end

		function [obj] = init_fnc(obj, hc)
			%TODO understand which one I actually really need
			obj.v_inc_st = 0;
			obj.v_dec_st = 0;
			obj.pw_storage = 0;
			obj.pw_adapt = 0;
			obj.pw_old = [];
			obj.pw_old{1} = 1*hc.Nc;
			obj.pw_old{2} = 1*hc.Nc;
			obj.pw_ad_step_count = 0;
			obj.v_inc_step_count = 0; %v_inc_steps;
			obj.v_dec_step_count = 0; %v_dec_steps;
			obj.v_inc_steps = obj.Tper+3;
			obj.v_hys_act = 0;
			obj.v_hys_step_count = 0; %v_hys_steps;
			obj.hys_p_count = 0;
			obj.pbold = 0;
			obj.pbc = 0;
			obj.ex_count = 0;
			obj.derT = zeros(hc.Nc,1);
			obj.T_old = hc.x_init(1);
			obj.wl = [ones(hc.Nc,1) zeros(hc.Nc, hc.ipl -1)];
			obj.VCred = zeros(hc.Nc,1);
			obj.Vn = hc.V_min*ones(hc.vd,1);
			obj.T_target = ones(hc.Nc, 1)*hc.core_crit_temp;
		end		
		function [F,V, obj] = ctrl_fnc(obj, hc, target_index, pvt, i_pwm, i_wl)

			obj.ex_count = obj.ex_count + 1;

			f_ref = hc.frtrc(min(target_index, size(hc.frtrc,1)),:)';
			p_budget = hc.tot_pw_budget(min(target_index, length(hc.tot_pw_budget)));
			if p_budget~=obj.pbold
				obj.pbold = p_budget;
				obj.pbc = 1;
				obj.hys_p_count = 0;
				obj.v_hys_act = 0;
				obj.v_inc_st = 0;
				obj.v_dec_st = 0;
				obj.v_inc_step_count = 0; %v_inc_steps;
				obj.v_dec_step_count = 0; %v_dec_steps;
				obj.v_hys_step_count = 0;
				obj.pw_ad_step_count = obj.pw_ad_steps;
			else
				obj.pbc = 0;
				obj.pw_ad_step_count = obj.pw_ad_step_count - 1;
				obj.pw_ad_step_count = (obj.pw_ad_step_count<0)*0 + (obj.pw_ad_step_count>=0)*obj.pw_ad_step_count;
			end

			T = pvt{hc.PVT_T};
			process = pvt{hc.PVT_P};
			
			if (obj.ex_count>obj.Tper+1) && (mod(obj.ex_count,obj.Tper)==2)
				obj.derT = T - obj.T_old;
				obj.T_old = T;
			end

			obj.pw_storage = 0;

			% Process Workload
			obj.wl = obj.wl*(1-obj.alpha_wl) + i_wl*obj.alpha_wl;
			Ceff = obj.wl * (hc.dyn_ceff_k)';

			% Process Power Budget
			% (?)
	
			% Adapt Measured&Computed Power
			aup = obj.pw_ad_aup + obj.pw_ad_step_count * (obj.pw_ad_achange - obj.pw_ad_aup) / obj.pw_ad_steps;
			adown = obj.pw_ad_adown + obj.pw_ad_step_count * (obj.pw_ad_achange - obj.pw_ad_adown) / obj.pw_ad_steps;
			
			obj.pw_adapt = obj.cp_pw_adapt(obj.pw_adapt, i_pwm, obj.pw_old{1}, (obj.pbc | (obj.pw_ad_step_count>(obj.pw_ad_steps-5))), aup, adown);
			obj.pw_old{1} = obj.pw_old{2};

			% Choose Voltage
			FD = diag(f_ref)*hc.VDom;			
			V = obj.cp_voltage_choice(hc, FD);
			F = f_ref;

			% Control Temperature:
			if (mod(obj.ex_count,obj.Tper)==2)
				Tband = [45.0 65.0 80.0 85 ];
				Tband = Tband + 273.15;
				derBand = [0 0.5 1.0 2.0];
				derBand = derBand + hc.sensor_noise*hc.sensor_noise_amplitude(hc.PVT_T)/2;
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
			
			VDred = min(diag(obj.VCred)*hc.VDom);
			V = V + fix(VDred/2)'*0.05;

			%cap F:
			% maybe remove?
			fmaxi = hc.FV_table(sum(hc.VDom * V > ones(hc.Nc,1)*hc.FV_table(:,1)', 2) + 1, 3);
			F = F + (F>fmaxi).*(fmaxi - F);
			F = F + hc.VDom*(mod(VDred,2).*sign(VDred))'*0.05;

			% Compute Power
			pu = (Ceff.*F.*(hc.VDom*V) + hc.leak_vdd_k).*(hc.VDom*V) + process*hc.leak_process_k;
			
			% DispatchPower
			%TODO, quad power budget dispatching
			delta_p = sum(pu) - p_budget + obj.pw_adapt;
			if (delta_p > 0)
				dompw = hc.VDom'*pu;
				domT = ((hc.VDom'*T ./ sum(hc.VDom)') + max(T.*hc.VDom)' ) / 2; %here mean*2 / 3??
				%TODO this T_target conversion is wrong!
				[resdpw, ~] = obj.cp_pw_dispatcher(domT, hc.core_crit_temp, obj.T_target(1:hc.vd), delta_p, dompw, hc.min_pw_red); %todo: pid_target(1:obj.vd)
				deltapd = dompw-resdpw;
				for vdi=1:hc.vd
					pd = pu.*hc.VDom(:,vdi);
					pd = pd(pd>0);
					cd = Ceff.*hc.VDom(:,vdi);
					cd = cd(cd>0);
					cidx = [1:1:hc.Nc]' .* hc.VDom(:,vdi);
					cidx = cidx(cidx>0);
					[pd, ~] = obj.cp_pw_dispatcher_c(cd,deltapd(vdi), pd, hc.min_pw_red);
					pu(cidx) = pd;
				end
				%pw_storage = pw_storage + pws;	
			end

			if sum(isempty(pu)) || sum(pu<=0)
				disp("error, something wrong");
				F = hc.F_min * ones(hc.Nc,1);
				V = hc.V_min * ones(hc.vd,1);
				%TODO should I trigger some other cleanup?
				return;
			end

			% Compute Freq
			%% Newton-Rapson		
			for vi=1:hc.vd
				fp = f_ref.*hc.VDom(:,vi);	
				cidx = [1:1:hc.Nc]' .* (fp>0);
				cidx = cidx(cidx>0);
				lim_sup = fp(fp>0);		
				lim_inf = hc.F_min*ones(sum(hc.VDom(:,vi)),1);
				pup = pu.*hc.VDom(:,vi);
				pup = pup(pup>0);
				Cip = Ceff.*hc.VDom(:,vi);
				Cip = Cip(Cip>0);
				d_pk = (process*hc.leak_process_k).*hc.VDom(:,vi);
				d_pk = d_pk(d_pk>0);
				xi = obj.pw_function(lim_inf, pup, Cip, hc.leak_vdd_k,d_pk);			
				xs = obj.pw_function(lim_sup, pup, Cip, hc.leak_vdd_k,d_pk);
				
				fo = zeros(length(cidx),1);
				
				%boundaries:
				sn = xs>0;
				so = xi>0;
				c_lim_inf = lim_inf + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);
				c_lim_sup = lim_sup + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_sup);
				%fo = fo + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);
				
				for nri=1:16
					fo = (c_lim_inf+c_lim_sup)/2;
					xs = obj.pw_function(fo, pup, Cip, hc.leak_vdd_k,d_pk);
					if (sum(abs(xs)) <= 1e-3*length(cidx)) || (sum(abs(c_lim_sup-fo))<= hc.F_discretization_step*length(cidx))
						break;
					end
					sn = xs>0;
					so = xi>0;
					
					xi = xi + (sn==so).*(xs-xi);
					c_lim_inf = c_lim_inf + (sn==so).*(fo-c_lim_inf);
					c_lim_sup = c_lim_sup + (sn~=so).*(fo-c_lim_sup);				
				end
				
				%discretize
				if hc.F_discretization_step > 0
					fo = round(fo/hc.F_discretization_step) * hc.F_discretization_step;
				end
				F(cidx) = fo;
				V(vi) = hc.FV_table(sum(max(fo) > hc.FV_table(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue
			end	

			%% HYSTERESIS V-FILTER ACTIVATION
			V_inc = (V > (obj.Vn+1e-3));
			V_dec = (V < (obj.Vn-1e-3));
			
			obj.v_inc_st = obj.v_inc_st | V_inc;
			obj.v_inc_step_count = (obj.v_inc_steps.*V_inc) + ((~V_inc).*(obj.v_inc_step_count - obj.v_inc_st));
			
			obj.v_inc_st = (obj.v_inc_step_count>0); % | v_dec_st;
			
			obj.v_dec_st = obj.v_inc_st & ( obj.v_dec_st | V_dec);
			%if v_dec_st
			%	s=s;
			%end
			%v_dec_step_count = (v_dec_steps.*V_dec) + ((~V_dec).*(v_dec_step_count - v_dec_st));
			obj.v_dec_step_count = ((obj.v_dec_step_count>0).*(obj.v_dec_step_count-1)) + (obj.v_dec_step_count<=0).*(obj.v_dec_st.*obj.v_dec_steps);
			
			obj.v_dec_st = (obj.v_dec_step_count>0);
			
			v_hys_act_hold = obj.v_hys_act;
			
			obj.v_hys_act = obj.v_hys_act | (obj.v_dec_st & V_inc);
			
			%Initialize
			obj.v_hys_step_count = ((obj.v_hys_step_count>0).*obj.v_hys_step_count) + (obj.v_hys_step_count<=0).*(obj.v_hys_act.*obj.v_hys_steps);
			%Reduce
			%pull = (Ceff.*F.*(obj.VDom*V) + obj.leak_vdd/1000).*(obj.VDom*V) + d_p*obj.leak_process/1000;
			v_hys_red = (i_pwm - p_budget) < -(p_budget*(0.125+obj.hys_p_count/10));
			%(sum(pull) - p_budget + pw_adapt) < -6 ;
			%(delta_p < (p_budget - 6)); %TODO: 6 depends on the #of cores, #of domains, and the ratio between the two.
			%TODO: also depends on the power budget itself, and the PMAX, PMIN
			obj.v_hys_step_count = (obj.v_hys_act & obj.v_hys_step_count) .* ...
				(obj.v_hys_step_count - v_hys_red);
			
			obj.v_hys_act = (obj.v_hys_step_count>0);
			
			obj.hys_p_count = obj.hys_p_count + (v_hys_act_hold ~= obj.v_hys_act) .* obj.v_hys_act;
			
			%{
			ppap(s,[1:obj.vd]) = v_inc_st;
			ppap(s,[obj.vd+1:obj.vd*2]) = v_inc_step_count;
			ppap(s,[obj.vd*2+1:obj.vd*3]) = v_dec_st;
			ppap(s,[obj.vd*3+1:obj.vd*4]) = v_dec_step_count;
			ppap(s,[obj.vd*4+1:obj.vd*5]) = v_hys_act;
			ppap(s,[obj.vd*5+1:obj.vd*6]) = v_hys_step_count;
			ppap(s,[obj.vd*6+1:obj.vd*7]) = V_inc;
			ppap(s,[obj.vd*7+1:obj.vd*8]) = V_dec;
			%}		
	
			V = V.*((~obj.v_hys_act)|(V<=obj.Vn)) + obj.Vn.*(obj.v_hys_act & (V>obj.Vn));
			
			obj.Vn = V;

			%pwap(s) = pw_adapt;

			% Save OG freq
			%F_og = F;
	
			% Process Freq
			%	Check vs maxF, Temp hysteresis, etc.
			fmaxi = hc.FV_table(sum(hc.VDom * V > ones(hc.Nc,1)*hc.FV_table(:,1)'+1e-6, 2) + 1, 3);
			F = F + (F>fmaxi).*(fmaxi - F);
			F = F + (F<hc.F_min).*(hc.F_min*ones(hc.Nc,1) - F);
			
			%TEST: TODO REMOVE
			pu = (Ceff.*F.*(hc.VDom*V) + hc.leak_vdd_k).*(hc.VDom*V) + process*hc.leak_process_k;
			obj.pw_old{2} = sum(pu);


		end
		function [obj] = cleanup_fnc(obj, hc)
		end
		function [obj] = plot_fnc(obj, hc)
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

