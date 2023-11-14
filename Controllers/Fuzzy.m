classdef Fuzzy < CP
	%FUZZY Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		%global var
		pw_storage = 0;
		pw_adapt = 0;
		pw_old = 0;
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

		function [obj] = init_fnc(obj, hpc_class)
			%TODO understand which one I actually really need
			obj.v_inc_st = 0;
			obj.v_dec_st = 0;
			obj.pw_storage = 0;
			obj.pw_adapt = 0;
			obj.pw_old = 0;
			obj.pw_ad_step_count = 0;
			obj.v_inc_step_count = 0; %v_inc_steps;
			obj.v_dec_step_count = 0; %v_dec_steps;
			obj.v_hys_act = 0;
			obj.v_hys_step_count = 0; %v_hys_steps;
			obj.hys_p_count = 0;
		end		
		function [F,V, obj] = ctrl_fnc(obj, hpc_class, target_index, pvt, i_pwm, i_wl)

			f_ref
			pbc
			s 
			T
			derT
			p_budget


			obj.pw_storage = 0;

			% Process Workload
			obj.wl = obj.wl*(1-obj.alpha_wl) + i_wl*obj.alpha_wl;
			Ceff = obj.wl * (hpc_class.dyn_ceff_k)';

			% Process Power Budget
			% (?)
	
			% Adapt Measured&Computed Power
			aup = obj.pw_ad_aup + obj.pw_ad_step_count * (obj.pw_ad_achange - obj.pw_ad_aup) / obj.pw_ad_steps;
			adown = obj.pw_ad_adown + obj.pw_ad_step_count * (obj.pw_ad_achange - obj.pw_ad_adown) / obj.pw_ad_steps;
			
			obj.pw_adapt = obj.cp_pw_adapt(obj.pw_adapt, pw_m, obj.pw_old{1}, (pbc | (obj.pw_ad_step_count>(obj.pw_ad_steps-5))), aup, adown);
			obj.pw_old{1} = obj.pw_old{2};

			% Choose Voltage
			FD = diag(f_ref)*hpc_class.VDom;			
			V = obj.cp_voltage_choice(hpc_class, FD);
			F = f_ref;

			% Control Temperature:
			if (mod(s,obj.Tper)==2)
				Tband = [45.0 65.0 80.0 85 ];
				Tband = Tband + 273.15;
				derBand = [0 0.5 1.0 2.0];
				derBand = derBand + hpc_class.measure_noise*hpc_class.T_noise_max/2;
				redmat =	[2	2	1	0	-1;
						 	2	1	1	0	-1;
						 	1	1	0	-1	-2;
						 	1	0	-1	-2	-3;
						 	0	0	-2	-3	-4];
	
				tr(:,2) = sum(T > Tband,2);
				tr(:,1) = sum(derT > derBand,2);
				tr = tr + 1; %indexing
	
				obj.VCred = obj.VCred + redmat((tr(:,2)-1)*size(redmat,1) + tr(:,1));
				%saturate: (no turbo atm)
				obj.VCred = obj.VCred + (obj.VCred > 0).*(0-obj.VCred);
			end
			
			VDred = min(diag(obj.VCred)*hpc_class.VDom);
			V = V + fix(VDred/2)'*0.05;

			%cap F:
			% maybe remove?
			fmaxi = hpc_class.FV_table(sum(hpc_class.VDom * V > ones(hpc_class.Nc,1)*hpc_class.FV_table(:,1)', 2) + 1, 3);
			F = F + (F>fmaxi).*(fmaxi - F);
			F = F + hpc_class.VDom*(mod(VDred,2).*sign(VDred))'*0.05;

			% Compute Power
			pu = (Ceff.*F.*(hpc_class.VDom*V) + hpc_class.leak_vdd_k).*(hpc_class.VDom*V) + d_p*hpc_class.leak_process_k;
			
			% DispatchPower
			%TODO, quad power budget dispatching
			delta_p = sum(pu) - p_budget + obj.pw_adapt;
			if (delta_p > 0)
				dompw = hpc_class.VDom'*pu;
				domT = ((hpc_class.VDom'*T ./ sum(hpc_class.VDom)') + max(T.*hpc_class.VDom)' ) / 2; %here mean*2 / 3??
				[resdpw, ~] = obj.cp_pw_dispatcher(domT, hpc_class.core_crit_temp, pid_target(1:obj.vd), delta_p, dompw, hpc_class.min_pw_red); %todo: pid_target(1:obj.vd)
				deltapd = dompw-resdpw;
				for vdi=1:hpc_class.vd
					pd = pu.*hpc_class.VDom(:,vdi);
					pd = pd(pd>0);
					cd = Ceff.*hpc_class.VDom(:,vdi);
					cd = cd(cd>0);
					cidx = [1:1:hpc_class.Nc]' .* hpc_class.VDom(:,vdi);
					cidx = cidx(cidx>0);
					[pd, ~] = obj.cp_pw_dispatcher_c(cd,deltapd(vdi), pd, hpc_class.min_pw_red);
					pu(cidx) = pd;
				end
				%pw_storage = pw_storage + pws;	
			end

			% Compute Freq
			%% Newton-Rapson		
			for vi=1:hpc_class.vd
				fp = f_ref.*hpc_class.VDom(:,vi);	
				cidx = [1:1:hpc_class.Nc]' .* (fp>0);
				cidx = cidx(cidx>0);
				lim_sup = fp(fp>0);		
				lim_inf = hpc_class.F_min*ones(sum(hpc_class.VDom(:,vi)),1);
				pup = pu.*hpc_class.VDom(:,vi);
				pup = pup(pup>0);
				Cip = Ceff.*hpc_class.VDom(:,vi);
				Cip = Cip(Cip>0);
				d_pk = (d_p*hpc_class.leak_process_k).*hpc_class.VDom(:,vi);
				d_pk = d_pk(d_pk>0);
				xi = pw_function(lim_inf, pup, Cip, hpc_class.leak_vdd_k,d_pk);			
				xs = pw_function(lim_sup, pup, Cip, hpc_class.leak_vdd_,d_pk);
				
				fo = zeros(length(cidx),1);
				
				%boundaries:
				sn = xs>0;
				so = xi>0;
				c_lim_inf = lim_inf + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);
				c_lim_sup = lim_sup + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_sup);
				%fo = fo + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);
				
				for nri=1:16
					fo = (c_lim_inf+c_lim_sup)/2;
					xs = pw_function(fo, pup, Cip, hpc_class.leak_vdd_k,d_pk);
					if (sum(abs(xs)) <= 1e-3*length(cidx)) || (sum(abs(c_lim_sup-fo))<= hpc_class.F_discretization_step*length(cidx))
						break;
					end
					sn = xs>0;
					so = xi>0;
					
					xi = xi + (sn==so).*(xs-xi);
					c_lim_inf = c_lim_inf + (sn==so).*(fo-c_lim_inf);
					c_lim_sup = c_lim_sup + (sn~=so).*(fo-c_lim_sup);				
				end
				
				%discretize
				if hpc_class.F_discretization_step > 0
					fo = round(fo/hpc_class.F_discretization_step) * hpc_class.F_discretization_step;
				end
				F(cidx) = fo;
				V(vi) = hpc_class.FV_table(sum(max(fo) > hpc_classobj.FV_table(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue
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
			v_hys_red = (pw_m - p_budget) < -(p_budget*(0.125+obj.hys_p_count/10));
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
			fmaxi = obj.FV_table(sum(obj.VDom * V > ones(obj.Nc,1)*obj.FV_table(:,1)'+1e-6, 2) + 1, 3);
			F = F + (F>fmaxi).*(fmaxi - F);
			F = F + (F<obj.F_min).*(obj.F_min*ones(obj.Nc,1) - F);
			
			%TEST: TODO REMOVE
			pu = (Ceff.*F.*(obj.VDom*V) + obj.leak_vdd/1000).*(obj.VDom*V) + d_p*obj.leak_process/1000;
			obj.pw_old{2} = sum(pu);


		end
		[] = cleanup_fnc(obj, hpc_class)
		[] = plot_fnc(obj, hpc_class)

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

