classdef Fuzzy < CP
	%FUZZY Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		%global var
		pw_storage = 0;
		pw_adapt = 0;
		pu_old = 0;
		pw_ad_step_count = 0;
		wl;
		VCred;

		%parameters
		pw_ad_aup = 0.015; %0.05; %0.025;
		pw_ad_adown = 0.02; %0.075; %0.025; %0.1;
		pw_ad_achange = 0.5; %0.5;
		pw_ad_steps = 20;
	end
	
	methods
		function obj = Fuzzy()
			%FUZZY Construct an instance of this class
			%   Detailed explanation goes here
			

			% 
		end
		
		function [F,V, obj] = ctrl_fnc(obj, hpc_class, target_index, pvt, i_pwm, i_wl)

			f_ref
			pbc
			s 
			Tper
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
			
			obj.pw_adapt = obj.cp_pw_adapt(obj.pw_adapt, pw_m, obj.pu_old, (pbc | (obj.pw_ad_step_count>(obj.pw_ad_steps-5))), aup, adown);
			obj.pu_old = sum(pu);

			% Choose Voltage
			FD = diag(f_ref)*hpc_class.VDom;			
			V = obj.cp_voltage_choice(hpc_class, FD);
			F = f_ref;

			% Control Temperature:
			if (mod(s,Tper)==2)
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
					td = T.*hpc_class.VDom(:,vdi);
					td = td(td>0);
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

