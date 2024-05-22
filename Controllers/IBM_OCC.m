classdef IBM_OCC < controller
	%CP Summary of this class goes here
	%   Detailed explanation goes here
	
	properties

		occ_dom_active = 0;
		occ_thermal_all = 1; %0
		occ_thermal_tick = 2; %16 % 2
		occ_p_gain = 1000;

		occ_speed_step = 10;	
		gf_max = 3.45; %obj.F_Max; or 3.45?

		dps_alpha_up = 0.999; 
		dps_alpha_down = 0.999;
		dps_step_up = 1000;
		dps_step_down = 8;
		% DPS favoring energy savings:
		%dps_alpha_up = 0.98; 
		%dps_alpha_down = 0.98;
		%dps_step_up = 8;
		%dps_step_down = 8;
		
	end

	properties(SetAccess=protected, GetAccess=public)	
		lf_max;
		cf_max;
		gf_min;

		core_dps;
		core_util;

		occ_thermal_tick_idx;
		cpu_speed;
		th_freq;
		st_proc_pcap_vote;
		proc_ghz_per_watt
		%Number of watts power must be below the node power cap before raising ppb_fmax
		pdrop_thresh;
		frequency_step_pstate;

		pmin;
		ptl;
		
		psm;
		node_ghz_per_watt;
		ghz_per_pstate;

	end
	
	methods
		function obj = IBM_OCC()
			%CP Construct an instance of this class
			%   Detailed explanation goes here
			
		end
		function [obj] = init_fnc(obj, hc, Nsim)
			%TODO understand which one I actually really need
			%TODO understand which needs to go out and which in

			obj.ex_count = 0;
			obj.T_target = ones(hc.Nc, 1)*hc.core_limit_temp;

			obj.lf_max = obj.gf_max;
			obj.cf_max = obj.gf_max;
			obj.gf_min = hc.F_min;

			obj.core_dps = 1000*ones(hc.Nc,1); %here they initialize at 2000, but they saturate at 1000 right away. dunno why....
			obj.core_util = ones(hc.Nc,1);

			obj.occ_thermal_tick_idx = obj.occ_thermal_tick;
			obj.cpu_speed = 1000*ones(hc.Nc,1);
			obj.th_freq = obj.gf_max*ones(hc.Nc,1);
			obj.st_proc_pcap_vote = obj.gf_max;
			obj.proc_ghz_per_watt = 1.118 / hc.Nc; % but with leakage this is really not good
			obj.pdrop_thresh = 10/15*hc.Nc; % normalized with the number of core!
			obj.frequency_step_pstate = hc.F_discretization_step; %Comes from refclk/dpll_divider attributes.

			%%%%%%%%%%%%%%%%%%5
			%function [ptl, pmin] = create_pstate_table(obj, fmax, fmin, frequency_step_pstate)	
				% assumed pmax == 0
				% assumed pmin is computed and not given
				obj.pmin = length(obj.gf_max:-obj.frequency_step_pstate:obj.gf_min);
				obj.ptl = zeros(obj.pmin+1, 3);
				obj.ptl(:,1) = 0:1:obj.pmin;
				obj.ptl(:,2) = [obj.gf_max:-obj.frequency_step_pstate:obj.gf_min obj.gf_min];
				
				tt = sum(obj.ptl(:,2).*ones(obj.pmin+1,hc.FV_levels) > hc.FV_table(:,3)',2) + 1;
				obj.ptl(:,3) = hc.FV_table(tt,1);
			%end
			%%%%%%%%%%%%%%%%%%5
				obj.psm = obj.pmin; %psm = 0 + (gf_max - gf_min) / frequency_step_pstate; %minimum pstate; 
				obj.node_ghz_per_watt = obj.proc_ghz_per_watt; %only one proc
				
				obj.ghz_per_pstate = hc.F_discretization_step; %frequency_step_khz/1000, in the main %todo


		end
		function [F,V, obj] = ctrl_fnc(obj, hc, target_index, pvt, i_pwm, i_wl)

			obj.ex_count = obj.ex_count + 1;

			f_ref = hc.frtrc(min(target_index, size(hc.frtrc,1)),:)';
			p_budget = hc.tot_pw_budget(min(target_index, length(hc.tot_pw_budget)));

			T = pvt{hc.PVT_T};
			process = pvt{hc.PVT_P};

			% ignore gpu power control
			% ignore memory power control
			pw_cpu = i_pwm; %since there is no other things that is consuming
	
			%%%%%%%%%%%%%%%
			% amec_pcap_calc: Calculate the pcap for the proc and the power capping limit for nominal cores
	
			n_avail_pw = p_budget - i_pwm;
			proc_fraction = pw_cpu / i_pwm; %we can ignore: *2^16
	
			p_budget_cpu = i_pwm + (proc_fraction* n_avail_pw); %we can ignore: /2^16
			% so basically it is p_budget_cpu = p_budget, since there is no other
			% things that is consuming and pw_cpu = i_pwm
	
			% if p_budget is not set, no pcapping will happen... ok?
			%if p_budget < obj.P_Max * obj.Nc % this is wrong! : n_avail_pw < 0
				pcap_fmin = obj.gf_min;
			%else
			%	pcap_fmin = F; %which F ?? max()?
			%end
	
			%%%%%%%%%%%%%%%
			% amec_pcap_controller: Calculate voting box input freq for staying with the current pcap
	
			avail_pw = p_budget_cpu - pw_cpu;
	
			if (obj.st_proc_pcap_vote > pcap_fmin)	
				% always this
				obj.st_proc_pcap_vote = obj.cf_max + (obj.proc_ghz_per_watt * avail_pw);
			else
				obj.st_proc_pcap_vote = obj.st_proc_pcap_vote + (obj.proc_ghz_per_watt * avail_pw);
			end
			%saturate:
			obj.st_proc_pcap_vote = obj.st_proc_pcap_vote + (obj.st_proc_pcap_vote>obj.gf_max)*(obj.gf_max-obj.st_proc_pcap_vote) + ...
								(obj.st_proc_pcap_vote<obj.gf_min)*(obj.gf_min-obj.st_proc_pcap_vote);


			%%%%%%%%%%%%%%%
			% ppb_fmax_calc: Calculate the performance preserving bounds voting box input freq
			%n_avail_pw = p_budget - pw_m;
			if n_avail_pw <= 0
				drop = obj.node_ghz_per_watt*(-n_avail_pw);
				% Need to shed power, make sure we drop by at least 1 Pstate
				drop = drop + (drop<obj.ghz_per_pstate).*(obj.ghz_per_pstate-drop);
				% prevent overflow
				if max(drop) < obj.lf_max
					obj.lf_max = obj.lf_max-max(drop);
				else
					obj.lf_max = obj.gf_min;
				end
			elseif (n_avail_pw > obj.pdrop_thresh)
				obj.lf_max = obj.lf_max + obj.ghz_per_pstate;
			end
			%saturate
			obj.lf_max = obj.lf_max + (obj.lf_max>obj.gf_max)*(obj.gf_max-obj.lf_max) + ...
						(obj.lf_max<obj.gf_min)*(obj.gf_min-obj.lf_max);
					
			%%%%%%%%%%%%%%%
			% amec_controller_proc_thermal: thermal
			obj.occ_thermal_tick_idx = obj.occ_thermal_tick_idx - 1;
			if obj.occ_thermal_tick_idx <= 0
				obj.occ_thermal_tick_idx = obj.occ_thermal_tick;
	
				% Get hottest core temperature in OCC processor
				if obj.occ_dom_active
					occ_temp = obj.VDom*max(T.*hc.VDom)';
				else
					occ_temp = max(T)*ones(hc.Nc,1);
				end
				if obj.occ_thermal_all
					occ_temp = T;
				end
	
				% Occ thermal setpoint is 85°C
				%	it gets lowered if there are some errors.
				error = obj.T_target - occ_temp;
	
				% Proportional Controller for the thermal control loop
				residue = error * obj.occ_p_gain;
	
				% stuff
				olimit = floor(65535/4/obj.occ_speed_step);
				tc = floor((residue*10)/(2^16));
				%saturation
				tc = tc + (tc>olimit).*(olimit-tc) + (tc< -olimit).*(-olimit-tc);
	
				%Calculate the new thermal CPU speed request
				obj.cpu_speed = obj.cpu_speed + tc*obj.occ_speed_step;
	
				%Proceed with residue summation to correctly follow set-point
				if residue < 0
					obj.cpu_speed = obj.cpu_speed + obj.occ_speed_step;
				end
	
				%saturate
				%cpu_speed = cpu_speed + (cpu_speed>1000).*(1000-cpu_speed) + (cpu_speed<400).*(400-cpu_speed);
				obj.cpu_speed = obj.cpu_speed + (obj.cpu_speed>1000).*(1000-obj.cpu_speed) + (obj.cpu_speed<200).*(200-obj.cpu_speed);
	
				% speed2freq
				% to handle max freq changing (i.e. mode change)
				obj.th_freq = (obj.cpu_speed >= 1000).*obj.gf_max + (obj.cpu_speed < 1000).*obj.gf_max .* obj.cpu_speed ./ 1000;
	
			end %thermal

			%%%%%%%%%%%%%%%
			% amec_mst_gen_soft_freq: generate soft frequency boundaries, softmin and softmax
			% these 2 values are supposed to be "for each part", but looking at the
			% code, there is ONLY ONE PART which is 0.
			% I will treat them in the same way as the processor vs domain.
	
			%Use the frequency delta sent by the customer
			%fdelta = ??
			%https://github.com/open-power/occ/blob/16131c38c2e593206447403c3ac79b18fb86e21a/src/occ_405/amec/amec_master_smh.c#L726
	
			% ftarget = ??
			% it is based on dps and on modes
	
			%soft_fmin = ftarget - fdelta;
			%soft_fmax = ftarget + fdelta;
	
			%%%%%%%%%%%%%%%
			% DPS
	
			% for each core
			% Update moving average of util_slack and util_active for all cores
			% They use a "moving average" with a memory of 64 units
			% i will use real MA:
			occ_alpha_util = 0.1; % 0.9^64 = 0.0012
			cutl = sum(i_wl(:,2:end),2);
			obj.core_util = obj.core_util*(1-occ_alpha_util) + cutl*occ_alpha_util;
			obj.core_util = obj.core_util + (obj.core_util>1).*(1-obj.core_util);
	
			obj.core_dps = obj.core_dps + (obj.core_util > obj.dps_alpha_up)*obj.dps_step_up - ...
							(obj.core_util < obj.dps_alpha_down)*obj.dps_step_down;
						
			%saturate
			obj.core_dps = obj.core_dps + (obj.core_dps > 1000).*(1000-obj.core_dps) + ...
				(obj.core_dps < 400).*(400-obj.core_dps);
	
			dps_freq = (obj.core_dps >= 1000).*obj.gf_max + ...
				(obj.core_dps < 1000).*obj.gf_max .* obj.core_dps ./ 1000;
			%we are not using dps, cuz we are measuring performance and not
			%energy without any wait time for idling
			dps_freq = obj.gf_max*ones(hc.Nc,1);
			
		
			%%%%%%%%%%%%%%%
			% Voting box:
			l_freq = min(obj.gf_max, obj.lf_max);
			% not found:
			% clip controller: pmax_clip_freq
	
			% thermal
			core_freq = min(obj.th_freq, l_freq);
			
			% for every core in cpu (or domain occ_dom_active)
			%check dps freq request
			core_freq = core_freq + (dps_freq<core_freq).*(dps_freq-core_freq);
			%Adjust frequency based on soft frequency boundaries
			%core_freq = core_freq + (soft_fmin(d)>core_freq)*(soft_fmin(d)-core_freq) + ...
			%			 (soft_fmax(d)<core_freq)*(soft_fmax(d)-core_freq);
	
			% why st_proc_pcap_vote here????
			% IPS FREQUENCY slv_ips_freq_request / amec_mst_ips_main
			% no need to do IPS (Idle Power Saving) because it will shut down
			% the whole CPU / domain, and we never want that
			%ips_freq = gf_max;		
			core_freq = min(core_freq,obj.st_proc_pcap_vote);
			
			core_freq = min(core_freq, f_ref);
	
			% AutoSlew: auto-slewing of frequency based on manufacturing
			% parameters: amec_master_auto_slew
			% so I will not do it
			%if foverride
			%	core_freq = foverride;
			%end
			%if pstate_foverride
			%	core_freq = pstate_foverride;
			%end
	
			%sat on f_min
			core_freq = core_freq + (core_freq<obj.gf_min).*(obj.gf_min-core_freq);
			%end
	
			obj.cf_max = max(core_freq);
			
			%%%%%%%%%%%%%%%
			% amec_slv_freq_smh: Frequency state machine:
			%{
			pmax_chip = 0;
			for vidx=1:obj.vd
				pmax = 0; % max
				fc = core_freq.*obj.VDom(:,vidx);
				fc = fc(fc>0);
				
				pstate = 0 + floor((gf_max - fc) / frequency_step_pstate);
				pstate = pstate + (pstate<0).*(0-pstate) + (pstate<psm).*(psm-pstate);
				pmax = max(pstate);
				
				%F = F + fc;			
			end
			%}
			% I'm redoing it in a way that makes sense for us
			F = core_freq;
			%V = zeros(obj.Nc,1);
			%if obj.ctrl_fixedv
			%	V = obj.ctrl_fixedv_value*ones(obj.vd, 1);
			%else
				for vidx=1:hc.vd
					fc = F.*hc.VDom(:,vidx);
					fcm = max(fc);
					V(vidx,:) = hc.FV_table(sum(fcm > hc.FV_table(:,3))+1,1);
					%V = V + vc*obj.VDom(:,vidx);			
				end
			%end
			
			%ignore manufacturing comands
	
			%ignore memory voting box		
					
			%
			%
		end
		function [obj] = cleanup_fnc(obj, hc)
		end
		function [obj] = plot_fnc(obj, hc, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop)
		end
		
	end
		
end

