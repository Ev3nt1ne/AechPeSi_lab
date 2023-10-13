function [cpxplot, cpuplot, cpfplot, cpvplot, wlop] = launch_occ_sim(obj, Tall, Ttick, robust, show)
%LAUNCH_OCC_SIM Summary of this function goes here
%   Detailed explanation goes here

	if (nargin < 2) || isempty(Tall)
			Tall = 0;
	end
	if (nargin < 3) || isempty(Ttick)
		Ttick = 16;
	end
	if (nargin < 4) || isempty(robust)
		robust = 0;
	end
	if (nargin < 5) || isempty(show)
		show = 1;
	end
	
	%% OCC INIT
	occ_dom_active = 0;
	occ_thermal_all = Tall; %0
	occ_thermal_tick = Ttick; %16 % 2
	occ_p_gain = 1000;
	%if obj.ctrl_fixedv
	%	occ_p_gain = occ_p_gain*30;
	%end
	occ_speed_step = 10;	
	gf_max = 3.45; %obj.F_Max; or 3.45?
	
	lf_max = gf_max;
	cf_max = gf_max;
	gf_min = obj.F_min;
	
	%
	core_dps = 1000*ones(obj.Nc,1); %here they initialize at 2000, but they saturate at 1000 right away. dunno why....
	core_util = ones(obj.Nc,1);
	dps_alpha_up = 0.999; 
	dps_alpha_down = 0.999;
	dps_step_up = 1000;
	dps_step_down = 8;
	% DPS favoring energy savings:
	%dps_alpha_up = 0.98; 
	%dps_alpha_down = 0.98;
	%dps_step_up = 8;
	%dps_step_down = 8;
	occ_thermal_tick_idx = occ_thermal_tick;
	cpu_speed = 1000*ones(obj.Nc,1);
	th_freq = gf_max*ones(obj.Nc,1);
	st_proc_pcap_vote = gf_max;
	proc_ghz_per_watt = 1.118 / obj.Nc; % but with leakage this is really not good
	%Number of watts power must be below the node power cap before raising ppb_fmax
	pdrop_thresh = 10/15*obj.Nc; % normalized with the number of core!
	frequency_step_pstate = obj.F_discretization_step; %Comes from refclk/dpll_divider attributes.	
%%%%%%%%%%%%%%%%%%5
%function [ptl, pmin] = create_pstate_table(obj, fmax, fmin, frequency_step_pstate)	
	% assumed pmax == 0
	% assumed pmin is computed and not given
	pmin = length(gf_max:-frequency_step_pstate:gf_min);
	ptl = zeros(pmin+1, 3);
	ptl(:,1) = 0:1:pmin;
	ptl(:,2) = [gf_max:-frequency_step_pstate:gf_min gf_min];
	
	tt = sum(ptl(:,2).*ones(pmin+1,obj.FV_levels) > obj.FV_table(:,3)',2) + 1;
	ptl(:,3) = obj.FV_table(tt,1);
%end
%%%%%%%%%%%%%%%%%%5
	psm = pmin; %psm = 0 + (gf_max - gf_min) / frequency_step_pstate; %minimum pstate; 
	node_ghz_per_watt = proc_ghz_per_watt; %only one proc
	
	ghz_per_pstate = obj.F_discretization_step; %frequency_step_khz/1000, in the main %todo
	
	%%
	d_is = [ones(obj.Nc,1) zeros(obj.Nc,obj.ipl-1)];
	wl = d_is;
	d_p = ones(obj.Nc,1);
	obj.usum = ((obj.core_Max_power-0.5)*obj.Ni_c-obj.Ni_c)*ones(1+obj.vd,1);

	T = obj.C(1:obj.Nc,:)*obj.x_init;
	x = obj.x_init + (rand(obj.Ns,1) - 0.5*ones(obj.Ns,1));
	sim_mul = ceil(obj.Ts_ctrl/obj.Ts);
	Nsim = ceil(obj.tsim / obj.Ts_ctrl);
	div_comms = ceil(obj.Ts_input / obj.Ts_ctrl);
	ci_index = 1;
	f_ref = obj.frplot(ci_index,:)';
	p_budget = obj.tot_pw_budget(ci_index);
	F = obj.F_min*ones(obj.Nc,1);
	V = obj.V_min*ones(obj.vd,1);
	Adl_true = obj.Ad_true;
	Bdl_true = obj.Bd_true;
	%
	pid_target = obj.core_crit_temp;
	if robust
		pid_target = pid_target - ones(obj.Nc, 1)*obj.T_margin;
	end

	cpuplot = zeros(Nsim*sim_mul+1,obj.Ni_c);
	cpxplot = zeros(Nsim*sim_mul+1,obj.Ns);
	cpxplot(1,:) = x;
	cpuplot(1,:) = NaN;
	cpfplot = zeros(Nsim+1,obj.Nc);
	cpvplot = zeros(Nsim+1,obj.vd);
	cpfplot(1,:) = F;
	cpvplot(1,:) = V;

	pw_ms = sum(obj.min_pw_red);
	
	obj = obj.init_compute_model(Adl_true, Bdl_true);
	
	pbp = 0;

	% LOOOP
	for s=1:Nsim

		%Read T, d_i, p_budget, f_ref
		if (mod(s, div_comms) == 0)
			ci_index = ci_index + 1;
			f_ref = obj.frplot(ci_index,:)';
			p_budget = obj.tot_pw_budget(min(ci_index, length(obj.tot_pw_budget)));
			if p_budget~=pbp
				pbp = p_budget;
			end
		end
		if (s~=1)
			%x = cpxplot(index+sim_mul,:)';
			T = obj.C(1:obj.Nc,:)*cpxplot(index+sim_mul,:)';
			%add noise:
			if obj.measure_noise
				nn = (rand(obj.Nc,1) - 0.5)*2 * obj.T_noise_max;
				T = T + nn;
			end
		end
		d_i = d_is;
		pw_m = pw_ms;

		%Compute model:
		index = 1+(s-1)*sim_mul;
		[cpuplot(index+1:index+sim_mul,:), cpxplot(index+1:index+sim_mul,:), d_is, pw_ms, obj] = obj.compute_model(sim_mul, cpxplot(index,:)', V, F, d_p);	

		%
		%
		
		% ignore gpu power control
		% ignore memory power control
		pw_cpu = pw_m; %since there is no other things that is consuming

		%%%%%%%%%%%%%%%
		% amec_pcap_calc: Calculate the pcap for the proc and the power capping limit for nominal cores

		n_avail_pw = p_budget - pw_m;
		proc_fraction = pw_cpu / pw_m; %we can ignore: *2^16

		p_budget_cpu = pw_m + (proc_fraction* n_avail_pw); %we can ignore: /2^16
		% so basically it is p_budget_cpu = p_budget, since there is no other
		% things that is consuming and pw_cpu = pw_m

		% if p_budget is not set, no pcapping will happen... ok?
		%if p_budget < obj.P_Max * obj.Nc % this is wrong! : n_avail_pw < 0
			pcap_fmin = gf_min;
		%else
		%	pcap_fmin = F; %which F ?? max()?
		%end

		%%%%%%%%%%%%%%%
		% amec_pcap_controller: Calculate voting box input freq for staying with the current pcap

		avail_pw = p_budget_cpu - pw_cpu;

		if (st_proc_pcap_vote > pcap_fmin)	
			% always this
			st_proc_pcap_vote = cf_max + (proc_ghz_per_watt * avail_pw);
		else
			st_proc_pcap_vote = st_proc_pcap_vote + (proc_ghz_per_watt * avail_pw);
		end
		%saturate:
		st_proc_pcap_vote = st_proc_pcap_vote + (st_proc_pcap_vote>gf_max)*(gf_max-st_proc_pcap_vote) + ...
							(st_proc_pcap_vote<gf_min)*(gf_min-st_proc_pcap_vote);

		%%%%%%%%%%%%%%%
		% ppb_fmax_calc: Calculate the performance preserving bounds voting box input freq
		%n_avail_pw = p_budget - pw_m;
		if n_avail_pw <= 0
			drop = node_ghz_per_watt*(-n_avail_pw);
			% Need to shed power, make sure we drop by at least 1 Pstate
			drop = drop + (drop<ghz_per_pstate).*(ghz_per_pstate-drop);
			% prevent overflow
			if max(drop) < lf_max
				lf_max = lf_max-max(drop);
			else
				lf_max = gf_min;
			end
		elseif (n_avail_pw > pdrop_thresh)
			lf_max = lf_max + ghz_per_pstate;
		end
		%saturate
		lf_max = lf_max + (lf_max>gf_max)*(gf_max-lf_max) + ...
					(lf_max<gf_min)*(gf_min-lf_max);
				
		%%%%%%%%%%%%%%%
		% amec_controller_proc_thermal: thermal
		occ_thermal_tick_idx = occ_thermal_tick_idx - 1;
		if occ_thermal_tick_idx <= 0
			occ_thermal_tick_idx = occ_thermal_tick;

			% Get hottest core temperature in OCC processor
			if occ_dom_active
				occ_temp = obj.VDom*max(T.*obj.VDom)';
			else
				occ_temp = max(T)*ones(obj.Nc,1);
			end
			if occ_thermal_all
				occ_temp = T;
			end

			% Occ thermal setpoint is 85°C
			%	it gets lowered if there are some errors.
			error = pid_target - occ_temp;

			% Proportional Controller for the thermal control loop
			residue = error * occ_p_gain;

			% stuff
			olimit = floor(65535/4/occ_speed_step);
			tc = floor((residue*10)/(2^16));
			%saturation
			tc = tc + (tc>olimit).*(olimit-tc) + (tc< -olimit).*(-olimit-tc);

			%Calculate the new thermal CPU speed request
			cpu_speed = cpu_speed + tc*occ_speed_step;

			%Proceed with residue summation to correctly follow set-point
			if residue < 0
				cpu_speed = cpu_speed + occ_speed_step;
			end

			%saturate
			%cpu_speed = cpu_speed + (cpu_speed>1000).*(1000-cpu_speed) + (cpu_speed<400).*(400-cpu_speed);
			cpu_speed = cpu_speed + (cpu_speed>1000).*(1000-cpu_speed) + (cpu_speed<200).*(200-cpu_speed);

			% speed2freq
			% to handle max freq changing (i.e. mode change)
			th_freq = (cpu_speed >= 1000).*gf_max + (cpu_speed < 1000).*gf_max .* cpu_speed ./ 1000;

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
		cutl = sum(d_i(:,2:end),2);
		core_util = core_util*(1-occ_alpha_util) + cutl*occ_alpha_util;
		core_util = core_util + (core_util>1).*(1-core_util);

		core_dps = core_dps + (core_util > dps_alpha_up)*dps_step_up - ...
						(core_util < dps_alpha_down)*dps_step_down;
					
		%saturate
		core_dps = core_dps + (core_dps > 1000).*(1000-core_dps) + (core_dps < 400).*(400-core_dps);

		dps_freq = (core_dps >= 1000).*gf_max + (core_dps < 1000).*gf_max .* core_dps ./ 1000;
		%we are not using dps, cuz we are measuring performance and not
		%energy without any wait time for idling
		dps_freq = gf_max*ones(obj.Nc,1);
		
	
		%%%%%%%%%%%%%%%
		% Voting box:
		l_freq = min(gf_max, lf_max);
		% not found:
		% clip controller: pmax_clip_freq

		% thermal
		core_freq = min(th_freq, l_freq);
		
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
		core_freq = min(core_freq,st_proc_pcap_vote);
		
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
		core_freq = core_freq + (core_freq<gf_min).*(gf_min-core_freq);
		%end

		cf_max = max(core_freq);
		
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
		if obj.ctrl_fixedv
			V = obj.ctrl_fixedv_value*ones(obj.vd, 1);
		else
			for vidx=1:obj.vd
				fc = F.*obj.VDom(:,vidx);
				fcm = max(fc);
				V(vidx) = obj.FV_table(sum(fcm > obj.FV_table(:,3))+1,1);
				%V = V + vc*obj.VDom(:,vidx);			
			end
		end
		
		%ignore manufacturing comands

		%ignore memory voting box		
				
		%
		%
		cpfplot(s+1,:) = F;
		cpvplot(s+1,:) = V;

	end %for loop

	% PLOTs
	if show
		t1 = obj.Ts*[0:Nsim*sim_mul]';
		%dimf = t1(end)/max(obj.Ts_input,obj.Ts_ctrl);
		%t2 = obj.Ts_input*[0:dimf]';
		t2 = obj.Ts_ctrl*[0:Nsim]';

		obj.xutplot(cpxplot,cpuplot);
		obj.powerconstrplot(cpuplot);
		obj.tempconstrplot(cpxplot);
		obj.perfplot(cpfplot, obj.wl_index);

		obj.fvplot(t2, cpfplot,cpvplot);
	end
	
	wlop = obj.wl_index / (size(obj.wrplot,3)-1) * 100;

end

%%


