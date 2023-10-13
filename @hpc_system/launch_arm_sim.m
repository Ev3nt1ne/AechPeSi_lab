function [cpxplot, cpuplot, cpfplot, cpvplot, wlop] = launch_arm_sim(obj, robust, mode, show)
	if (nargin < 2) || isempty(robust)
		robust = 0;
	end
	if (nargin < 3) || isempty(mode)
		mode = 0;
		% 0: P->T (ControlPulp)
		% 1: T->P (Inverse)
		% 2: T&P (Voting)
	end
	if (nargin < 4) || isempty(show)
		show = 1;
	end
	
	
	%%%%%%%%%%% 
	% ARM INIT
	
	arm_one_th = 1;
	
	switch_on_temperature = 50 + 273.15;
	tdp = 10;
	pid_ki = 1;
	integral_cutoff = 0;
	integral_max = 100;
	kp_overshoot = 1;
	kp_undeshoot = 1;	
	integral_error = zeros(obj.Nc,1);
	thermal_allocatable_power = obj.VDom *(tdp*sum(obj.VDom))';
	
	granted_power = tdp*ones(obj.Nc,1);
	carry_over_power = zeros(obj.vd,1);
	
	gf_max = 3.45; %obj.F_Max; or 3.45?
	perf_lim = obj.FV_table(sum(gf_max > obj.FV_table(:,3)+1e-6)+1,1) * ones(obj.vd,1);
	
	
	%%%%%%%%%%%%%%

	d_is = [ones(obj.Nc,1) zeros(obj.Nc,obj.ipl-1)];
	wl = d_is;
	d_p = ones(obj.Nc,1);
	%obj.usum = ((obj.core_Max_power-0.5)*obj.Ni_c-obj.Ni_c)*ones(1+obj.vd,1);

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
	F_MA = zeros(obj.Nc,1);
	pid_target = ones(obj.Nc, 1)*obj.core_crit_temp;
	if robust
		pid_target = pid_target - ones(obj.Nc, 1)*obj.T_margin;
	end
	pid_integ = zeros(obj.Nc,1);

	cpuplot = zeros(Nsim*sim_mul+1,obj.Ni_c);
	cpxplot = zeros(Nsim*sim_mul+1,obj.Ns);
	cpxplot(1,:) = x;
	cpuplot(1,:) = NaN;
	cpfplot = zeros(Nsim+1,obj.Nc);
	cpvplot = zeros(Nsim+1,obj.vd);
	cpfplot(1,:) = F;
	cpvplot(1,:) = V;
	F_og = F;

	pw_adapt = 0;
	pu = obj.min_pw_red;
	pw_ms = sum(pu);
	pbc = 0;
	
	fref_changed = zeros(obj.Nc,1);
	
	obj = obj.init_compute_model(Adl_true, Bdl_true);

	% LOOOP
	for s=1:Nsim

		%Read T, d_i, p_budget, f_ref
		if (mod(s, div_comms) == 0)
			ci_index = ci_index + 1;
			fref_old = f_ref;
			f_ref = obj.frplot(min(ci_index, size(obj.frplot,1)),:)';
			fref_changed = f_ref~=fref_old;				
			pold = p_budget;
			p_budget = obj.tot_pw_budget(min(ci_index, length(obj.tot_pw_budget)));
			if p_budget~=pold
				pbc = 1;
			else
				pbc = 0;
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
		
		% pre-do stuff for below:
		% Process Workload
		[wl, Ceff] = obj.cp_wl_process(d_i, wl);
		%https://github.com/ARM-software/SCP-firmware/blob/master/module/scmi_perf/src/perf_plugins_handler.c#L331
		%they take the max perf, the min max_limit, and the max min_limit
		
		%https://github.com/ARM-software/SCP-firmware/blob/master/module/traffic_cop/src/mod_traffic_cop.c#L326
		% there is a limiter depending on the number of active cores. I
		% don't care about this because all cores are active
		
		%F = f_ref;
		%reset on change:
		tt = sum(fref_changed.*obj.VDom)';
		perf_lim = perf_lim + ...
			(tt|pbc) .* ( obj.FV_table(sum(gf_max > obj.FV_table(:,3)+1e-6)+1,1) - perf_lim);
		
		fmaxi = obj.FV_table(sum(obj.VDom * perf_lim > ones(obj.Nc,1)*obj.FV_table(:,1)'+1e-6, 2) + 1, 3);
		F = min(f_ref,fmaxi);
		
		FD = diag(F)*obj.VDom; %diag(f_ref)*obj.VDom;			
		V = obj.cp_voltage_choice(FD);

		demand_power = (Ceff.*F.*(obj.VDom*V) + obj.leak_vdd/1000).*(obj.VDom*V) + d_p*obj.leak_process/1000;

		%%%%%%%%%%%%%%%
		%control update
		
		%update when data are available (?)
		
		%%%%%%%%%%%%
		%pi_control (for each dev_id = THERMAL DOMAIN)
		for vdi=1:obj.vd

			%{
			details The temperature that the system will achive once stabilised.
				  Due to the PI nature of the controller, some
				  overshoot/undershoot may occur. Note that the controller can
				  only limit the temperature by placing a limit to the power to
				  the heat source. It has no direct control on the heat source
				  itself and therefore only the upper limit can be controlled.
			%}
			% their pid target control temperature.
			%control_temperature = 60;
			%could not find but I assume is max(T)
			tdpd = tdp*sum(obj.VDom(:,vdi));
			if arm_one_th
				td = max(T.*obj.VDom(:,vdi))';
			else
				td = (T.*obj.VDom(:,vdi));
				td = td(td>0);
			end
			cidx = [1:1:obj.Nc]' .* obj.VDom(:,vdi);
			cidx = cidx(cidx>0);
			
			err = pid_target(1:length(td)) - td;
			pid_kp = (err<0).*kp_overshoot + (err>=0).*kp_undeshoot;

			%https://github.com/ARM-software/SCP-firmware/blob/master/module/thermal_mgmt/src/mod_thermal_mgmt.c#L37
			%if (T<switch_on_temperature)
				%integral_error = 0
				%thermal_allocatable_power = tdp;
			%end
			% pid will work only if T>=switch_on_temperature
			ttpd = td>=switch_on_temperature;
			integral_error(cidx) = integral_error(cidx).*ttpd; %zeros it
			thermal_allocatable_power(cidx) = thermal_allocatable_power(cidx) + (~ttpd).*(tdpd-thermal_allocatable_power(cidx));

			terr = integral_error(cidx) + err;

			tt = (err<integral_cutoff) && (terr(1:length(err)) < integral_max);
			tt = tt & ttpd;
			integral_error(cidx) = integral_error(cidx) + tt.*(terr-integral_error(cidx));

			pi_power = (pid_kp * err) + (pid_ki * integral_error(cidx));

			tt = (pi_power + tdpd) > 0;
			%tt = tt & ttpd;
			thermal_allocatable_power(cidx) = thermal_allocatable_power(cidx) + ...
				(tt&ttpd).*(pi_power + tdpd - thermal_allocatable_power(cidx)) + ...
				(ttpd&~tt).*(0- thermal_allocatable_power(cidx));

			tot_spare_power = 0;

			%%%%%%%%%%%%%%%%%%%%
			%thermal protection
			%not doing it

			%%%%%%%%%%%%%%%%%%%
			%distribute power
			%idle_power = 0;
			%tot_weighted_demand_power = 0;
			%tot_spare_power = 0;
			%tot_power_deficit = 0;

			%STEP 0:
			%Initialise the actors' demand power.
			%get power: done above:

			%for each core:
				%activity is idle vs full, I will just put 1
				activity = 1024;
				prev_used_power = (granted_power(cidx) .* activity) / 1024;
				idle_power = sum(granted_power(cidx) - prev_used_power);

				%so:
				idle_power = 0;

				% weight is the same for all
				weight = 1;
				tot_weighted_demand_power = sum(weight .* demand_power(cidx));

			allocatable_power = thermal_allocatable_power + idle_power;


			%STEP 1:
			%The power available is allocated in proportion to the actors' weight and
			%	their power demand.

			%for each core:
				%allocate_power
				granted_power(cidx) = (weight .* demand_power(cidx)) .* (allocatable_power(cidx) / tot_weighted_demand_power);

				tt = granted_power(cidx) > demand_power(cidx);			
				spare_power = tt.*(granted_power(cidx) - demand_power(cidx)); %otherwise 0;
				power_deficit = ~tt.*(demand_power(cidx) - granted_power(cidx)); %otherwise 0
				granted_power(cidx) = granted_power(cidx) + tt.*(demand_power(cidx)-granted_power(cidx));

				tot_spare_power = sum(spare_power);
				tot_power_deficit = sum(power_deficit);

			%STEP 2:
			%A further allocation based on actors' power deficit, if any spare power
			%	left.
			%	Finally, get the corresponding performance level and place a new limit.
			tot_spare_power = tot_spare_power + carry_over_power(vdi);
			carry_over_power(vdi) = 0;
			%for each core:
				%re_allocate_power
				if (tot_spare_power > 0)
					tt = (power_deficit > 0);
					%The actor has been given less than requested, and it may still take
					%	some power
					granted_power(cidx) = granted_power(cidx) + ...
							tt.*(power_deficit .* tot_spare_power ./ (~tt|tot_power_deficit));

					tr = granted_power(cidx) > demand_power(cidx);
					carry_over_power(vdi) = carry_over_power(vdi) + ...
									sum( (tr.*tt).*(granted_power(cidx) - demand_power(cidx)) );
					granted_power(cidx) = granted_power(cidx) + (tr.*tt).*(demand_power(cidx)-granted_power(cidx));
					%else
					%The actor has received the power it requested. The amount of
					%	power left can be used in the next fast-loop.
					carry_over_power(vdi) = carry_over_power(vdi) + sum( (~tr.*tt).*spare_power );
				end
				
				%power2freqV
				F(cidx) = (((granted_power(cidx) - d_p(cidx)*obj.leak_process/1000) ./ V(vdi)) - obj.leak_vdd/1000) ./ (V(vdi)) ./ Ceff(cidx);
				
				F = F + (F<obj.F_min).*(obj.F_min*ones(obj.Nc,1) - F);
				
				%If we have been granted the power we requested (which is at most the
				%	limit already placed by other plugins), then there's no need to limit
				%	further.
				new_perf_limit = max(F(cidx));
				V_n = obj.FV_table(sum(new_perf_limit > obj.FV_table(:,3)+1e-6)+1,1);
				tt = max(granted_power(cidx) < demand_power(cidx)) && (V_n < perf_lim(vdi));
				if tt
					V(vdi) = V_n;
					perf_lim(vdi) = V_n;
				end				
				
		end % for vdi=1:obj.vd
			

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

		obj.xutplot(t1, cpxplot,cpuplot);
		obj.powerconstrplot(t1,cpuplot);
		obj.tempconstrplot(cpxplot);
		obj.perfplot(cpfplot, obj.wl_index);

		obj.fvplot(t2, cpfplot,cpvplot);
	end
	
	wlop = obj.wl_index / (size(obj.wrplot,3)-1) * 100;
end
