function [cpxplot,cpuplot, cpfplot, cpvplot, xlplot, wlop] = launch_nmpc_sim(obj, obs)
%LAUNCH_NMPC_SIM Summary of this function goes here
%   Detailed explanation goes here

	if (nargin < 2) || isempty(obs)
		obs = 1;
	end

	d_is = [ones(obj.Nc,1) zeros(obj.Nc,obj.ipl-1)];
	wl = d_is;
	d_p = ones(obj.Nc,1);
	
	% TODO:
	other_states = 4;
	air_pos = 0;
	mb_pos = 1;
	pcb_pos = 2;
	al_pos = 3;

	mul_comms = floor(obj.Ts_mpc/obj.Ts_input);
	div_comms = ceil(obj.Ts_input / obj.Ts_ctrl);
	T = obj.C(1:obj.Nc,:)*obj.x_init;
	Tobs = obj.x_init;
	x = obj.x_init + (rand(obj.Ns,1) - 0.5*ones(obj.Ns,1));
	sim_mul = ceil(obj.Ts_mpc/obj.Ts);
	Nsim = ceil(obj.tsim / obj.Ts_mpc);
	f_ref = obj.frplot(1,:)';
	ci_index = 1;
	F = obj.F_min*ones(obj.Nc,1);
	V = obj.V_min*ones(obj.vd,1);
	Adl_true = obj.Ad_true;
	Bdl_true = obj.Bd_true;
	%
	cpuplot = zeros(Nsim*sim_mul+1,obj.Ni_c);
	cpxplot = zeros(Nsim*sim_mul+1,obj.Ns);
	cpxplot(1,:) = x;
	cpuplot(1,:) = NaN;
	cpfplot = zeros(Nsim+1,obj.Nc);
	cpvplot = zeros(Nsim+1,obj.vd);
	cpfplot(1,:) = F;
	cpvplot(1,:) = V;
	%V_s = V;
	%F_s = F;

	pw_adapt = 0;
	pu = obj.min_pw_red;
	pw_ms = sum(pu);
	pu_prev = pu;

	obj = obj.lin_mpc_setup2();						
	failed = 0;
	
	% Observer
	Adl_obs = obj.Ad_obs;
	Bdl_obs = obj.Bd_obs;
	poles = ones(obj.Ns,1);
	poles(1:2:obj.Ns-other_states) = obj.Obs_poles(1);
	poles(2:2:obj.Ns-other_states) = obj.Obs_poles(2);
	poles(end-other_states+1:end) = 0.2;
	%
	if obs
		LK = place(Adl_obs', obj.C(1:obj.Nc,:)', poles)';
	else
		LK = 0;
	end
	%obs_mul = round(obj.Ts_mpc/obj.Ts_obs);
	%obs_div = sim_mul / obs_mul;
	C_obs = eye(obj.Ns);
	C_obs(2:2:end-other_states,:) = 0;
	C_obs(end-al_pos,:) = 0;
	C_obs(end-pcb_pos,:) = 0;
	%
	xlplot = zeros(Nsim+1, obj.Ns);
	xl = obj.x_init;
	xlplot(1,:) = xl;
	
	index = 1;
	
	p_budget = obj.tot_pw_budget(ci_index);
	pold = p_budget;
	pbc = 0;
	pw_ad_steps = 20;
	pw_ad_step_count = 0;
	pw_ad_aup = 0.015; %0.05; %0.025;
	pw_ad_adown = 0.02; %0.075; %0.025; %0.1;
	pw_ad_achange = 0.5; %0.5;
	pu_old = pw_ms;
	
	% wl stuff
	%quantum_storage = obj.quantum_instr*ones(obj.Nc, 1);
	%wl_index = ones(obj.Nc, 1);
	
	obj = obj.init_compute_model(Adl_true, Bdl_true);
	
	pwap = zeros(Nsim, 1);
	
	[h0, ~, ~] = obj.createPsApprox();
	T0v = ones(obj.Nc,length(obj.mpc_h0_T))*diag(obj.mpc_h0_T );
	F0v = ones(obj.Nc,length(obj.mpc_h0_F))*diag(obj.mpc_h0_F );
	
	% LOOOP
	for s=1:Nsim

		%Read T, d_i, p_budget, f_ref
		if mul_comms > 0
			f_ref = obj.frplot(min(s*mul_comms,size(obj.frplot,1)),:)';
			
			p_budget = obj.tot_pw_budget(min(s*mul_comms, length(obj.tot_pw_budget)));
			
			for ni=1:obj.Nhzn
				obj.usum(ni,:) = [obj.tot_pw_budget(min((s+ni-1)*mul_comms,size(obj.tot_pw_budget,1))) ...
							obj.quad_pw_budget(min((s+ni-1)*mul_comms,size(obj.quad_pw_budget,1)),:)];
			end
		else 
			if (mod(s, div_comms) == 0)
				ci_index = ci_index + 1;
				f_ref = obj.frplot(min(ci_index, size(obj.frplot,1)),:)';
				
				p_budget = obj.tot_pw_budget(min(ci_index, length(obj.tot_pw_budget)));
				
				for ni=1:obj.Nhzn
				obj.usum(ni,:) = [obj.tot_pw_budget(min(ci_index+ni-1,size(obj.tot_pw_budget,1))) ...
							obj.quad_pw_budget(min(ci_index+ni-1,size(obj.quad_pw_budget,1)),:)];
				end
			end
		end
		if p_budget~=pold
			pold = p_budget;
			pbc = 1;
			pw_ad_step_count = pw_ad_steps;
		else
			pbc = 0;
		end
		
		if (s~=1)
			%x = cpxplot(index+sim_mul,:)';
			T = obj.C(1:obj.Nc,:)*cpxplot(index+sim_mul,:)';
			%add noise:
			if obj.measure_noise
				nn = (rand(obj.Nc,1) - 0.5)*2 * obj.T_noise_max;
				T = T + nn;
			end	
			Tobs = cpxplot(index+sim_mul,:)';
		end
		d_i = d_is;
		pw_m = pw_ms;

		%Compute model:
		index = 1+(s-1)*sim_mul;
		[cpuplot(index+1:index+sim_mul,:), cpxplot(index+1:index+sim_mul,:), d_is, pw_ms, obj] = obj.compute_model(sim_mul, cpxplot(index,:)', V, F, d_p);	

		%new stuf
		pw_storage = 0;

		% Process Workload
		[wl, Ceff] = obj.cp_wl_process(d_i, wl);

		% Process Power Budget
		% (?)

		% Adapt Measured&Computed Power
		aup = pw_ad_aup + pw_ad_step_count * (pw_ad_achange - pw_ad_aup) / pw_ad_steps;
		adown = pw_ad_adown + pw_ad_step_count * (pw_ad_achange - pw_ad_adown) / pw_ad_steps;
		
		%pw_adapt = obj.cp_pw_adapt(pw_adapt, pw_m, pu_old, (pbc | (pw_ad_step_count>(pw_ad_steps-5))), aup, adown);
		%pw_adapt = obj.cp_pw_adapt(pw_adapt, pw_m, pu_old, pbc);
		pu_old = sum(pu); %sum(pu.*(F./F_og));
		
		[~, Tidx] = min(abs(T0v - (T-273.15)),[],2);
		% before overwriting F
		[~, Fidx] = min(abs(F0v - F),[],2);
		
		for i=1:obj.Nc
			papx(i,1) = obj.mpc_h0_app(Tidx(i), Fidx(i));
		end
		h0v = h0 + papx;
		
		% Choose Voltage
		FD = diag(f_ref)*obj.VDom;
		V = obj.cp_voltage_choice(FD);
		F = f_ref;
		if obj.ctrl_fixedv
			V = obj.V_Max*ones(obj.vd, 1);
		end

		% Compute Power
		%leak_store = obj.exp_leakage;
		%obj.exp_leakage = 0;
		%pu = obj.power_compute(F,obj.VDom*V,T,wl,d_p, 0);
		if obj.separate_omega
			pu = F .* (obj.VDom*(V .* V));
		else
			pu = F .* (obj.VDom*(V .* V)) .* Ceff;
		end
		%obj.exp_leakage = leak_store;
		
		%if (s==320)
		%	s=320;
		%end

		% MPC
		if obs
			xl = Adl_obs*xl + Bdl_obs*[pu_prev;obj.temp_amb*1000] + LK*(T - obj.C(1:obj.Nc,:)*xl);
			xlplot(s+1,:) = xl;
			rm = 1;
			rk = 1;
			rn = size(T);
			Tmpc = reshape([reshape(T,rm,[]);zeros(rk,rn(1)/rm*rn(2))],[],rn(2));
			Tmpc = [Tmpc; zeros(other_states, 1)];
			pu = obj.lin_mpc2((eye(obj.Ns)-C_obs)*xl + Tmpc, obj.temp_amb*1000, pu, obj.usum);
		else
			xlplot(s+1,:) = 0;
			ppur(s+1,:) = pu;
			%res = obj.lin_mpc2(Tobs, Ceff, obj.temp_amb*1000, pu, obj.usum);
			res = obj.lin_mpc2(Tobs, Ceff,obj.temp_amb*1000, h0v, pu, obj.usum);
			pu = res{1};
			ppumpc(s+1,:) = pu;
			tmpc(s+1, :) = res{2}';
			pceff(s+1,:) = Ceff';
		end
		if isnan(pu)
			%pu = obj.min_pw_red;
			pu = pu_prev;
			failed = failed + 1;
		%else
			%prev_pu = pu;
		end
		
		pu_min = 0.1;
		if sum(pu<pu_min)>0
			pu((pu<pu_min)) = 0.8;
			failed = failed + 1;
		end
		%cpdebug(s+1,:) = upid;
		
		% Compute Freq	
		%TODO
		%{
		alp = 0.15; %0.2995;
		F = pu.^(1 / (2*alp + 1));
		FD = diag(F)*obj.VDom;
		V = obj.cp_voltage_choice(FD);
		if obj.ctrl_fixedv
			V = obj.V_Max*ones(obj.vd, 1);
		end
		%}		
	
		%% Newton-Rapson
		obj.min_pw_red = 0.1*ones(obj.Nc,1);
		%pu = pu + (pu<obj.min_pw_red).*(obj.min_pw_red-pu);
		for vi=1:obj.vd
			fp = f_ref.*obj.VDom(:,vi);	
			cidx = [1:1:obj.Nc]' .* (fp>0);
			cidx = cidx(cidx>0);
			lim_sup = fp(fp>0);		
			lim_inf = obj.F_min*ones(sum(obj.VDom(:,vi)),1);
			pup = pu.*obj.VDom(:,vi);
			pup = pup(pup>0);
			Cip = Ceff.*obj.VDom(:,vi);
			Cip = Cip(Cip>0);
			d_pk = (d_p*obj.leak_process/1000).*obj.VDom(:,vi);
			d_pk = d_pk(d_pk>0);
			T_pk = T.*obj.VDom(:,vi);
			T_pk = T_pk(T_pk>0);
			xi = pw_function(obj, lim_inf, pup, Cip, obj.leak_vdd/1000,d_pk,T_pk);			
			xs = pw_function(obj, lim_sup, pup, Cip, obj.leak_vdd/1000,d_pk,T_pk);

			fo = zeros(length(cidx),1);

			%boundaries:
			sn = xs>0;
			so = xi>0;
			c_lim_inf = lim_inf + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);
			c_lim_sup = lim_sup + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_sup);
			%fo = fo + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);

			for nri=1:16
				fo = (c_lim_inf+c_lim_sup)/2;
				xs = pw_function(obj, fo, pup, Cip, obj.leak_vdd/1000,d_pk,T_pk);
				if (sum(abs(xs)) <= 1e-3*length(cidx)) || (sum(abs(c_lim_sup-fo))<= obj.F_discretization_step*length(cidx))
					break;
				end
				sn = xs>0;
				so = xi>0;

				xi = xi + (sn==so).*(xs-xi);
				c_lim_inf = c_lim_inf + (sn==so).*(fo-c_lim_inf);
				c_lim_sup = c_lim_sup + (sn~=so).*(fo-c_lim_sup);				
			end

			%discretize
			if obj.F_discretization_step > 0
				fo = round(fo/obj.F_discretization_step) * obj.F_discretization_step;
			end
			F(cidx) = fo;
			V(vi) = obj.FV_table(sum(max(fo) > obj.FV_table(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue
		end %for
		%}
		
		%%
		pwap(s) = pw_adapt;

		% Process Freq
		%	Check vs maxF, Temp hysteresis, etc.					
		fmaxi = obj.FV_table(sum(obj.VDom * V > ones(obj.Nc,1)*obj.FV_table(:,1)'+1e-6, 2) + 1, 3);
		F = F + (F>fmaxi).*(fmaxi - F);
		F = F + (F<obj.F_min).*(obj.F_min*ones(obj.Nc,1) - F);
		% discretize
		
		%TEST: TODO REMOVE
		pu = (Ceff.*F.*(obj.VDom*V) + obj.leak_vdd/1000).*(obj.VDom*V) + d_p*obj.leak_process/1000;

		%
		cpfplot(s+1,:) = F;
		cpvplot(s+1,:) = V;
		

	end %for loop

	% PLOTs
	t1 = obj.Ts*[0:Nsim*sim_mul]';
	%dimf = t1(end)/max(obj.Ts_input,obj.Ts_ctrl);
	%t2 = obj.Ts_input*[0:dimf]';
	t2 = obj.Ts_mpc*[0:Nsim]';

	obj.xutplot(cpxplot,cpuplot);
	obj.powerconstrplot(cpuplot);
	obj.tempconstrplot(cpxplot);
	obj.perfplot(cpfplot,obj.wl_index);
	obj.fvplot(t2, cpfplot,cpvplot);
	
	%obj.obsplot(xlplot, cpxplot);
	
	%figure();
	%plot(obj.Ts*[1:Nsim]', pwap);
	
	figure();
	%plot(obj.Ts*sim_mul*[1:Nsim]', pceff(2:end,:)*20+293, 'g'); hold on;
	plot(obj.Ts*[1:Nsim*sim_mul]', cpxplot(2:end,:) - 273, 'b'); hold on; grid on;
	plot(obj.Ts*sim_mul*[1:Nsim]', tmpc(2:end,:) - 273, 'm'); hold on;
	xlabel("Time [s]");
	ylabel("Temperature [T]");
	
	figure();
	plot(obj.Ts*sim_mul*[1:Nsim]', ppur(2:end,:), 'g'); hold on;
	plot(obj.Ts*sim_mul*[1:Nsim]', ppumpc(2:end,:), 'b'); hold on;
	
	% DATA
	disp(strcat('[MPC] number of times the optimization algorithm failed: ',int2str(failed), '/', int2str(Nsim)));
	%failed
	
	wlop = obj.wl_index / (size(obj.wrplot,3)-1) * 100;

end

function [r] = pw_function(obj, f, pu, Ci, ks1, ks2, T)
	vdd_alpha = 0.3095; %0.2995
	vdd_offset = 0.07;
	softmax_alpha = 10;
	maxfv = sum(f.*exp(softmax_alpha*f)) / sum(exp(softmax_alpha*f));
	
	if obj.separate_omega
		r = (max(f)*vdd_alpha + vdd_offset)^2.*f - pu;
	else
		r = (max(f)*vdd_alpha + vdd_offset)^2.*f.*Ci - pu;
	end
	%r = (maxfv^2*vdd_alpha^2*Ci.*f + ks1*maxfv*vdd_alpha*ones(length(f),1) + ks2) - pu;
end


