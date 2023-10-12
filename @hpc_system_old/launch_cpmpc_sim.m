function [cpxplot,cpuplot, cpfplot, cpvplot, xlplot, wlop] = launch_cpmpc_sim(obj, obs)

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
	F_MA = zeros(obj.Nc,1);

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
	F_og = F;

	pw_adapt = 0;
	pu = obj.min_pw_red;
	pw_ms = sum(pu);
	pu_prev = pu;

	obj = obj.lin_mpc_setup();						
	failed = 0;
	
	% Observer
	Adl_obs = obj.Ad_obs;
	Bdl_obs = obj.Bd_obs;
	poles = ones(obj.Ns,1);
	poles(1:2:obj.Ns-other_states) = obj.Obs_poles(1);
	poles(2:2:obj.Ns-other_states) = obj.Obs_poles(2);
	poles(end-other_states+1:end) = 0.2;
	%
	LK = place(Adl_obs', obj.C(1:obj.Nc,:)', poles)';
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
	pbc = 0;
	
	% wl stuff
	%quantum_storage = obj.quantum_instr*ones(obj.Nc, 1);
	%wl_index = ones(obj.Nc, 1);
	
	obj = obj.init_compute_model(Adl_true, Bdl_true);
	
	% LOOOP
	for s=1:Nsim

		%Read T, d_i, p_budget, f_ref
		if mul_comms > 0
			f_ref = obj.frplot(min(s*mul_comms,size(obj.frplot,1)),:)';
			
			pold = p_budget;
			p_budget = obj.tot_pw_budget(min(s*mul_comms, length(obj.tot_pw_budget)));
			if p_budget~=pold
				pbc = 1;
			else
				pbc = 0;
			end
			
			for ni=1:obj.Nhzn
				obj.usum(ni,:) = [obj.tot_pw_budget(min((s+ni-1)*mul_comms,size(obj.tot_pw_budget,1))) ...
							obj.quad_pw_budget(min((s+ni-1)*mul_comms,size(obj.quad_pw_budget,1)),:)];
			end
		else 
			if (mod(s, div_comms) == 0)
				ci_index = ci_index + 1;
				f_ref = obj.frplot(min(ci_index, size(obj.frplot,1)),:)';
				
				pold = p_budget;
				p_budget = obj.tot_pw_budget(min(ci_index, length(obj.tot_pw_budget)));
				if p_budget~=pold
					pbc = 1;
				else
					pbc = 0;
				end
				
				for ni=1:obj.Nhzn
				obj.usum(ni,:) = [obj.tot_pw_budget(min(ci_index+ni-1,size(obj.tot_pw_budget,1))) ...
							obj.quad_pw_budget(min(ci_index+ni-1,size(obj.quad_pw_budget,1)),:)];
				end
			end
		end
		
		if (s~=1)
			%x = cpxplot(index+sim_mul,:)';
			T = obj.C(1:obj.Nc,:)*cpxplot(index+sim_mul,:)';
			Tobs = cpxplot(index+sim_mul,:)';
		end
		d_i = d_is;
		pw_m = pw_ms;
		
		if mod(s,100) == 0
			lalla = 1;
		end

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
		pw_adapt = obj.cp_pw_adapt(pw_adapt, pw_m, sum(pu_prev.*(F./F_og)), pbc);

		% Choose Voltage
		FD = diag(f_ref - F_MA)*obj.VDom;
		V = obj.cp_voltage_choice(FD);
		F = f_ref;
		if obj.ctrl_fixedv
			V = obj.V_Max*ones(obj.vd, 1);
		end

		% Compute Power
		%leak_store = obj.exp_leakage;
		%obj.exp_leakage = 0;
		pu = obj.power_compute(F,obj.VDom*V,T,wl,d_p, 0);
		%obj.exp_leakage = leak_store;

		% DispatchPower
		%TODO, quad power budget dispatching
		%delta_p = sum(pu) - obj.tot_pw_budget + pw_adapt;
		%if (delta_p > 0)
		%	[pu, pws] = obj.cp_pw_dispatcher(T, obj.core_crit_temp, delta_p, pu);
		%	pw_storage = pw_storage + pws;
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
			pu = obj.lin_mpc((eye(obj.Ns)-C_obs)*xl + Tmpc, obj.temp_amb*1000, pu, obj.usum);
		else
			xlplot(s+1,:) = 0;
			pu = obj.lin_mpc(Tobs, obj.temp_amb*1000, pu, obj.usum);
		end
		if isnan(pu)
			%pu = obj.min_pw_red;
			pu = pu_prev;
			failed = failed + 1;
		%else
			%prev_pu = pu;
		end
		%cpdebug(s+1,:) = upid;

		% Compute Freq
		%Ceff = d_i * (obj.ceff_pw_coeff / 1000)';
		if obj.exp_leakage
			maxT = 125+273.15;
			pex = exp((obj.VDom*V*obj.exp_leak_coeff(1) + (min(T,ones(length(T),1)*maxT)-273.15)*obj.exp_leak_coeff(2) - ones(length(T),1)*obj.exp_leak_coeff(3))/1000);
		else
			pex = 1;
		end
		
		if ~obj.iterative_fv
			F = (((pu - d_p*obj.leak_process/1000) ./ (obj.VDom*V) ./ pex) - obj.leak_vdd/1000) ./ (obj.VDom*V) ./ Ceff;
		else			
			%% Newton-Rapson
			obj.min_pw_red = 0.7*ones(obj.Nc,1);
			pu = pu + (pu<obj.min_pw_red).*(obj.min_pw_red-pu);
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
					if (sum(abs(xs)) <= 5/obj.vd) || (sum(abs(c_lim_sup-xs))<= obj.F_discretization_step*length(cidx))
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
			
		end %if iterative_fv
		
		%%

		% Save OG freq
		F_og = F;

		% Process Freq
		%	Check vs maxF, Temp hysteresis, etc.					
		fmaxi = F;
		if ~obj.ctrl_fixedv
			fmaxi = obj.FV_table(sum(obj.VDom * V > ones(obj.Nc,1)*obj.FV_table(:,1)'+1e-6, 2) + 1, 3);
		else
			fmaxi = f_ref+0.001;
		end
		F = F + (F>fmaxi).*(fmaxi - F);
		F = F + (F<obj.F_min).*(obj.F_min*ones(obj.Nc,1) - F);
		% discretize
		if obj.F_discretization_step > 0
			F = fix(F/obj.F_discretization_step) * obj.F_discretization_step;
		end

		alpha_MA =  0.1; %0.01;
		if obj.ctrl_MA
			%if powerbudget has changed || freq changed
			% interpolate a parabola
			%	parabola depends on the changes (powe >> freq, and Delta of changes)
			% f_MA = f_ref - f_app --> sat 0
			%	Sat 0 is a choice, without we have TurboBoost!
			%TODO Turbo boosting is not working!
			turbo_boost = zeros(obj.Nc,1);
			tt = (f_ref - F_og)>-turbo_boost;
			F_MA = F_MA*(1-alpha_MA) + alpha_MA*(tt.*(f_ref - F_og) + ~tt.*0); %This is 0 and not turbo_boost!
		end
		
		pu_prev = pu;

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

	obj.obsplot(xlplot, cpxplot);
	
	obj.fvplot(t2, cpfplot,cpvplot);
	
	% DATA
	disp(strcat('[MPC] number of times the optimization algorithm failed: ',int2str(failed), '/', int2str(Nsim)));
	%failed
	
	wlop = obj.wl_index / (size(obj.wrplot,3)-1) * 100;

end

function [r] = pw_function(obj, f, pu, Ci, ks1, ks2, T)
	vdd_alpha = 0.3095; %0.2995
	softmax_alpha = 10;
	maxfv = sum(f.*exp(softmax_alpha*f)) / sum(exp(softmax_alpha*f));
	
	if obj.exp_leakage
		maxT = 125+273.15;
		pex = exp((max(f)*vdd_alpha*obj.exp_leak_coeff(1) + (min(T,ones(length(T),1)*maxT)-273.15)*obj.exp_leak_coeff(2) - ones(length(T),1)*obj.exp_leak_coeff(3))/1000);
	else
		pex = 1;
	end
	
	r = (max(f)^2*vdd_alpha^2*Ci.*f + ks1*max(f)*vdd_alpha*ones(length(f),1).*pex + ks2) - pu;
	%r = (maxfv^2*vdd_alpha^2*Ci.*f + ks1*maxfv*vdd_alpha*ones(length(f),1) + ks2) - pu;
end

