function [cpxplot, cpuplot, cpfplot, cpvplot, wlop] = launch_ncp_sim(obj, robust, show)
	if (nargin < 2) || isempty(robust)
		robust = 0;
	end
	if (nargin < 3) || isempty(show)
		show = 1;
	end

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
	pid_target = ones(obj.Nc, 1)*obj.core_limit_temp;
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
	
	obj = obj.init_compute_model(Adl_true, Bdl_true);
	
	obj.min_pw_red = 0.6;
	derT = zeros(obj.Nc,1);
	VCred = zeros(obj.Nc,1);
	Tper=8;
	pbc = 0;
	pold = p_budget;
	
	Vn = V;
	v_inc_st = 0;
	v_dec_st = 0;
	v_inc_steps = Tper+3;
	v_dec_steps = 20;
	v_inc_step_count = 0; %v_inc_steps;
	v_dec_step_count = 0; %v_dec_steps;
	v_hys_act = 0;
	v_hys_steps = 5; %15
	v_hys_step_count = 0; %v_hys_steps;
	hys_p_count = 0;
	
	pwap = zeros(Nsim, 1);
	ppap = zeros(Nsim, (6+2)*obj.vd);
	
	pw_ad_steps = 20;
	pw_ad_step_count = 0;
	pw_ad_aup = 0.015; %0.05; %0.025;
	pw_ad_adown = 0.02; %0.075; %0.025; %0.1;
	pw_ad_achange = 0.5; %0.5;
	pu_old = pw_ms;
	
	% LOOOP
	for s=1:Nsim

		%Read T, d_i, p_budget, f_ref
		if (mod(s, div_comms) == 0)
			ci_index = ci_index + 1;
			f_ref = obj.frplot(min(ci_index, size(obj.frplot,1)),:)';		
			p_budget = obj.tot_pw_budget(min(ci_index, length(obj.tot_pw_budget)));
		end
		if p_budget~=pold
			pold = p_budget;
			pbc = 1;
			hys_p_count = 0;
			v_hys_act = 0;
			v_inc_st = 0;
			v_dec_st = 0;
			v_inc_step_count = 0; %v_inc_steps;
			v_dec_step_count = 0; %v_dec_steps;
			v_hys_step_count = 0;
			pw_ad_step_count = pw_ad_steps;
		else
			pbc = 0;
			pw_ad_step_count = pw_ad_step_count - 1;
			pw_ad_step_count = (pw_ad_step_count<0)*0 + (pw_ad_step_count>=0)*pw_ad_step_count;
		end
		if (s~=1)
			%x = cpxplot(index+sim_mul,:)';
			T = obj.C(1:obj.Nc,:)*cpxplot(index+sim_mul,:)';
			%add noise:
			if obj.measure_noise
				nn = (rand(obj.Nc,1) - 0.5)*2 * obj.T_noise_max;
				T = T + nn;
			end			
			if (s>Tper+1) && (mod(s,Tper)==2)
				derT = T - obj.C(1:obj.Nc,:)*cpxplot(index-sim_mul*Tper+1,:)';
			end
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
		
		pw_adapt = obj.cp_pw_adapt(pw_adapt, pw_m, pu_old, (pbc | (pw_ad_step_count>(pw_ad_steps-5))), aup, adown);
		pu_old = sum(pu); %sum(pu.*(F./F_og));
		
		% Choose Voltage
		FD = diag(f_ref)*obj.VDom;			
		V = obj.cp_voltage_choice(FD);
		F = f_ref;
		
		%if s==300 %1300
		%	s=s;
		%end
		
		% Control Temperature:
		if (mod(s,Tper)==2)
			Tband = [45.0 65.0 80.0 85 ];
			Tband = Tband + 273.15;
			derBand = [0 0.5 1.0 2.0];
			derBand = derBand + obj.measure_noise*obj.T_noise_max/2;
			%{
			redmat =	[1 1 1 0 -1;
						 1 1 0 0 -1;
						 0 0 0 -1 -2;
						 0 0 -1 -2 -3;
						 0 -1 -2 -3 -4];
			%}
			redmat =	[2	2	1	0	-1;
						 2	1	1	0	-1;
						 1	1	0	-1	-2;
						 1	0	-1	-2	-3;
						 0	0	-2	-3	-4];

			tr(:,2) = sum(T > Tband,2);
			tr(:,1) = sum(derT > derBand,2);
			tr = tr + 1; %indexing

			VCred = VCred + redmat((tr(:,2)-1)*size(redmat,1) + tr(:,1));
			%saturate: (no turbo atm)
			VCred = VCred + (VCred > 0).*(0-VCred);
		end
		
		VDred = min(diag(VCred)*obj.VDom);
		V = V + fix(VDred/2)'*0.05;		
		
		%F = F + obj.VDom*(VDred.*mod(VDred,2))'*0.05;
		
		%cap F:
		% maybe remove?
		fmaxi = obj.FV_table(sum(obj.VDom * V > ones(obj.Nc,1)*obj.FV_table(:,1)', 2) + 1, 3);
		F = F + (F>fmaxi).*(fmaxi - F);
		
		F = F + obj.VDom*(mod(VDred,2).*sign(VDred))'*0.05;

		% Compute Power
		%leak_store = obj.exp_leakage;
		%obj.exp_leakage = 0;
		%pu = obj.power_compute(F,obj.VDom*V,T,wl,d_p, 0);
		%obj.exp_leakage = leak_store;
		pu = (Ceff.*F.*(obj.VDom*V) + obj.leak_vdd/1000).*(obj.VDom*V) + d_p*obj.leak_process/1000;
				
		% DispatchPower
		%TODO, quad power budget dispatching
		delta_p = sum(pu) - p_budget + pw_adapt;
		if (delta_p > 0)
			dompw = obj.VDom'*pu;
			domT = ((obj.VDom'*T ./ sum(obj.VDom)') + max(T.*obj.VDom)' ) / 2; %here mean*2 / 3??
			[resdpw, ~] = obj.cp_pw_dispatcher(domT, pid_target(1:obj.vd), delta_p, dompw); %todo: pid_target(1:obj.vd)
			deltapd = dompw-resdpw;
			for vdi=1:obj.vd
				td = T.*obj.VDom(:,vdi);
				td = td(td>0);
				pd = pu.*obj.VDom(:,vdi);
				pd = pd(pd>0);
				cd = Ceff.*obj.VDom(:,vdi);
				cd = cd(cd>0);
				cidx = [1:1:obj.Nc]' .* obj.VDom(:,vdi);
				cidx = cidx(cidx>0);
				[pd, ~] = obj.cp_pw_dispatcher_c(cd,deltapd(vdi), pd);
				pu(cidx) = pd;
			end
			%[pu, pws] = obj.cp_pw_dispatcher(T, pid_target, delta_p, pu);
			%pw_storage = pw_storage + pws;	
		end
		
		if sum(isempty(pu)) || sum(pu<=0)
			disp("error, something wrong");
			[cpuplot(index+1:index+sim_mul,:), cpxplot(index+1:index+sim_mul,:), d_is, pw_ms, obj] = obj.compute_model(sim_mul, cpxplot(index,:)', V, F, d_p);	
			break;
		end
		
		% Compute Freq
		%% Newton-Rapson		
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
			xi = pw_function(lim_inf, pup, Cip, obj.leak_vdd/1000,d_pk);			
			xs = pw_function(lim_sup, pup, Cip, obj.leak_vdd/1000,d_pk);
			
			fo = zeros(length(cidx),1);
			
			%boundaries:
			sn = xs>0;
			so = xi>0;
			c_lim_inf = lim_inf + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);
			c_lim_sup = lim_sup + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_sup);
			%fo = fo + (sn==so).*(sn.*lim_inf + ~sn.*lim_sup - lim_inf);
			
			for nri=1:16
				fo = (c_lim_inf+c_lim_sup)/2;
				xs = pw_function(fo, pup, Cip, obj.leak_vdd/1000,d_pk);
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
		end		
		
		%% HYSTERESIS V-FILTER ACTIVATION
		V_inc = (V > (Vn+1e-3));
		V_dec = (V < (Vn-1e-3));
		
		v_inc_st = v_inc_st | V_inc;
		v_inc_step_count = (v_inc_steps.*V_inc) + ((~V_inc).*(v_inc_step_count - v_inc_st));
		
		v_inc_st = (v_inc_step_count>0); % | v_dec_st;
		
		v_dec_st = v_inc_st & ( v_dec_st | V_dec);
		%if v_dec_st
		%	s=s;
		%end
		%v_dec_step_count = (v_dec_steps.*V_dec) + ((~V_dec).*(v_dec_step_count - v_dec_st));
		v_dec_step_count = ((v_dec_step_count>0).*(v_dec_step_count-1)) + (v_dec_step_count<=0).*(v_dec_st.*v_dec_steps);
		
		v_dec_st = (v_dec_step_count>0);
		
		v_hys_act_hold = v_hys_act;
		
		v_hys_act = v_hys_act | (v_dec_st & V_inc);
		
		%Initialize
		v_hys_step_count = ((v_hys_step_count>0).*v_hys_step_count) + (v_hys_step_count<=0).*(v_hys_act.*v_hys_steps);
		%Reduce
		%pull = (Ceff.*F.*(obj.VDom*V) + obj.leak_vdd/1000).*(obj.VDom*V) + d_p*obj.leak_process/1000;
		v_hys_red = (pw_m - p_budget) < -(p_budget*(0.125+hys_p_count/10));
		%(sum(pull) - p_budget + pw_adapt) < -6 ;
		%(delta_p < (p_budget - 6)); %TODO: 6 depends on the #of cores, #of domains, and the ratio between the two.
		%TODO: also depends on the power budget itself, and the PMAX, PMIN
		v_hys_step_count = (v_hys_act & v_hys_step_count) .* ...
			(v_hys_step_count - v_hys_red);
		
		v_hys_act = (v_hys_step_count>0);
		
		hys_p_count = hys_p_count + (v_hys_act_hold ~= v_hys_act) .* v_hys_act;
		
		%if v_hys_act
		%	s=s;
		%end
		
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

		
		V = V.*((~v_hys_act)|(V<=Vn)) + Vn.*(v_hys_act & (V>Vn));
		
		Vn = V;
		
		%%
		pwap(s) = pw_adapt;

		% Save OG freq
		%F_og = F;

		% Process Freq
		%	Check vs maxF, Temp hysteresis, etc.
		fmaxi = obj.FV_table(sum(obj.VDom * V > ones(obj.Nc,1)*obj.FV_table(:,1)'+1e-6, 2) + 1, 3);
		F = F + (F>fmaxi).*(fmaxi - F);
		F = F + (F<obj.F_min).*(obj.F_min*ones(obj.Nc,1) - F);
		
		%TEST: TODO REMOVE
		pu = (Ceff.*F.*(obj.VDom*V) + obj.leak_vdd/1000).*(obj.VDom*V) + d_p*obj.leak_process/1000;

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
		
		%figure();
		%plot(obj.Ts*[1:Nsim]', pwap);
		
		%{
		figure();
		subplot(2,3,1);
		plot(ppap(:,[1:obj.vd*2]));hold on;
		subplot(2,3,2);
		plot(ppap(:,[obj.vd*2+1:obj.vd*4]));hold on;
		subplot(2,3,3);
		plot(ppap(:,[obj.vd*4+1:obj.vd*6]));hold on;
		subplot(2,3,4);
		plot(ppap(:,[obj.vd*6+1:obj.vd*7]));hold on;
		subplot(2,3,5);
		plot(ppap(:,[obj.vd*7+1:obj.vd*8]));hold on;
		%}
		
	end
	

	
	wlop = obj.wl_index / (size(obj.wrplot,3)-1) * 100;
end

function [r] = pw_function(f, pu, Ci, ks1, ks2)
	vdd_alpha = 0.3095; %0.2995
	vdd_offset = 0.07;
	softmax_alpha = 10;
	maxfv = sum(f.*exp(softmax_alpha*f)) / sum(exp(softmax_alpha*f));
	
	r = ( (max(f)*vdd_alpha + vdd_offset)^2*Ci.*f + ks1*(max(f)*vdd_alpha+vdd_offset)*ones(length(f),1) + ks2) - pu;
	%r = (maxfv^2*vdd_alpha^2*Ci.*f + ks1*maxfv*vdd_alpha*ones(length(f),1) + ks2) - pu;
end
