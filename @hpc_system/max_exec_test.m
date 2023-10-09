function [wlop] = max_exec_test(obj)
	
	d_is = [ones(obj.Nc,1) zeros(obj.Nc,obj.ipl-1)];
	wl = d_is;
	d_p = ones(obj.Nc,1);
	
	T = obj.C(1:obj.Nc,:)*obj.x_init;
	x = obj.x_init + (rand(obj.Ns,1) - 0.5*ones(obj.Ns,1));
	sim_mul = ceil(obj.Ts_ctrl/obj.Ts);
	Nsim = ceil(obj.tsim / obj.Ts_ctrl);
	div_comms = ceil(obj.Ts_input / obj.Ts_ctrl);
	ci_index = 1;
	f_ref = obj.frplot(ci_index,:)';
	F = obj.F_min*ones(obj.Nc,1);
	V = obj.V_min*ones(obj.vd,1);
	Adl_true = obj.Ad_true;
	Bdl_true = obj.Bd_true;
	
	obj = obj.init_compute_model(Adl_true, Bdl_true);

	cpxplot = ones(sim_mul+1, obj.Ns) .* x';
	
	% LOOOP
	for s=1:Nsim

		%Read T, d_i, p_budget, f_ref
		if (mod(s, div_comms) == 0)
			ci_index = ci_index + 1;
			f_ref = obj.frplot(min(ci_index, size(obj.frplot,1)),:)';		
		end

		%Compute model:
		[~, cpxplot(1:sim_mul,:), ~, ~, obj] = obj.compute_model(sim_mul, cpxplot(sim_mul,:)', V, F, d_p);	

		FD = diag(f_ref)*obj.VDom;			
		V = obj.cp_voltage_choice(FD);
		F = f_ref;

	end %for loop

	% PLOTs	
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


