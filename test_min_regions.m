	V = [0.5:0.1:1.1];
	d_p = 1;
	T = [20:5:90];
	%T = [20:5:90];
	
	Power_static = V*hpc.leak_vdd_k + d_p*hpc.leak_process_k;
	
	res = [];
	plane = [];
	for i=1:length(T)
		res(:,i) = (Power_static .* exp(hpc.leak_exp_vdd_k*V + T(i)*hpc.leak_exp_t_k + hpc.leak_exp_k))';
		plane(:,i) = Power_static;
	end
	
	F = hpc.FV_table( [sum(V>hpc.FV_table(:,1))+1]',3);
	
	ci = hpc.dyn_ceff_k;
	Ceff_low = [1 zeros(1,hpc.ipl-1)] * ci';
	Ceff_high = [zeros(1,hpc.ipl-1) 1] * ci';
	Power_dyn_low = Ceff_low .* F .* (V' .* V');
	Power_dyn_high = Ceff_high .* F .* (V' .* V');

	[T,C] = create_regions(V, T, res, 4, 1)

	xlabel("Voltage [V]")
	ylabel("Temperature [C]")

	%%

	[k0, k1, k2] = hpc.pws_ls_approx([0.5 1.2], [20 90], 0.9, 0.2355, -0.1, 1)

