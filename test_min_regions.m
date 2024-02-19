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



	%%













	%%

	disp("Power, Exceed: 95p(W) -- Av(W) -- Time(s)");
	a = [res.power.ex95W ...
		res.power.exAvW	 ...
		res.power.exTime];
	disp(a);

	disp("Temp: Max -- p95 -- Av");
	a = [res.temp.Max	...
		res.temp.p95	...
		res.temp.Av];
	disp(a);

	disp("Exceed Crit: Mean95p -- Mean80p -- MeanAv -- TotTime");
	a = [res.temp.exCr.Mean95	...
		res.temp.exCr.Mean80	...
		res.temp.exCr.MeanAv	...
		res.temp.exCr.TotTime];
	disp(a);

	disp("Exceed Control Limit: Mean95p -- Mean80p -- MeanAv -- TotTime");
	a = [res.temp.exMn.Mean95 ...
		res.temp.exMn.Mean80  ...
		res.temp.exMn.MeanAv  ...
		res.temp.exMn.TotTime];
	disp(a);

	disp("freq: Av -- Std Av");
	a = [res.freq.Av ...
		res.freq.StdAv];
	disp(a);

	disp("Vdd: Av -- Std Av");
	a = [res.vdd.Av ...
		res.vdd.StdAv];
	disp(a);

	disp('Perf fd: l2norm -- Av -- StdAv');
	a = [res.perf.fd.l2norm	...
	res.perf.fd.Av			...
	res.perf.fd.StdAv];
	disp(a);
	
	disp("Perf wl: Max -- Av -- min -- Std -- Skew -- Kurt");
	a = [res.perf.wl.Max	...
	res.perf.wl.Av			...
	res.perf.wl.min			...
	res.perf.wl.Std			...
	res.perf.wl.Skew		...
	res.perf.wl.Kurt];
	disp(a);




