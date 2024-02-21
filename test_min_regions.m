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


	%%

	V = [0.5:0.05:1.2];

	alp = 0.3716; %0.3095;
	alp0 = 0.306; %0.01;
	k1 = 1;
	k2 = 0;

	alpV0 = mean((V - alp0).^(1/alp) - V.^(1/alp))

	F = V.^alp + alp0;

	%{
	fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-20],...
               'Upper',[20],...
               'StartPoint',[1]);
	%}

	%% (ORIGINAL + Offset)
	f = fit(hpc.FV_table(:,3), hpc.FV_table(:,1), fittype("x^a+b"))
	alp = 0.4768;
	alp0 = -0.7277;
	k1 = 1;
	k2 = 0;

	%% (ORIGINAL + Offset + k)
	f = fit(hpc.FV_table(:,3), hpc.FV_table(:,1), fittype("c*x^a+b"))
	% not working

	%% Original (THIS IS JUST WRONG!)
	f = fit(hpc.FV_table(:,3), hpc.FV_table(:,1), fittype("x^a"))
	alp = 0.2502;
	alp0 = 0;
	k1 = 1;
	k2 = 0;

	%% Voltage
	f = fit(hpc.FV_table(:,1), hpc.FV_table(:,3), fittype("x^a+b"))
	alp = 3.497;
	alp0 = 1.93;
	k1 = 1;
	k2 = 0;
	
	%%
	f = fit(hpc.FV_table(:,3), hpc.FV_table(:,1), 'exp1')
	k1 = 0.306;
	k2 = 0;
	alp = 0.3716;
	alp0 = 0;

	%%
	f = fit(hpc.FV_table(:,3), hpc.FV_table(:,1), 'poly1')
	k1 = 0.2995;
	k2 = 0;
	alp = 1;
	alp0 = 0.05411;

	%%
	f = fit(hpc.FV_table(:,3), hpc.FV_table(:,1), 'poly2')
	k1 = 0.05337;
	k2 = 0.02734;
	alp = 2;
	alp0 = 0.3731;

	%% Voltage
	f = fit(hpc.FV_table(:,1), hpc.FV_table(:,3), 'poly2')
	k1 = -1.979;
	k2 = 6.659;
	alp = 2;
	alp0 = -1.48;

	%%

	figure();
	plot(hpc.FV_table(:,1), hpc.FV_table(:,3));
	hold on; grid on;
	plot(k2*hpc.FV_table(:,3) + k1*hpc.FV_table(:,3).^alp + alp0, hpc.FV_table(:,3));
	%plot(hpc.FV_table(:,3).^alp + alp0, hpc.FV_table(:,3));
	%plot(hpc.FV_table(:,1), hpc.FV_table(:,1).^(1/alp) + alpV0);

	%%
	figure();
	plot(hpc.FV_table(:,1), hpc.FV_table(:,3));
	hold on; grid on;
	%plot(hpc.FV_table(:,1), hpc.FV_table(:,1).^alp + alp0);
	plot(hpc.FV_table(:,1), k2*hpc.FV_table(:,1) + k1*hpc.FV_table(:,1).^alp + alp0);
	%plot(hpc.FV_table(:,1), hpc.FV_table(:,1).^(1/alp) + alpV0);




	%% Analysis of FV^2

	figure()
	plot(hpc.FV_table(:,1), hpc.FV_table(:,3).*hpc.FV_table(:,1).*hpc.FV_table(:,1));
	title("FV^2 wrt V");

	figure()
	plot(hpc.FV_table(:,3), hpc.FV_table(:,3).*hpc.FV_table(:,1).*hpc.FV_table(:,1));
	title("FV^2 wrt F");



	%%
	
	iter = [ ];
	for i=1:length(ctrl.solver_stats)
		iter = [iter; ctrl.solver_stats(i).solveroutput.info.iter];
	end

	

	






