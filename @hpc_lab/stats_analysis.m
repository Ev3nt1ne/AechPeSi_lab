function res = stats_analysis(obj,ctrl,chip,simres,index)
%STATS_ANALYSIS Summary of this function goes here
%   Detailed explanation goes here

    x = simres.cpxplot;
    u = simres.cpuplot;
    f = simres.cpfplot;
    v = simres.cpvplot;
    w = simres.wlop;

	xcp = x(2:end,:)*chip.Cc' -273.15;
	ucp = u(2:end,:);
	fcp = f(2:end,:);
	vcp = v(2:end,:);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%
	% Power
	
	pwct = max(round(size(ucp,1)/max(length(obj.tot_pw_budget)-1,1)),1);
	pwtb = repelem(obj.tot_pw_budget(min(2, length(obj.tot_pw_budget)):end), pwct ,1);

    maxPbudget = chip.core_max_power*chip.Nc;
    ptime = sum(obj.tot_pw_budget < maxPbudget)*obj.Ts_target;
	
	% here it is difficult, because if I make it ceil, it will always be ok for
	% pw budget going down, but never for power budget going up. If I take
	% floor it will be the opposite. Hope that round will be ok!
	delay = round((chip.delay_F_max+chip.delay_V_max) / chip.Ts);
	sidx = delay + pwct;
	%
	pwgr = sum(ucp(1+sidx:end,:),2) - pwtb(1:end-sidx);
	pwgrP = pwgr ./ pwtb(1:end-sidx);
	
	pwgr = pwgr(pwgr > 0);
	pwgrP = pwgrP(pwgrP > 0);
	
	if isempty(pwgr)
		pwgr = 0;
	end
	if isempty(pwgrP)
		pwgrP = 0;
	end
	
	power.exMaxW = max(pwgr);
	power.exMaxP = max(pwgrP);
	power.ex95W = prctile(pwgr, 95);
	power.ex95P = prctile(pwgrP, 95);
	power.exAvW = mean( pwgr );
	power.exAvP = mean( pwgrP );
	power.exSDW = std(pwgr);
	power.exSkew = skewness(pwgr);
	power.exKurt = kurtosis(pwgr);
	%
	power.exTime = sum( sum(ucp(1+sidx:end,:),2) > pwtb(1:end-sidx) ) *chip.Ts;
    power.exTimeP = power.exTime / ptime;
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%
	% Temperature
	
	crt = chip.core_crit_temp-273.15;
	mnt = chip.core_limit_temp-273.15;
	
	temp.Max = max(xcp(:));
	temp.p95 = prctile(xcp(:),95);
	temp.Av = mean(xcp(:));
	temp.SDAv = std(mean(xcp));
	temp.SkewAv = skewness(mean(xcp));
	temp.KurtAv = kurtosis(mean(xcp));
	temp.MeanSkew = mean(skewness(xcp));
	temp.MeanKurt = mean(kurtosis(xcp));
	%
	exTC = (xcp-crt) .* (xcp > crt);
	%exTC = exTC(exTC > 0);
	if isempty(exTC)
		exTC = 0;
	end
	temp.exCr.Max = max( exTC, [], "all" );
	temp.exCr.Mean95 = mean(prctile( exTC ,95));
	temp.exCr.Mean80 = mean(prctile( exTC ,80));
	temp.exCr.MeanAv = mean( mean(exTC) );
	temp.exCr.StdAv = std(mean( exTC ));
	temp.exCr.SkewAv = skewness(mean( exTC ));
	temp.exCr.KurtAv = kurtosis(mean( exTC ));
	%
    %temp.exCr.MaxExceed = temp.exCr.Max - crt;
    %
	exTn = (xcp-mnt) .* (xcp > mnt);
	%exTn = exTn(exTn > 0);
	if isempty(exTn)
		exTn = 0;
	end
	temp.exMn.Max = max( exTn, [], "all" );
	temp.exMn.Mean95 = mean(prctile( exTn ,95));
	temp.exMn.Mean80 = mean(prctile( exTn ,80));
	temp.exMn.MeanAv = mean(mean( exTn ));
	temp.exMn.StdAv = std(mean( exTn ));
	temp.exMn.SkewAv = skewness(mean( exTn ));
	temp.exMn.KurtAv = kurtosis(mean( exTn ));
    %
    %TODO not considering variable limit
    %   I should do: max(t-mnt(t))
    %temp.exMn.MaxExceed = temp.exMn.Max - mnt;
	%
	temp.exCr.TotTime = sum(sum( xcp > crt) *chip.Ts);
	temp.exCr.StdTime = std(sum( xcp > crt) *chip.Ts);
	temp.exCr.SkewTime = skewness(sum( xcp > crt) *chip.Ts);
	temp.exCr.KurtTime = kurtosis(sum( xcp > crt) *chip.Ts);
	%
	temp.exMn.TotTime = sum(sum( xcp > mnt) *chip.Ts);
	temp.exMn.StdTime = std(sum( xcp > mnt) *chip.Ts);
	temp.exMn.SkewTime = skewness(sum( xcp > mnt) *chip.Ts);
	temp.exMn.KurtTime = kurtosis(sum( xcp > mnt) *chip.Ts);
	% no need these two: I can divide per chip.Nc
	temp.exCr.MeanTime = mean(sum( xcp > crt) *chip.Ts);
	temp.exMn.MeanTime = mean(sum( xcp > mnt) *chip.Ts);
    temp.exMn.MeanTimeP = temp.exMn.MeanTime / obj.tsim;
    temp.exCr.MeanTimeP = temp.exCr.MeanTime / obj.tsim;
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%
	% freq
	
	freq.Max = max(fcp(:));
	freq.min = min(fcp(:));
	freq.Av = mean(fcp(:));
	freq.StdAv = std(mean(fcp));
	freq.SkewAv = skewness(mean(fcp));
	freq.KurtAv = kurtosis(mean(fcp));
	freq.MeanSkew = mean(skewness(fcp));
	freq.MeanKurt = mean(kurtosis(fcp));
	
	% oscillations:
	%want to study for 1ms, 10ms, 100ms, 1s
	mdim = [1e-3 1e-2 1e-1 1] / ctrl.Ts_ctrl;
	tt = mdim > obj.tsim / ctrl.Ts_ctrl;
	mdim(tt) = obj.tsim / ctrl.Ts_ctrl;
	mdim(mdim<1) = 1;
	
	for md = 1:length(mdim)
		dd = mdim(md);
		ve = zeros(size(fcp) - [dd 0]);
		for i=1:length(fcp)-dd+1
			idx = i+dd-1;
			ve(i,:) = std(fcp(i:idx,:));
		end
	
		freq.Osc.Max(md) = max(ve(:));
		freq.Osc.MeanMax(md) = mean(max(ve));
		freq.Osc.StdMax(md) = std(max(ve));
		freq.Osc.SkewMax(md) = skewness(max(ve));
		freq.Osc.KurtMax(md) = kurtosis(max(ve));
		freq.Osc.Av(md) = mean(ve(:));
		freq.Osc.StdAv(md) = std(mean(ve));
		freq.Osc.MeanStd(md) = mean(std(ve));
		freq.Osc.SkewAv(md) = skewness(mean(ve));
		freq.Osc.KurtAv(md) = kurtosis(mean(ve));
		freq.Osc.MeanSkew(md) = mean(skewness(ve));
		freq.Osc.MeanKurt(md) = mean(kurtosis(ve));	
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%
	% Voltage
	
	vdd.Max = max(vcp(:));
	vdd.min = min(vcp(:));
	vdd.Av = mean(vcp(:));
	vdd.StdAv = std(mean(vcp));
	vdd.SkewAv = skewness(mean(vcp));
	vdd.KurtAv = kurtosis(mean(vcp));
	vdd.MeanSkew = mean(skewness(vcp));
	vdd.MeanKurt = mean(kurtosis(vcp));
	
	% oscillations:
	%want to study for 1ms, 10ms, 100ms, 1s
	mdim = [1e-3 1e-2 1e-1 1] / ctrl.Ts_ctrl;
	tt = mdim > obj.tsim / ctrl.Ts_ctrl;
	mdim(tt) = obj.tsim / ctrl.Ts_ctrl;
	mdim(mdim<1) = 1;
	
	for md = 1:length(mdim)
		dd = mdim(md);
		ve = zeros(size(vcp) - [dd 0]);
		for i=1:length(vcp)-dd+1
			idx = i+dd-1;
			ve(i,:) = std(vcp(i:idx,:));
		end
	
		vdd.Osc.Max(md) = max(ve(:));
		vdd.Osc.MeanMax(md) = mean(max(ve));
		vdd.Osc.StdMax(md) = std(max(ve));
		vdd.Osc.SkewMax(md) = skewness(max(ve));
		vdd.Osc.KurtMax(md) = kurtosis(max(ve));
		vdd.Osc.Av(md) = mean(ve(:));
		vdd.Osc.StdAv(md) = std(mean(ve));
		vdd.Osc.MeanStd(md) = mean(std(ve));
		vdd.Osc.SkewAv(md) = skewness(mean(ve));
		vdd.Osc.KurtAv(md) = kurtosis(mean(ve));
		vdd.Osc.MeanSkew(md) = mean(skewness(ve));
		vdd.Osc.MeanKurt(md) = mean(kurtosis(ve));	
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%
	% Perf
	smref = round((size(fcp,1))/(size(obj.frtrc{index}(2:end,:),1)));
	smf = round((size(obj.frtrc{index}(2:end,:),1))/(size(fcp,1)));
	
	fd = repelem(obj.frtrc{index}(2:end,:),max(smref,1),1) - repelem(fcp,max(smf,1),1);
	
	n2 = reshape(fd,[],1);
	perf.fd.l2norm = norm(n2) / sqrt(chip.Nc) / sqrt(ceil(obj.tsim / ctrl.Ts_ctrl));
	%perf.fdAv2norm = mean(norm
	perf.fd.Max = max(fd(:));
	perf.fd.min = min(fd(:));
	perf.fd.Av = mean(fd(:));
	perf.fd.Sum = sum(fd(:));
	perf.fd.StdAv = std(mean(fd));
	perf.fd.SkewAv = skewness(mean(fd));
	perf.fd.KurtAv = kurtosis(mean(fd));
	
	perf.fd.MeanSkew = mean(skewness(fd));
	perf.fd.MeanKurt = mean(kurtosis(fd));
	perf.fd.MeanStd = mean(std(fd));
	
	
	perf.wl.Max = max(w);
	perf.wl.Av = mean(w);
	perf.wl.min = min(w);
	perf.wl.Std = std(w);
	perf.wl.Skew = skewness(w);
	perf.wl.Kurt = kurtosis(w);
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%
	% out
	
	res.power = power;
	res.temp = temp;
	res.freq = freq;
	res.vdd = vdd;
	res.perf = perf;

end

