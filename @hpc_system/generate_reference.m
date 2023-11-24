function obj = generate_reference(obj, show, alpha_u, alpha_f, alpha_w)
	if (nargin < 5) || isempty(alpha_w)
		alpha_w = 1;
	end
	if (nargin < 4) || isempty(alpha_f)
		alpha_f = 0.5; %0.02;
	end
	if (nargin < 3) || isempty(alpha_u)
		alpha_u = 0.05;
	end
	if (nargin < 2) || isempty(show)
		show = 0;
	end

	%Tsm = min(obj.Ts_ctrl, obj.Ts_mpc);
	inp_mul = ceil(obj.Ts_input/obj.Ts);
	Ninputs = ceil(obj.tsim / obj.Ts_input);
	wl_dim = length(obj.ceff_pw_coeff);
	freq_min_period = 0.25;
	freq_mul = round(freq_min_period/obj.Ts_input);
	%freq_mul = ceil(Ninputs/pp);

	%init
	%obj.urplot = zeros(Ninputs+1, obj.Ni_c);
	%obj.zrplot = zeros(Ninputs*inp_mul+1, obj.Ni_c);
	obj.frplot = zeros(Ninputs+1, obj.Nc);
	%obj.wrplot = zeros(obj.Ni_c, wl_dim, Ninputs*inp_mul);

	%noise
	%Covr = chol(obj.pn_sigma);
	%z = 0; %repmat(obj.pn_mu,obj.Ni_c,1) + randn(obj.Ni_c,1)*Covr;
	%obj.zrplot(1,:) = z;

	%uref
	%u_ref = (obj.core_Max_power-0.5) * ones(obj.Ni_c,1);
	%obj.urplot(1,:) = u_ref;

	%fref
	f_ref = ones(obj.Nc,1)*(obj.F_Max/5*4);
	obj.frplot(1,:) = f_ref;
	alpha_f1 = alpha_f;
	alpha_f2 = alpha_f / 10;

	%% OTHERS
	
	for s=1:Ninputs
		%{
		for ssim=1:inp_mul
			index = (s-1)*inp_mul + ssim;
			z = z*(1-obj.pn_alpha) + obj.pn_alpha*(repmat(obj.pn_mu,obj.Ni_c,1) + randn(obj.Ni_c,1)*Covr);

			obj.zrplot(index+1,:) = z;			
		end

		%Power
		changing_u = randn(obj.Ni_c)+1.5*ones(obj.Ni_c);
		u_ref = u_ref*(1-alpha_u) + alpha_u*changing_u(:,1)*obj.core_Max_power/2;

		%saturate:
		u_ref = u_ref + (u_ref < obj.core_min_power).*(-u_ref+obj.core_min_power);
		u_ref = u_ref + (u_ref > obj.core_Max_power).*(-u_ref+obj.core_Max_power);

		obj.urplot(s+1,:) = u_ref;
		%}
		if mod(s,freq_mul) == 0
			%Freq
			changing_f = (rand(obj.Nc,1)*1.5)*(obj.F_Max-obj.F_min) + obj.F_min*ones(obj.Nc,1);
			change = rand(obj.Nc, 1) - (1-(1/obj.freq_fact)-0.5);
			alpha_f = (change > 0.51)*alpha_f1 + ((change > 0.31)&(change <= 0.51))*alpha_f2;
			f_ref = f_ref.*(1-alpha_f) + alpha_f.*changing_f;

			%saturate:
			f_ref = f_ref + (f_ref < obj.F_min).*(-f_ref+obj.F_min);
			f_ref = f_ref + (f_ref > obj.F_Max).*(-f_ref+obj.F_Max);
		end

		obj.frplot(s+1,:) = f_ref;
	end	%for

	%% wlref
	
	%ll = Ninputs*inp_mul;

	% Init
	obj.quantum_us = max(obj.Ts*1e6, obj.quantum_us);
	wl_min_exec_qt = ceil(obj.wl_min_exec_us / obj.quantum_us);
	wl_mean_exec_qt = ceil(obj.wl_mean_exec_us / obj.quantum_us) - wl_min_exec_qt;
	wl_mean_exec_qt(wl_mean_exec_qt==0) = 1;
	ll = ceil(obj.tsim*1e6/obj.quantum_us);
	obj.wrplot = zeros(obj.Ni_c, wl_dim, ll);
	%check on normalization
	obj.wl_prob = obj.wl_prob/sum(obj.wl_prob);

	% create distribution graph
	cov = wl_mean_exec_qt/3; %sqrt(input)*1.8;
	x={};
	y={};
	for i=1:obj.ipl
		x{i}(:) = 1:1:wl_mean_exec_qt(i)*2;
		yn = exp(-(x{i}-wl_mean_exec_qt(i)).^2 ./ (2*cov(i).^2));
		y{i}(:) = yn ./ sum(yn);
		%figure(), plot(x{i},y{i})
	end

	ppwl = zeros(obj.Nc, ll);
	% I tried, but it's hard to make it vectorial for each core
	for c=1:obj.Nc
	%init
	new_wl = zeros(1,obj.ipl);
	primary_wl = 1;
	wi = 1;
	while wi<=ll

		% Add min execution requirements		
		for nl=1:wl_min_exec_qt(primary_wl)		
			%call function for wls
			sp = obj.wl_prob_compute(primary_wl, new_wl);
			
			% Update new_wl
			new_wl = new_wl + ((sp>0) .* (1-(new_wl>0))).*wl_min_exec_qt;
			new_wl = new_wl - (new_wl>0);
			new_wl = new_wl + (new_wl<0).*(-new_wl);
			
			% Add
			obj.wrplot(c,:,wi) = sp;

			if show
				ppwl(c,wi) = primary_wl;
			end

			% Keep track of wi
			wi = wi + 1;
			if (wi > ll)
				break;
			end
			
			%clc;
			
		end %internal for
		if (wi > ll)
			break;
		end

		%%%%%%%%%

		stop = 0;
		counter = 0;
		% do not clear new_wl
		while stop==0		
			%call function for secondary wls
			sp = obj.wl_prob_compute(primary_wl, new_wl);
			
			% Update new_wl
			new_wl = new_wl + ((sp>0) .* (1-(new_wl>0))).*wl_min_exec_qt;
			new_wl = new_wl - (new_wl>0);	
			new_wl = new_wl + (new_wl<0).*(-new_wl);
			
			% Add
			obj.wrplot(c,:,wi) = sp;

			if show
				ppwl(c,wi) = primary_wl;
			end

			% Keep track of wi
			wi = wi + 1;
			if (wi > ll)
				break;
			end
			
			% get out check
			cv = rand(1);
			counter = counter + 1;
			p = 1 - y{primary_wl}(counter)/(1-sum(y{primary_wl}(1:counter-1)));
			if counter == wl_mean_exec_qt(primary_wl)*2
				p = 0;
			end
			stop = cv > p;
			
			%clc;
			
		end %while stop==0
		if (wi > ll)
			break;
		end

		% Choose new primary workload
		% do not clear new_wl
		sp = rand(1,obj.ipl);
		ttsp = (sp>(1-obj.wl_prob));
		if sum(ttsp) > 0
			sp = ttsp .* sp;
		else
			sp = sp .* obj.wl_prob;
		end
		primary_wl = sum((sp==max(sp)) .* [1:1:obj.ipl]);

	end %external while
	end %for core
	
	
	%% Show PLOTS
	if show
		%{
		figure('Name', 'Show Power Reference');
		ax1=subplot(2,1,1);
		plot(obj.Ts*inp_mul*[0:Ninputs]', obj.urplot);
		grid on,xlabel('Time [s]'),ylabel('Input Ref [W]'), hold on;
			yline(obj.core_Max_power, '--', 'Max Power', 'LineWidth',1,  'Color', 'b');
			yline(obj.core_min_power, '--','min Power', 'LineWidth',1,  'Color', 'b');

		ax2=subplot(2,1,2);
		plot(obj.Ts*[0:Ninputs*inp_mul]', obj.zrplot);
		grid on,xlabel('Time [s]'),ylabel('Input Noise [W]'), hold on;

		linkaxes([ax1,ax2],'x');
		%}

		%

		figure('Name', 'Show Freq Reference');
		plot(obj.Ts*inp_mul*[0:Ninputs]', obj.frplot);
		grid on,xlabel('Time [s]'),ylabel('Input Ref [GHz]'), hold on;
			yline(obj.F_Max, '--', 'Max Freq', 'LineWidth',1,  'Color', 'b');
			yline(obj.F_min, '--','min Freq', 'LineWidth',1,  'Color', 'b');

		%

		figure('Name', 'Show Workload Reference');
		axx = [];
		for i=1:obj.ipl
			hold on;
			axx = [axx subplot(obj.ipl+1, 1, i)];
			plot([1:1:ll]'*obj.quantum_us*1e-6,squeeze(obj.wrplot(:,i,:)));
		end
		axx = [axx subplot(obj.ipl+1, 1, obj.ipl+1)];
		plot([1:1:ll]'*obj.quantum_us*1e-6,ppwl'); grid on;
		xlabel("Time [s]");
		linkaxes(axx, 'x');
		
		npp = 4; %number per plot
		ppf = 3; %plot per figure
		nf = ceil(obj.Nc/(npp*ppf)); %number of figures
		axx2 = [];
		stop = 0;
		for fi=1:nf
			figure('Name', strcat("Wl Reference core: ", int2str((fi-1)*(npp*ppf)+1), " - ", int2str((fi-1)*(npp*ppf)+(npp*ppf)) ) );
			hold on;
			for sp=1:ppf
				idx = (fi-1)*(npp*ppf) + (sp-1)*npp +1;
				idxf = idx+npp-1;
				if idxf >= obj.Nc
					stop = 1;
					idxf = obj.Nc;
				end
				for pp=1:obj.ipl
					hold on;
					axx2 = [axx2 subplot(obj.ipl+1, ppf, (pp-1)*ppf + sp)];
					plot([1:1:ll]'*obj.quantum_us*1e-6,squeeze(obj.wrplot(idx:idxf,pp,:)));
					xlim([0 obj.tsim]);
				end
				axx2 = [axx2 subplot(obj.ipl+1, ppf, obj.ipl*ppf + sp)];
				plot([1:1:ll]'*obj.quantum_us*1e-6,ppwl(idx:idxf,:)'); grid on;
				xlabel(strcat("Time [s] - Cores: ", int2str(idx), " - ", int2str(idxf)));
				xlim([0, obj.tsim]);
				if stop
					break;
				end
			end %for sp
			if stop
				break;
			end
		end %for fi
		
		linkaxes(axx2, 'x');
		
	end %if show

end
		