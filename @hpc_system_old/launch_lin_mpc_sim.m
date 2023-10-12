function [xplot, uplot] = launch_lin_mpc_sim(obj, robust, obs, comp)
	if (nargin < 2) || isempty(robust)
		robust = 0;
	end
	if (nargin < 3) || isempty(obs)
		if obj.observable == 1
			obs = 1;
		else
			obs = 0;
		end
	end
	if (nargin < 4) || isempty(comp)
		comp = 0;
	end

	obj = obj.lin_mpc_setup();

	sim_mul = ceil(obj.Ts_mpc/obj.Ts);
	Nsim = ceil(obj.tsim / obj.Ts_mpc);

	failed = 0;

	uplot = zeros(Nsim*sim_mul+1, obj.Ni_c);
	xplot = zeros(Nsim*sim_mul+1, obj.Ns);	

	x = obj.x_init + (rand(obj.Ns,1) - 0.5*ones(obj.Ns,1));
	div_comms = ceil(obj.Ts_input / obj.Ts_mpc);
	uindex = 1;
	u_ref = obj.urplot(uindex,:)';
	%uprev = obj.urplot(1);
	xplot(1,:) = x;
	uplot(1,:) = NaN;

	Adl_true = obj.Ad_true;
	Bdl_true = obj.Bd_true;

	if obs == 1

		Adl_obs = obj.Ad_obs;
		Bdl_obs = obj.Bd_obs;
		poles = ones(obj.Ns,1);
		poles(1:2:obj.Ns) = obj.Obs_poles(1);
		poles(2:2:obj.Ns-1) = obj.Obs_poles(2);

		LK = place(Adl_obs', obj.C', poles)';

		obs_mul = round(obj.Ts_mpc/obj.Ts_obs);
		obs_div = sim_mul / obs_mul;

		xlplot = zeros(Nsim*sim_mul+1, obj.Ns);

		xl = obj.x_init;
		xlplot(1,:) = xl;
		
		C_obs = eye(obj.Ns);
		C_obs(2:2:end,:) = 0;
		C_obs(end,:) = 0;
	end

	for s=1:Nsim
		if (mod(s, div_comms) == 0)
			uindex = uindex + 1;
			u_ref = obj.urplot(uindex,:)';
		end
		if obs == 1
			u = obj.lin_mpc(C_obs*x + (eye(obj.Ns)-C_obs)*xl, obj.temp_amb*1000, u_ref, obj.usum);
		else
			u = obj.lin_mpc(x, obj.temp_amb*1000, u_ref, obj.usum);
		end
		if isnan(u)
			u = obj.core_min_power*ones(obj.Ni_c,1);
			failed = failed + 1;
		%else
		%	uprev = u;
		end

		for ssim=1:sim_mul
			index = (s-1)*sim_mul + ssim;
			u_d = u+obj.zrplot(index,:)';

			%saturate:
			u_d = u_d + (u_d < obj.core_min_power).*(-u_d+obj.core_min_power);
			u_d = u_d + (u_d > obj.core_Max_power*1.2).*(-u_d+obj.core_Max_power*1.2);

			%x = A_nom*x + B_nom*[u;temp_amb*1000];
			x = Adl_true*x + Bdl_true*[u_d;obj.temp_amb*1000];

			if obs == 1
				if mod((ssim-1),obs_div) == 0
					xl = Adl_obs*xl + Bdl_obs*[u;obj.temp_amb*1000] + LK*(obj.C*x - obj.C*xl);
				end
				xlplot(index+1,:) = xl;
			end

			uplot(index+1,:) = u_d;
			xplot(index+1,:) = x;
		end	%for ssim		
	end %for s

	% PLOTs
	t1 = obj.Ts*[0:Nsim*sim_mul]';
	t2 = obj.Ts_input*[0:ceil(obj.tsim/obj.Ts_input)]';

	obj.xutplot(xplot,uplot);
	obj.powerplot(t1,uplot);
	obj.powerconstrplot(uplot);
	obj.tempconstrplot(xplot);
	obj.perfuplot(uplot, t2);			
	if obs == 1
		obj.obsplot(xlplot, xplot);
	end

	% DATA
	disp(strcat('[MPC] number of times the optimization algorithm failed: ',int2str(failed)));
	%failed

	% COMPARISON
	if comp == 1
		x = xplot(1,:)'; % = x_init, but with consistent initial state noise

		soauplot = zeros(Nsim*sim_mul+1,obj.Ni_c);
		soaxplot = zeros(Nsim*sim_mul+1,obj.Ns);

		soaxplot(1,:) = x;
		soauplot(1,:) = NaN;

		Integ = zeros(obj.Nc, 1);

		sim_mul2 = ceil(obj.Ts_ctrl/obj.Ts);
		%TODO: ceil? floor?
		sim_div2 = obj.Ts_mpc/obj.Ts_ctrl; %sim_mul / sim_mul2
		div_comms2 = ceil(obj.Ts_input / obj.Ts_ctrl);
		uindex = 1;
		u_ref = obj.urplot(uindex,:)';

		pid_target = obj.core_crit_temp*ones(obj.Nc,1);
		%if robust == 1
		pid_target = pid_target - obj.Cty(end,1:end-1)';

		for s=1:(Nsim*sim_div2)

			if (mod(s, div_comms2) == 0)
				uindex = uindex + 1;
				u_ref = obj.urplot(uindex,:)';
			end
			%u_ref = obj.urplot(floor(s/sim_div2)+1,:)';

			% PID
			[upid, Integ] = obj.cp_pid((obj.C(1:end-1,:)*x),pid_target,Integ,u_ref);

			% Controller Output
			u = u_ref + upid;

			for ssim=1:sim_mul2

				index = (s-1)*sim_mul2 + ssim;

				%noise 
				u_d = u+obj.zrplot(index,:)';

				%x = A_nom*x + B_nom*[u;temp_amb*1000];
				x = Adl_true*x + Bdl_true*[u_d;obj.temp_amb*1000];		

				soauplot(index+1,:) = u_d;
				soaxplot(index+1,:) = x;
			end	%ssim
		end %s	

		% PLOTs
		t1 = obj.Ts*[0:Nsim*sim_mul]';
		t2 = obj.Ts_input*[0:ceil(obj.tsim/obj.Ts_input)]';

		obj.xutplot(soaxplot,soauplot);
		obj.powerplot(t1,uplot);
		obj.powerconstrplot(uplot, soauplot);
		obj.tempconstrplot(xplot,soaxplot);
		obj.perfuplotcomp(uplot, soauplot, t2);

	end %comp==1

end	%launch_lin_mpc_sim		
		