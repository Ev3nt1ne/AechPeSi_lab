function [cpxplot, cpuplot, cpfplot, cpvplot, wlop] = ...
	simulation(obj, ctrl, show)
%SIMULATION Summary of this function goes here
%   Detailed explanation goes here

	% To Optimize execution, istead of having several functions called and
	%	several ifs statement, we use ctrl.Ts_ctrl as the main time step of
	%	the simulation.
	Nsim = ceil(obj.tsim / ctrl.Ts_ctrl);

	% System:
	% initial offset:
	sys_ts_offset = 0;
	% frequency
	sys_mul = round(ctrl.Ts_ctrl / obj.Ts);
	% output completed
	%ctrl_output_mul = round(obj.ctrl_Texec / obj.Ts);

	% Inputs
	% initial offset:
	target_ts_offset = 0;
	% frequency
	target_mul = round(obj.Ts_target / ctrl.Ts_ctrl);
	target_counter = target_mul + target_ts_offset;
	target_index = 0;

	% INIT
	obj = obj.init_compute_model(obj.Ad_true, obj.Bd_true);

	%TODO
	x = obj.x_init + (rand(obj.Ns,1) - 0.5*ones(obj.Ns,1));
	F = obj.F_min*ones(obj.Nc,1);
	V = obj.V_min*ones(obj.vd,1);
	process = ones(obj.Nc,1);
	cpuplot = zeros(Nsim+1,obj.Nc);
	cpxplot = zeros(Nsim+1,obj.Ns);
	cpxplot(1,:) = x;
	cpuplot(1,:) = NaN;
	cpfplot = zeros(Nsim+1,obj.Nc);
	cpvplot = zeros(Nsim+1,obj.vd);
	cpfplot(1,:) = F;
	cpvplot(1,:) = V;
	pwm = 1;
	wl  = 1;

	ctrl.init_fnc(obj);

	% LOOOP
	for s=1:Nsim

		% Input step managing
		target_counter = (target_counter-1)*(target_counter>0) + (target_counter<=0)*(target_mul-1);
		target_index = target_index + (target_counter==target_mul-1);

		ctrl_pwm = pwm;
		ctrl_wl = wl;
		
		%Compute model:
		index = 1+(s-1)*sys_mul;
		[cpuplot(index+1:index+sys_mul,:), cpxplot(index+1:index+sys_mul,:), wl, pwm, obj] = obj.compute_model(sys_mul, cpxplot(index,:)', V, F, process);	


		% Sim output managing
		T = obj.C*cpxplot(index,:)';
		% Noise:
		T = T + (obj.sensor_noise)*( (rand(size(T)) - 0.5)*2 * obj.sensor_noise_amplitude(obj.PVT_T) );	

		%if (mod(s-1 + ctrl_ts_offset, ctrl_mul) == 0)
			pvt = {process, [], T};
			[F, V, ctrl] = ctrl.ctrl_fnc(obj, target_index, pvt, ctrl_pwm, ctrl_wl);
		%end

		cpfplot(s+1,:) = F;
		cpvplot(s+1,:) = V;
	end

	wlop = obj.wl_index / (size(obj.wrplot,3)-1) * 100;
	ctrl.cleanup_fnc(obj);

	if show
		obj.xutplot(cpxplot,cpuplot);
		obj.powerconstrplot(cpuplot);
		obj.tempconstrplot(cpxplot);
		obj.perfplot(cpfplot,obj.wl_index);		
		obj.fvplot(cpfplot,cpvplot);

		ctrl.plot_fnc(obj);
	end

	
end

