function [cpxplot, cpuplot, cpfplot, cpvplot, wlop] = ...
	simulation(obj, ctrl, show)
%SIMULATION Summary of this function goes here
%   Detailed explanation goes here	

	obj.anteSimCheckTM();
	obj.anteSimCheckLab();
	obj.anteSimCheckPM();

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
	ts_div = obj.Ts_target / ctrl.Ts_ctrl;
	target_mul = floor(ts_div);
	ctrl_mul = floor(1/ts_div) -1;
	ctrl_mul = (ctrl_mul<0)*0 + (ctrl_mul>=0)*ctrl_mul;
	target_counter = target_mul + target_ts_offset;
	target_index = 0;
	% Manage non-integer part:
	% TODO: This does not work well with y.xx when xx> 60 and y = [2,3,4,5] 
	%	because it removes too much. instead with 1, and 6, 7... works well 
	%	I should investigate the math further.
	if target_mul>ctrl_mul
		ntd = (ts_div-target_mul)/target_mul/target_mul;
		nip_sign = -1;
	else
		ntd = (1/ts_div)-(ctrl_mul+1);
		nip_sign = +1;
	end
	nip_mul = ceil(1/ntd);
	nip_counter = nip_mul + target_ts_offset;

	% INIT
	obj.init_compute_model(obj.Ad_true, obj.Bd_true);

	%TODO
	x = obj.t_init + (rand(obj.Ns,1) - 0.5*ones(obj.Ns,1));
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
	wl  = [ones(obj.Nc,1) zeros(obj.Nc, obj.ipl -1)];
	pvt{obj.PVT_P} = process;
	pvt{obj.PVT_V} = [];
	pvt{obj.PVT_T} = ctrl.C*x;

	ctrl = ctrl.init_fnc(obj, Nsim);

	% LOOOP
	for s=1:Nsim

		% Input step managing
		target_counter = (target_counter-1)*(target_counter>0) + (target_counter<=0)*(target_mul-1);
		nip_counter = (nip_counter-1)*(nip_counter>0) + (nip_counter<=0)*(nip_mul-1);
		target_index = target_index + (target_counter==target_mul-1) + nip_sign*(nip_counter==nip_mul-1) + ctrl_mul;

		ctrl_pwm = pwm;
		ctrl_wl = wl;
		
		%Compute model:
		index = 1+(s-1)*sys_mul;
		[cpuplot(index+1:index+sys_mul,:), cpxplot(index+1:index+sys_mul,:), wl, pwm, obj] = obj.compute_model(sys_mul, cpxplot(index,:)', V, F, process);	

		% Sim output managing
		T = ctrl.C*cpxplot(index,:)';
		% Noise:
		dim = min(length(T), obj.Nc);
		T = T + [(obj.sensor_noise)*( (rand(dim,1) - 0.5)*2 * obj.sensor_noise_amplitude(obj.PVT_T) ); zeros(length(T)-dim,1)];	

		%if (mod(s-1 + ctrl_ts_offset, ctrl_mul) == 0)
			pvt{obj.PVT_P} = process;
			pvt{obj.PVT_V} = [];
			pvt{obj.PVT_T} = T;
			[F, V, ctrl] = ctrl.ctrl_fnc(obj, target_index, pvt, ctrl_pwm, ctrl_wl);
		%end

		cpfplot(s+1,:) = F;
		cpvplot(s+1,:) = V;
	end

	wlop = obj.wl_index / (size(obj.wltrc,3)-1) * 100;
	ctrl = ctrl.cleanup_fnc(obj);

	if show
		% Pause because it is bugged on Linux
		pause(0.5);
		obj.xutplot(cpxplot,cpuplot);
		pause(0.5);
		obj.powerconstrplot(cpuplot);
		pause(0.5);
		obj.tempconstrplot(cpxplot);
		pause(0.5);
		obj.perfplot(cpfplot,obj.wl_index);	
		pause(0.5);
		obj.fvplot(cpfplot,cpvplot);

		t1 = obj.Ts*[1:Nsim*sys_mul]';
		t2 = obj.Ts*sys_mul*[1:Nsim]';

		ctrl = ctrl.plot_fnc(obj, t1, t2, cpxplot, cpuplot, cpfplot, cpvplot, wlop);
	end

end

