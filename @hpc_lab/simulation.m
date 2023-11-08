function [cpxplot, cpuplot, cpfplot, cpvplot, wlop] = ...
	simulation(obj, ctrl, ctrl_fcn)
%SIMULATION Summary of this function goes here
%   Detailed explanation goes here

	Nsim = ceil(obj.tsim / obj.Ts);

	% Controller:
	% initial offset:
	ctrl_ts_offset = 0;
	% frequency
	ctrl_mul = round(ctrl.Ts_ctrl / obj.Ts);
	% output completed
	%ctrl_output_mul = round(obj.ctrl_Texec / obj.Ts);

	% Inputs
	% initial offset:
	target_ts_offset = 0;
	% frequency
	target_mul = round(obj.Ts_target / obj.Ts);
	target_counter = target_mul + target_ts_offset;
	target_index = 0;

	% INIT
	obj = obj.init_compute_model(obj.Ad_true, obj.Bd_true);

	%TODO
	x = obj.x_init + (rand(obj.Ns,1) - 0.5*ones(obj.Ns,1));
	F = obj.F_min*ones(obj.Nc,1);
	V = obj.V_min*ones(obj.vd,1);
	process = ones(obj.Nc,1);
	cpuplot = zeros(Nsim+1,obj.Ni);
	cpxplot = zeros(Nsim+1,obj.Ns);
	cpxplot(1,:) = x;
	cpuplot(1,:) = NaN;
	cpfplot = zeros(Nsim+1,obj.Nc);
	cpvplot = zeros(Nsim+1,obj.vd);
	cpfplot(1,:) = F;
	cpvplot(1,:) = V;

	% LOOOP
	for s=1:Nsim

		% Input step managing
		target_counter = (target_counter-1)*(target_counter>0) + (target_counter<=0)*(target_mul-1);
		target_index = target_index + (target_counter==target_mul-1);

		% Sim output managing
		T = obj.C*cpxplot(s,:)';
		% Noise:
		T = T + (obj.sensor_noise)*( (rand(size(T)) - 0.5)*2 * obj.sensor_noise_amplitude(obj.PVT_T) );	

		%if (mod(s-1 + ctrl_ts_offset, ctrl_mul) == 0)
		%	pvt = {process, [], T};
		%	[F, Vc] = ctrl_fcn(obj, ctrl, target_index, pvt, pws);
		%end

		%V = obj.VDom * Vc;
		
		%Compute model:
		[cpuplot(s+1,:), cpxplot(s+1,:), d_is, pw_ms, obj] = obj.compute_model(1, cpxplot(s,:)', V, F, process);	

	end


end

