classdef power_model < handle
	%POWER_MODEL Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		
		leak_exp (1,1) {mustBeNonnegative, mustBeNumericOrLogical} ...
			= 1;

		leak_vdd_k (1,1) {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 324.29e-3;
		leak_temp_k	(1,1) {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 1.0;
		leak_process_k (1,1) {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 528.04e-3;

		leak_exp_vdd_k (1,1) {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 4507.748e-3;
		leak_exp_t_k (1,1) {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 29.903e-3;
		leak_exp_k (1,1) {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= -6033.035e-3;
		 
		dyn_ceff_k {mustBeNumeric, mustBeNonempty, mustBeFinite, mustBeVector} ...
			= [306.741, 694.866, 1235.452, 1651.217, 1600.92]*1e-3;

		FV_table (:,3) {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [0.50, 0.40, 1.3500;
				0.55, 0.40, 1.6000;
				0.60, 1.35, 1.8000;
				0.65, 1.35, 2.0000;
				0.70, 1.35, 2.2000;
				0.75, 1.60, 2.4000;
				0.80, 1.80, 2.6000;
				0.85, 1.80, 2.7500;
				0.90, 2.00, 2.9000;
				0.95, 2.00, 3.0500;
				1.00, 2.00, 3.2000;
				1.05, 2.60, 3.3500;
				1.10, 2.60, 3.4500;
				1.15, 2.60, 3.5500;
				1.20, 2.60, 3.6600];

		pw_dev_per {mustBeNonnegative, mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 1;
			 
		%
		delay_F_mean (1,1) {mustBeNonnegative, mustBeNumeric, mustBeFinite} ...
			= 1e-5;			% Mean Frequency application Delay
		delay_V_mean (1,1) {mustBeNonnegative, mustBeNumeric, mustBeFinite} ...
			= 1e-5;			% Mean Voltage application Delay
		delay_F_max (1,1) {mustBeNonnegative, mustBeNumeric, mustBeFinite} ...
			= 5e-5;			% Max Frequency application Delay
		delay_V_max (1,1) {mustBeNonnegative, mustBeNumeric, mustBeFinite} ...
			= 2e-5;			% Max Voltage application Delay
		%
		F_discretization_step (1,1) {mustBeNonnegative, mustBeNumeric, mustBeFinite} ...
			= 0.05;	

		% Power Noise
		%TODO 
		pw_gmean = 0.02;
		pw_gvar = 1;
		pw_glim = [-0.02 0.05];
		pw_3sigma_on3 = 0.05;
		
	end

	properties(Dependent)
		% System Structure
		ipl;						% Instruction power levels
		FV_levels;

		% System Parameters
		static_pw_coeff;
		polyFV;
		
		V_max;
		V_min;
		F_max;
		F_min;
		core_max_power;
		core_min_power;	
	end
	
	methods
		function obj = power_model()
			%POWER_MODEL Construct an instance of this class
			%   Detailed explanation goes here
			%obj.Property1 = inputArg1 + inputArg2;
			obj.create_core_pw_noise();
		end
		function obj = anteSimCheckPM(obj)
			%TODO here obj.Nc is wrong
			if length(obj.pw_dev_per) ~= obj.Nc
				disp("recreating");
				% TODO here instead of recreating, add missing ones or
				%	reduce it
				obj.create_core_pw_noise();
			end
		end
	end

	methods
		function ps = ps_compute(obj, V,T,process, noise)

			% F,V,d_i,d_p can be any dim
			% T has to be dim = #cores
			dim = length(T);
			
			kps = [obj.leak_vdd_k, obj.leak_temp_k, obj.leak_process_k, ...
					obj.leak_exp_vdd_k, obj.leak_exp_t_k, obj.leak_exp_k];

			ps = V*kps(1) + process*kps(3);
			%TODO: is this better as an if, or as a "vectorized if"
			%computation?
			if obj.leak_exp
				maxT = 125+273.15;
				ps = ps .* ...
					exp(V*kps(4) + (min(T,ones(dim,1)*maxT)-273.15)*kps(5) + kps(6));
			end
		end
		function pd = pd_compute(obj, F,V,instr,noise)

			%TODO 
			% F,V,d_i,d_p can be any dim
			% T has to be dim = #cores
			dim = length(F);

			tt = dim > 1;
			npw = obj.pw_dev_per(1:dim) * (tt && noise) + ~(tt && noise);
			
			Ceff = instr * obj.dyn_ceff_k' .* npw;
			pd = Ceff .* F .* (V .* V);
		end
		function pu = power_compute(obj, F, V, T, instr, process, noise)

			if nargin < 7 || isempty(noise)
				noise = 0;
			end

			Power_static = obj.ps_compute(V,T,process, noise);
			Power_dyn = obj.pd_compute(F,V,instr,noise);			

			pu = Power_static + Power_dyn;		
		end
		%TODO: here I use obj.Nc!!! FIX!
		function obj = create_core_pw_noise(obj)
			Covr = chol(obj.pw_gvar);
			obj.pw_dev_per = ones(obj.Nc, 1);
			for c=1:obj.Nc
				stop = 0;
				while(~stop)	
					z = repmat(obj.pw_gmean,1,1) + randn(1,1)*Covr * obj.pw_3sigma_on3/3;	
					if (z>obj.pw_glim(1)) && (z<obj.pw_glim(2))
						stop = 1;
					end
				end %while
				obj.pw_dev_per(c) = obj.pw_dev_per(c) + z;
			end %for
		end %func
	end

	%% Dependent Variables
	methods
		function value = get.static_pw_coeff(obj)
			value = [obj.leak_vdd_k, obj.leak_temp_k, obj.leak_process_k];
		end
		function value = get.FV_levels(obj)
			value = size(obj.FV_table,1);
		end
		function value = get.V_max(obj)
			value = obj.FV_table(end,1);
		end
		function value = get.V_min(obj)
			value = obj.FV_table(1,1);
		end
		function value = get.F_min(obj)
			value = obj.FV_table(1,2);
		end
		function value = get.F_max(obj)
			value = obj.FV_table(end,3);
		end
		function value = get.core_max_power(obj)
			value = 10.5348;
			%TODO: here put the call to the function model
		end
		function value = get.core_min_power(obj)
			value = 0.8939;
			%TODO: here put the call to the function model
		end
		function value = get.ipl(obj)
			value = size(obj.dyn_ceff_k,2);
		end
		function value = get.polyFV(obj)
			value = polyfit(obj.FV_table(:,3),obj.FV_table(:,1), 3);
		end
	end
end

