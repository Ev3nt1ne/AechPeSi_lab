classdef power_model
	%POWER_MODEL Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		exp_leakage = 1;

		leak_vdd = 324.29;
		leak_temp = 1.0;
		leak_process = 528.04;
		
		ceff_pw_coeff = [306.741, 694.866, 1235.452, 1651.217, 1600.92];
		%[466.741, 694.866, 1235.452, 1651.217, 1600.92];
		exp_leak_coeff = [4507.748, 29.903, 6033.035];
		
		min_pw_red;

		FV_table = [0.50, 0.40, 1.3500;
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
			 
		% Power Noise
		%pn_mu = 0.1;
		%pn_sigma = 2.0; %0.5
		%pn_alpha = 0.2; %0.8
		pw_gmean = 0.02;
		pw_gvar = 1;
		pw_glim = [-0.02 0.05];
		pw_3sigma_on3 = 0.05;
		
		core_pw_noise_char;

		%
		delay_F_mean = 1e-5;		% Mean Frequency application Delay
		delay_V_mean = 1e-5;		% Mean Voltage application Delay
		delay_F_max = 5e-5;			% Max Frequency application Delay
		delay_V_max = 2e-5;			% Max Voltage application Delay
		%
		F_discretization_step = 0.05;

		%% WL Parameters
		
		wl_prob = [0.0345 0.3448 0.3448 0.2414 0.0345];			% Probability for each WL	
		wl_min_exec_us = [510 100 100 480 520];					% minimal execution in us
		wl_mean_exec_us = [520e2 210e1 180e1 560e2 1.2e6/1e1];	% mean execution in us

		wl_uniq_min = [0.9 0.3 0.4 0.6 0.9];					% minimal percentage of presence when Primary Wl
		wl_comb_array = [0 1 0 0 0;								% Combination for each wl with each other:
						 1 0 1 1 1;								%	1 = they can be present together, 0 = they cannot
						 0 1 0 1 1;								%	Matrix has to be symmetric, with 0 on the diagonal.
						 0 1 1 0 1;
						 0 1 1 1 0];
		max_num_wl = 3;											% max number of contemporary wls (primary [1] + secondary). It may not be always respected.
		mem_wl = [0.95 0.15 0.05 0.20 0];							% memory boundness
		
	end

	properties(Dependent)
		%% System Structure
		ipl;						% Instruction power levels
		FV_levels;
		quantum_instr;

		%% System Parameters
		static_pw_coeff;
		polyFV;
		
		V_Max;
		V_min;
		F_Max;
		F_min;
		core_Max_power;
		core_min_power;	
	end
	
	methods
		function obj = power_model(inputArg1,inputArg2)
			%POWER_MODEL Construct an instance of this class
			%   Detailed explanation goes here
			%obj.Property1 = inputArg1 + inputArg2;
		end
		
		function outputArg = method1(obj,inputArg)
			%METHOD1 Summary of this method goes here
			%   Detailed explanation goes here
			%outputArg = obj.Property1 + inputArg;
		end
	end
end

