classdef perf_model
	%PERF_MODEL Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
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
	
	methods
		function obj = perf_model(inputArg1,inputArg2)
			%PERF_MODEL Construct an instance of this class
			%   Detailed explanation goes here
			
		end
		
	end
end

