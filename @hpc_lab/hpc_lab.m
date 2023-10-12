classdef hpc_lab
	%HPC_LAB Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		tm;
	end
	
	%% Gerenal Methods
	methods
		function obj = hpc_lab()
			%HPC_LAB Construct an instance of this class
			%   Detailed explanation goes here
			obj.tm = thermal_model(9,3,3);
		end
		
	end

	%% Simulations
	methods
		[] = sim_tm_autonomous(obj, ts, Nsample, exp_gamma)
	end

	%% Graphs
	methods

	end
end

