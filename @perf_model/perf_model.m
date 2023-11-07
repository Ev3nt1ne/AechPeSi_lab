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
		wl_mem_weigth = [0.95 0.15 0.05 0.20 0];				% memory boundness

		quantum_us = 50;
		quantum_F = 1; %GHz
	end

	properties(Dependent)
		quantum_instr;
	end

	properties(SetAccess=protected, GetAccess=public)		
		% wl:
		% How much? coefficient
		dur_c = 5/6;
		min_secondary_space = 0.075;		
	end
	
	methods
		function obj = perf_model(inputArg1,inputArg2)
			%PERF_MODEL Construct an instance of this class
			%   Detailed explanation goes here
			
		end
		function [pwl, res] = quanta2wl(wlp, instr, mem_instr_conv)
			mem_wgt = sum(wlp.*obj.wl_mem_weigth,2).*(instr>0);
			cinstr = mem_wgt.*mem_instr_conv + (1-mem_wgt).*instr;
						
			pwl = min(cinstr, obj.qt_storage); %obj.qt_storage./cinstr;
			res = instr - ((pwl - mem_wgt.*mem_instr_conv) ./ (1-mem_wgt));
			% saturate to 0
			% instr = ..
		end

		
	end

	%% Dependent Variables
	methods
		function value = get.quantum_instr(obj)
			% quantum_instr has to be >= than the maximum reachable
			% frequency. This is needed because in "compute_model" we
			% assumed that the maximum number of quanta considered for the
			% computation of wl is 2. In case we compute quantum_instr with
			% a value < F_Max it is possible that in the qt_storage there
			% is only a small fraction of instructions, and so the F_s will
			% take the small fraction, a whole new quanta, plus a little
			% part of a third quanta. If the chosen value here is << F_Max
			% the number could be even greater.
	
			%This has to be fixed, since quanta are considered at the
			%benchmark executed frequency, so not always max one.
			%TODO
			value = obj.quantum_F * 1e9 * (obj.quantum_us * 1e-6);
		end
	end

	
end

