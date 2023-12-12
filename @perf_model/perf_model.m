classdef perf_model < handle
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
		dur_c = 5/6; %TODO: what is this?
		min_secondary_space = 0.075; % minimum percentage of wl to be considered as secondary. Otherwise discarded	
	end
	
	methods
		function obj = perf_model(inputArg1,inputArg2)
			%PERF_MODEL Construct an instance of this class
			%   Detailed explanation goes here
			
		end
		function [sp] = wl_prob_compute(obj, primary_wl, new_wl, ipl)
			
			hm_min = obj.wl_uniq_min(primary_wl);
			
			% Compute How Much?
			hmp = rand(1) + hm_min*obj.dur_c;
			hmp = hmp + (hmp>1)*(1-hmp) + (hmp<hm_min)*(hm_min-hmp);
			% Compute prob secondary wls
			%sp = (rand(1,5)-obj.wl_prob).*(obj.wl_comb_array(primary_wl,:)|(new_wl>0));
			sp = rand(1,ipl);
			ttsp = (sp>(1-obj.wl_prob));
			if sum(ttsp) > 0
				sp = ttsp .* sp;
			else
				sp = sp .* obj.wl_prob;
			end	
			sp = sp.*(obj.wl_comb_array(primary_wl,:)|(new_wl>0));
			sp = sp - (sp<=0).*sp + (sp<=0.001)*0.001;
			sp = sp + new_wl;
			% Remove excessive secondary wls
			while sum(sp>0) > obj.max_num_wl-1
				spmin = min(sp(sp>0&sp<1));
				if isempty(spmin)
					break;
				end
				sp = sp - sp.*(sp==spmin);
			end
			% Compute how much? secondary wls
			sp = sp - new_wl;
			% Fix case sp=0 and new_wl != 0
			if sum(sp) == 0
				sp = obj.wl_comb_array(primary_wl,:);
			end
			% Normalize on 1-hmp, and forzing hmp to be <1 in case new_wl != 0
			if sum(new_wl)>0
				max_hmp = max(hm_min, 1 - obj.min_secondary_space*min(obj.max_num_wl,sum(obj.wl_comb_array(primary_wl,:)))); 
				hmp = hmp - (hmp>max_hmp)*(hmp-(max_hmp));
			end
			sp = sp/sum(sp)*(1-hmp);
			
			sp(primary_wl) = sp(primary_wl) + hmp;
			
		end
		function wl_trace = generate_wl_trace(obj, Nc, ts, show)
			
			if (nargin < 2) || isempty(Nc)
				warning("[PM Error]: Nc missing or empty!");
				disp(strcat("[PM] Insert number of Cores (Nc)."));
				wl_trace = [];
				return;
			end
			if (nargin < 3) || isempty(ts)
				warning("[PM Error]: ts missing or empty!");
				disp(strcat("[PM] Insert the simulation time (ts) in seconds."));
				wl_trace = [];
				return;
			end
			if (nargin < 4) || isempty(show)
				show = 0;
			end

			% Init
			%TODO: perform a check somewhere
				% obj.quantum_us = max(obj.Ts*1e6, obj.quantum_us);
				% if (obj.Ts*1e6 > obj.quantum_us)
				%	disp error, and obj.quantum_us = obj.Ts*1e6;
			ll = ceil((ts*1e6)/obj.quantum_us);
			ipl = length(obj.wl_prob);
			wl_trace = zeros(Nc,ipl, ll);

			wl_min_exec_qt = ceil(obj.wl_min_exec_us / obj.quantum_us);
			wl_mean_exec_qt = ceil(obj.wl_mean_exec_us / obj.quantum_us) - wl_min_exec_qt;
			wl_mean_exec_qt(wl_mean_exec_qt<=0) = 0;

			%check on normalization
			obj.wl_prob = obj.wl_prob/sum(obj.wl_prob);
		
			% create distribution graph
			cov = wl_mean_exec_qt/3; %sqrt(input)*1.8;
			cov(cov<=0) = 1;
			x={};
			y={};
			%figure()
			for i=1:ipl
				x{i}(:) = 1:1:wl_mean_exec_qt(i)*2;
				yn = exp(-(x{i}-wl_mean_exec_qt(i)).^2 ./ (2*cov(i).^2));
				y{i}(:) = yn ./ sum(yn);
				%plot(x{i},y{i}), hold on
			end
		
			% ppwl: holding information on primary wl for showing graphs
			if show
				ppwl = zeros(Nc, ll);
			end
			% I tried, but it's hard to make it vectorial for each core
			for c=1:Nc
			%init
			new_wl = zeros(1,ipl);
			primary_wl = 1;
			wi = 1;
			while wi<=ll
		
				% Add min execution requirements		
				for nl=1:wl_min_exec_qt(primary_wl)		
					%call function for wls
					sp = obj.wl_prob_compute(primary_wl, new_wl, ipl);
					
					% Update new_wl
					new_wl = new_wl + ((sp>0) .* (1-(new_wl>0))).*wl_min_exec_qt;
					new_wl = new_wl - (new_wl>0);
					new_wl = new_wl + (new_wl<0).*(-new_wl);
					
					% Add
					wl_trace(c,:,wi) = sp;
		
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
					sp = obj.wl_prob_compute(primary_wl, new_wl, ipl);
					
					% Update new_wl
					new_wl = new_wl + ((sp>0) .* (1-(new_wl>0))).*wl_min_exec_qt;
					new_wl = new_wl - (new_wl>0);	
					new_wl = new_wl + (new_wl<0).*(-new_wl);
					
					% Add
					wl_trace(c,:,wi) = sp;
		
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
					if counter >= wl_mean_exec_qt(primary_wl)*2
						p = -1;
					else
						p = 1 - y{primary_wl}(counter)/(1-sum(y{primary_wl}(1:counter-1)));
					end					
					stop = cv > p;
					
					%clc;
					
				end %while stop==0
				if (wi > ll)
					break;
				end
		
				% Choose new primary workload
				% do not clear new_wl
				sp = rand(1,ipl);
				ttsp = (sp>(1-obj.wl_prob));
				if sum(ttsp) > 0
					sp = ttsp .* sp;
				else
					sp = sp .* obj.wl_prob;
				end
				primary_wl = sum((sp==max(sp)) .* [1:1:ipl]);
		
			end %external while
			end %for core


			% Show PLOTS
			if show		
				figure('Name', 'Show Workload Reference');
				axx = [];
				for i=1:ipl
					hold on;
					axx = [axx subplot(ipl+1, 1, i)];
					plot([1:1:ll]'*obj.quantum_us*1e-6,squeeze(wl_trace(:,i,:)));
				end
				axx = [axx subplot(ipl+1, 1, ipl+1)];
				plot([1:1:ll]'*obj.quantum_us*1e-6,ppwl'); grid on;
				xlabel("Time [s]");
				linkaxes(axx, 'x');
				
				npp = 4; %number per plot
				ppf = 3; %plot per figure
				nf = ceil(Nc/(npp*ppf)); %number of figures
				axx2 = [];
				stop = 0;
				for fi=1:nf
					figure('Name', strcat("Wl Reference core: ", int2str((fi-1)*(npp*ppf)+1), " - ", int2str((fi-1)*(npp*ppf)+(npp*ppf)) ) );
					hold on;
					for sp=1:ppf
						idx = (fi-1)*(npp*ppf) + (sp-1)*npp +1;
						idxf = idx+npp-1;
						if idxf >= Nc
							stop = 1;
							idxf = Nc;
						end
						for pp=1:ipl
							hold on;
							axx2 = [axx2 subplot(ipl+1, ppf, (pp-1)*ppf + sp)];
							plot([1:1:ll]'*obj.quantum_us*1e-6,squeeze(wl_trace(idx:idxf,pp,:)));
							xlim([0 ll*obj.quantum_us*1e-6]);
						end
						axx2 = [axx2 subplot(ipl+1, ppf, ipl*ppf + sp)];
						plot([1:1:ll]'*obj.quantum_us*1e-6,ppwl(idx:idxf,:)'); grid on;
						xlabel(strcat("Time [s] - Cores: ", int2str(idx), " - ", int2str(idxf)));
						xlim([0, ll*obj.quantum_us*1e-6]);
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
		function [pwl, res] = quanta2wl(obj, wlp, instr, mem_instr_conv)
			mem_wgt = sum(wlp.*obj.wl_mem_weigth,2).*(instr>0);
			cinstr = instr.*(mem_wgt + (1-mem_wgt).*mem_instr_conv); %mem_wgt.*mem_instr_conv.*instr + (1-mem_wgt).*instr;
						
			pwl = min(cinstr, obj.qt_storage); %obj.qt_storage./cinstr;
			res = round(instr - pwl.*(instr./cinstr));
			%((pwl - mem_wgt.*mem_instr_conv) ./ (1-mem_wgt));
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

