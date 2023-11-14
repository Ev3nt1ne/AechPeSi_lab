classdef hpc_lab < thermal_model & power_model & perf_model
	%HPC_LAB Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		%tm;
		tasim = 2; % Free Simulation Time in [s]
		tsim = 1;					% Simulation Time in [s]
		t_init;
		t_outside = 25+273;
		
		Ts_target (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 1e-3;			% Commands min Ts, old Ts_input

		%measure_noise = 1;
		sensor_noise (1,1) {mustBeNonnegative, mustBeNumericOrLogical} ...
			= 1;
		%T_noise_max = 1.0;
		sensor_noise_amplitude (3,1) {mustBeNonnegative, mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1.0; 1.0; 1.0];	% Process, Voltage, Temperature amplitudes.

		graph_show = 1;
		graph_save = 0;

		core_pm;

		x_init;						% Initial Conditions 
		%urplot;						% Input Reference Plot
		frtrc;						% Freq Reference Trace
		%zrplot;						% Power Noise Plot
		wltrc;						% Wokrload Trace

		core_crit_temp = 358.15;	

		%% All Controllers
		tot_pw_budget;
		quad_pw_budget;
		min_pw_red;		


	end

	properties(SetAccess=protected, GetAccess=public)	

		% Internal Variables for Simulation() fnc
	 	wl_index;
	 	qt_storage;

		F_cng_times;
		V_cng_times;
		F_cng_us;
		V_cng_us;
		F_cng_error;
		V_cng_error;

		V_T;
		F_T;

	 	V_s;
	 	F_s;

		delay_F_index;
		delay_V_index;

		delay_F_div;
		delay_V_div;

	 	A_s;
	 	B_s;
	end


	properties(SetAccess=immutable, GetAccess=public)
		customColormap;
		
		PVT_P = 1;
		PVT_V = 2;
		PVT_T = 3;
	end
	
	%% Gerenal Methods
	methods
		function obj = hpc_lab()
			%HPC_LAB Construct an instance of this class
			%   Detailed explanation goes here
			%obj.tm = thermal_model(9,3,3);

			% Object Initialization %%
			% Call superclass constructor before accessing object
			% You cannot conditionalize this statement
			%obj = obj@BaseClass1(args{:});

			% Post Initialization %%
			% Any code, including access to object
			load('-mat', "CustomColorMap.mat");
			obj.customColormap = CustomColormap;
		end
		
	end

	%% Simulations
	methods
		[] = sim_tm_autonomous(obj, ts, Nsample, exp_gamma)
		[cpxplot, cpuplot, cpfplot, cpvplot, wlop] = simulation(obj, ctrl, ctrl_fcn)
		function [obj] = init_compute_model(obj, A, B)
			obj.qt_storage = obj.quantum_instr*ones(obj.Nc, 1);
			obj.wl_index = ones(obj.Nc, 1);
			obj.F_cng_times = zeros(obj.Nc,1);
			obj.V_cng_times = zeros(obj.vd,1);
			obj.F_cng_us = zeros(obj.Nc,1);
			obj.V_cng_us = zeros(obj.vd,1);
			obj.F_cng_error = zeros(obj.Nc,1);
			obj.V_cng_error = zeros(obj.vd,1);

			obj.F_T = obj.F_min*ones(obj.Nc,1);
			obj.V_T = obj.V_min*ones(obj.vd,1);
			obj.F_s = obj.F_min*ones(obj.Nc,1);
			obj.V_s = obj.V_min*ones(obj.vd,1);
			obj.delay_F_index = zeros(obj.Nc,1);
			obj.delay_V_index = zeros(obj.vd,1);
			obj.delay_F_div = ceil(obj.delay_F_mean / obj.Ts);
			obj.delay_V_div = ceil(obj.delay_V_mean / obj.Ts);
			obj.A_s = A;
			obj.B_s = B;
			% cannot put it here, if I want it to be constant among several
			% runs
			%obj.core_pw_noise_char = ones(obj.Nc,1);
			%if (obj.model_variation)
			%obj = obj.create_core_pw_noise();
			%obj = obj.create_thermal_model_noise();
			%end
		end
		function [uplot, xplot, d_is, pw_ms, obj] = compute_model(obj, N, x, V, F, d_p)
			
			d_is = zeros(1,obj.ipl);
			%delay_F_index = zeros(obj.Nc,1);
			%delay_V_index = zeros(obj.vd,1);		
			pw_ms = 0;			
			
			uplot = zeros(N, obj.Nc);
			xplot = zeros(N, obj.Ns);

			%F
			ttd = (obj.delay_F_index > 0)|(obj.F_T ~= obj.F_s);
			cng = (~ttd).*(obj.F_T ~= F);

			obj.F_cng_times = obj.F_cng_times + cng;

			obj.delay_F_index = obj.delay_F_index.*ttd + cng.*obj.delay_F_div;
			obj.F_T = (~cng).*obj.F_T + cng.*F;
			obj.F_cng_error = obj.F_cng_error + ttd.*(obj.F_T ~= F);

			%V
			ttd = (obj.delay_V_index > 0)|(obj.V_T ~= obj.V_s);
			cng = (~ttd).*(obj.V_T ~= V);

			obj.V_cng_times = obj.V_cng_times + cng;

			obj.delay_V_index = obj.delay_V_index.*ttd + cng.*obj.delay_V_div;
			obj.V_T = (~cng).*obj.V_T + cng.*V;
			obj.V_cng_error = obj.V_cng_error + ttd.*(obj.V_T ~= V);
			
			for sim=1:N
				%index = (ix0-1)*N + sim;

				% Delay Application
				tv = obj.delay_V_index == 0;
				obj.V_s = obj.V_s.*(~tv) + obj.V_T.*tv;
				tf = obj.delay_F_index == 0;
				obj.F_s = obj.F_s.*(~tf) + obj.F_T.*tf;

				obj.F_cng_us = obj.F_cng_us + ~tf;
				obj.V_cng_us = obj.V_cng_us + ~tv;

				obj.delay_F_index = obj.delay_F_index - (obj.F_T~=obj.F_s);
				obj.delay_V_index = obj.delay_V_index - (obj.V_T~=obj.V_s);
				%

				%noise 

				instr = obj.F_s * 1e9 * (obj.Ts);
				wl = zeros(obj.Nc, obj.ipl);
				wld = zeros(obj.Nc, 1);
				while (sum(instr) > 0)
					vidx = max(mod(obj.wl_index, size(obj.wltrc,3)),1);
					[m,n,l] = size(obj.wltrc);
					idx = sub2ind([m,n,l],repelem(1:m,1,n),repmat(1:n,1,m),repelem(vidx(:).',1,n));
					wlp = reshape(obj.wltrc(idx),[],m).';
					
					[pwl, instr] = obj.quanta2wl(wlp, instr, (obj.F_min./obj.F_s));					

					% add to accumulator
					wld = wld + pwl;
					wl = wl + pwl .* wlp;

					% compute new values
					obj.qt_storage = obj.qt_storage - pwl;
					ttwl = (obj.qt_storage<=0);
					obj.qt_storage = obj.qt_storage + obj.quantum_instr .* ttwl;
					% obj.quantum_instr*(1-ttwl) + ttwl.*(obj.qt_storage - cinstr);
					obj.wl_index = obj.wl_index + ttwl;
				end

				wl = wl ./ wld;

				d_is = d_is + wl;

				pu_s = obj.power_compute(obj.F_s,obj.VDom*obj.V_s,obj.C(1:obj.Nc,:)*x,wl,d_p);
				pw_ms = pw_ms + sum(pu_s);
				%x = A_nom*x + B_nom*[u;temp_amb*1000];
				x = obj.A_s*x + obj.B_s*[pu_s;obj.temp_amb*1000];		

				uplot(sim,:) = pu_s;
				xplot(sim,:) = x;
			end	%ssim
			d_is = d_is / N;
			pw_ms = pw_ms / N;
		end
	end

	%% Graphs
	methods(Static)
		function [] = savetofile(fig, path_name, bmpres)
			
			% Mange screens
			graph_dpi = 416;
			file_ext = ".png";
			%bmpres = 1; %1.5;
						
			%save as fig
			saveas(fig, strcat(path_name,".fig"));
			%change res for bmp:
			fig.Visible = 'off';
			fig.Position = [1, 1, 1920*bmpres,1080*bmpres];
			exportgraphics(fig, strcat(path_name,file_ext), 'Resolution', graph_dpi);
			close(fig);			
		end
	end
	methods
		function [fig] = xutplot(obj,x,u)
			t = [0:1:length(x)-1]*obj.Ts;
			fig = figure('Name', 'T-P Evolution');
			movegui(fig, 'northwest');
			sbp = 1 + obj.full_model_layers + (obj.add_states>0);
			%TODO: make it parametric. Also make the names of add_state
			%		parametric in thermal_model
			ax1=subplot(sbp,1,1);
			plot(t, x(:,1+end-obj.add_states:end)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Others [°C]'), hold on;
			legend("Heat-sink", "PCB", "Motherboard", "Air");
			%
			ax2=subplot(sbp,1,2);
			plot(t, x(:,2:2:end-obj.add_states)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Heat Spreader (Cu) [°C]'),hold on;
			%
			ax3=subplot(sbp,1,3);
			plot(t, x(:,1:2:end-obj.add_states)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Cores (Si) [°C]'), hold on;
			if max(max(x(:,1:2:end-1))) >= (obj.core_crit_temp-(0.1*(obj.core_crit_temp-273)))
				yline(obj.core_crit_temp-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			end
			%
			ax4=subplot(sbp,1,sbp);
			plot(t, u);
			grid on,xlabel('Time [s]'),ylabel('Cores Power Output [W]'), hold on;
			if max(max(u)) >= (obj.core_max_power-(0.34*obj.core_max_power))
				yline(obj.core_max_power, '--', 'Max Power', 'LineWidth',1,  'Color', 'b');
			end
			if min(min(u)) <= (obj.core_min_power+(0.34*obj.core_max_power))
				yline(obj.core_min_power, '--','min Power', 'LineWidth',1,  'Color', 'b');
			end
			
			linkaxes([ax1,ax2,ax3,ax4],'x');
		end
	end

	%% Libraries
	methods(Static)
		[p,n] = numSubplots(n)
	end

	%% Dependent Variables
	methods

	end
end

