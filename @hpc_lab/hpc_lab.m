classdef hpc_lab < thermal_model
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
		urplot;						% Input Reference Plot
		frplot;						% Freq Reference Plot
		zrplot;						% Power Noise Plot
		wrplot;						% Wokrload Plot

		core_crit_temp = 358.15;	

		%% All Controllers
		tot_pw_budget;
		quad_pw_budget;
		min_pw_red;		


		% TO CANC
		%F_min = 0.8;
		%ipl = 5;
		core_max_power = 15;
		core_min_power = 0.4;


	end

	properties(Dependent)
		quantum_instr;
	end

	properties(SetAccess=protected, GetAccess=public)		
		% model
	 	wl_index;
	 	qt_storage;
	 	V_s;
	 	F_s;
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
			obj.F_s = obj.F_min*ones(obj.Nc,1);
			obj.V_s = obj.V_min*ones(obj.vd,1);
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
			delay_F_index = zeros(obj.Nc,1);
			delay_V_index = zeros(obj.vd,1);		
			pw_ms = 0;
			
			delay_F_div = ceil(obj.delay_F_mean / obj.Ts);
			delay_V_div = ceil(obj.delay_V_mean / obj.Ts);
			
			uplot = zeros(N, obj.Ni_c);
			xplot = zeros(N, obj.Ns);
			
			for sim=1:N
				%index = (ix0-1)*N + sim;

				% Delay Application
				% Application of V before F
				delay_F_index = delay_F_index + obj.VDom*(V==obj.V_s);
				delay_V_index = delay_V_index + (V~=obj.V_s);
				tv = delay_V_index > (delay_V_div*ones(obj.vd,1));
				obj.V_s = obj.V_s - obj.V_s.*tv + V.*tv;
				%
				tf = delay_F_index > (delay_F_div*ones(obj.Nc,1));
				obj.F_s = obj.F_s - obj.F_s.*tf + F.*tf;
				%noise 

				%wl
				wlp = zeros(obj.Nc, obj.ipl);
				wlpn = zeros(obj.Nc, obj.ipl);
				for c=1:obj.Nc
					wlp(c,:) = obj.wrplot(c,:,max(mod(obj.wl_index(c), size(obj.wrplot,3)),1));
					wlpn(c,:) = obj.wrplot(c,:,max(mod(obj.wl_index(c)+1, size(obj.wrplot,3)),1+1));
				end
				if (obj.F_s(1) > obj.F_min)
					asd = 1;
				end
				instr = obj.F_s * 1e9 * (obj.Ts);
				mem_instr = obj.F_min * 1e9 * (obj.Ts);
				a_citr = sum(wlp.*obj.mem_wl,2);
				a_citrn = sum(wlpn.*obj.mem_wl,2);
				citr = a_citr.*mem_instr + (1-a_citr).*instr;
				citrn = a_citrn.*mem_instr + (1-a_citrn).*instr;
				pwl = obj.qt_storage./citr;
				ttwl = (pwl>1);
				pwl = pwl + ttwl.*(1-pwl);

				wl = pwl .* wlp + ...
					(1-pwl) .* wlpn;

				obj.qt_storage = obj.quantum_instr*(1-ttwl) - (citr - obj.qt_storage);
				obj.wl_index = obj.wl_index + (1-ttwl);

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
			value = obj.F_max * 1e9 * (obj.quantum_us * 1e-6);
		end
	end
end

