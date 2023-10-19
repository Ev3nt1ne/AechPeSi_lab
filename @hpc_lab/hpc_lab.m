classdef hpc_lab < thermal_model
	%HPC_LAB Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		%tm;
		tasim = 2; % Free Simulation Time in [s]
		tsim = 1;					% Simulation Time in [s]
		t_init;
		t_outside = 25+273;

		graph_show = 1;
		graph_save = 0;

		core_pm;

		x_init;						% Initial Conditions 
		urplot;						% Input Reference Plot
		frplot;						% Freq Reference Plot
		zrplot;						% Power Noise Plot
		wrplot;						% Wokrload Plot


		% TO CANC
		%F_min = 0.8;
		%ipl = 5;
		core_crit_temp = 85+273;
		core_max_power = 15;
		core_min_power = 0.4;
	end

	properties(SetAccess=immutable, GetAccess=public)
		customColormap;
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
end

