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

		osunix;
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

			if isunix
				obj.osunix = 1;
			else
				obj.osunix = 0;
			end
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
					obj.wl_index = obj.wl_index + ttwl;
				end

				wl = wl ./ wld;

				d_is = d_is + wl;

				pu_s = obj.power_compute(obj.F_s,obj.VDom*obj.V_s,obj.C(1:obj.Nc,:)*x,wl,d_p);
				pw_ms = pw_ms + sum(pu_s);
				x = obj.A_s*x + obj.B_s*[pu_s;obj.temp_amb*1000];		

				uplot(sim,:) = pu_s;
				xplot(sim,:) = x;
			end	%ssim
			d_is = d_is / N;
			pw_ms = pw_ms / N;
		end
	end

	%% Approximation
	methods
		[k0, k1, k2] = pws_ls_approx(obj)
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
		function [fig] = powerconstrplot(obj,u1,u2)
						
			fig = figure('Name', 'Power Constraints Compliance');
			movegui(fig, 'southeast');
			title('Power Constraints analysis')
			
			t1 = [0:1:length(u1)-1]*obj.Ts;
			
			ax1 = subplot(3,4,[1:4]);
			gr = sum(u1,2);
			plot(t1,gr, '.');
			grid on,title('Total Power'),ylabel('[W]'),xlabel('Time [s]'),hold on;
			if ((nargin >=4) && (isempty(u2)==0))
				gr = sum(u2,2);
				plot(t1,gr);
			end
			%if max(gr) >= obj.usum(1)*0.2
			plot([0:max(length(obj.tot_pw_budget)-1,1)]*(obj.tsim/max(length(obj.tot_pw_budget)-1,1)),[obj.tot_pw_budget], '--', 'LineWidth',1,  'Color', '#CC0000');
			ax2 = subplot(3,4,[5:8]);
			p = plot(t1,u1*obj.VDom, '.');
			grid on,title('Domain Power'),ylabel('[W]'),xlabel('Time [s]'),hold on;
			if ((nargin >=4) && (isempty(u2)==0))
				p = plot(t1,u2*obj.VDom);
			end
			for pl=1:obj.vd
				plot([0:max(size(obj.quad_pw_budget,1)-1,1)]*(obj.tsim/max(size(obj.quad_pw_budget,1)-1,1)), [obj.quad_pw_budget], '--', 'LineWidth',1,  'Color', p(pl).Color);
			end
			
			u1 = u1(2:end,:);	
			dim = size(u1,1);
			%dim2 = dim
			
			pbt = repelem(obj.tot_pw_budget(min(2, length(obj.tot_pw_budget)):end),max(round(size(u1,1)/max(length(obj.tot_pw_budget)-1,1)),1),1);
			pbq = repelem(obj.quad_pw_budget(min(2, size(obj.quad_pw_budget,1)):end,:),max(round(size(u1,1)/max(size(obj.quad_pw_budget,1)-1,1)),1),1 );
			gr1 = [sum(sum(u1,2) > pbt)*obj.Ts; ...
					sum( (u1*obj.VDom) > pbq )'*obj.Ts];
			
			ttqt = [sum( (sum(u1,2) > pbt)>0 ); sum( ((u1*obj.VDom) > pbq )>0 )'];
			ttqt(ttqt<1) = 1;
			gp1 = [ sum((sum(u1,2) - pbt).*(sum(u1,2) > pbt)) / ttqt(1); ...
				sum( ((u1*obj.VDom) - pbq).*((u1*obj.VDom) > pbq))' ./  ttqt(2:end)];
			
			if (nargin < 4) || isempty(u2)
				BAR1 = [gr1(1) 0];
				BAR2 = [0 gp1(1)];
			else
				%TODO:
				%u2 = u2(2:end,:);
				%gr2 = [sum(sum(u2,2) > obj.tot_pw_budget)*obj.Ts; sum( (u2*obj.VDom)' > obj.quad_pw_budget, 2 )*obj.Ts];
				%BAR1 = [gr1(1) gr2(1) 0 0];
				%BAR2 = [0 0 gp1(1) gp2(1)];
			end
			
			subplot(3,4,9);
			%yyaxis left
			%axl = gca; % current axes
			b = bar(1,BAR1);
			%axl.Color = 'none';
			grid on,xlabel('Total'),ylabel('Time [s]');

			for pb = 1:2:floor(length(b)/2)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData/(dim*obj.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end			
			
			yyaxis right
			b2 = bar(1, BAR2);
			grid on,ylabel('\DeltaPw [W]');
			%axr = gca;
			%bar(BAR/(dim*obj.Ts)*100, 'FaceColor', "#0072BD");
			%axr.Color = 'none';
			%axr.YTick = axl.YTick/(dim*obj.Ts)*100;
			%ylabel('[%]');

			
			subplot(3,4,[10:12]);
			if (nargin < 4) || isempty(u2)
				BAR1 = [gr1(2:end) zeros(length(gr1)-1,1)];
				BAR2 = [zeros(length(gp1)-1,1) gp1(2:end)];
			else
				%TODO:
				%BAR = [gr1(2:end) gr2(2:end)];
			end
			%yyaxis left
			%axl = gca; % current axes
			b = bar(1:obj.vd, BAR1);
			%axl.Color = 'none';
			grid on,xlabel('Domain'),ylabel('Time [s]');
			
			for pb = 1:2:floor(length(b)/2)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData/(dim*obj.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			
			yyaxis right
			b2 = bar(1:obj.vd, BAR2);
			grid on,ylabel('\DeltaPw [W]');
			
			%axr = gca;
			%bar(BAR/(dim*obj.Ts)*100, 'FaceColor', "#0072BD");
			%axr.Color = 'none';
			%axr.YTick = axl.YTick/(dim*obj.Ts)*100;
			%ylabel('[%]');

		end
		function [fig] = tempconstrplot(obj,x1,x2)
			
			x1 = x1(2:end,:);
			dim = size(x1,1);
			%dim2 = dim
			gr1=obj.C(1:obj.Nc,:)*sum(x1 > obj.core_crit_temp)';
			tt1 = sum(x1 > obj.core_crit_temp);
			tt1(tt1==0) = 1;
			gt1=obj.C(1:obj.Nc,:)*(sum( (x1 - obj.core_crit_temp) .* (x1 > obj.core_crit_temp)) ./ tt1 )';

			if (nargin < 3) || isempty(x2)
				BAR1 = [sum(obj.Ts*gr1)/obj.Nc 0];
				BAR2 = [0 sum(gt1)/sum(tt1>1)];
			else
				x2 = x2(2:end,:);
				gr2=obj.C(1:obj.Nc,:)*sum(x2 > obj.core_crit_temp)';
				tt2 = sum(x2 > obj.core_crit_temp);
				tt2(tt2==0) = 1;
				gt2=obj.C(1:obj.Nc,:)*(sum( (x2 - obj.core_crit_temp) .* (x2 > obj.core_crit_temp)) ./ tt2 )';
				BAR1 = [sum(obj.Ts*gr1)/obj.Nc sum(obj.Ts*gr2)/obj.Nc 0 0 ];
				BAR2 = [0 0 sum(gt1)/sum(tt1>1) sum(gt2)/sum(tt2>1)];
			end			
			
			fig = figure('Name', 'Temperature Constraints Compliance');
			movegui(fig, 'southwest');
			title('Temperature Constraints analysis');
			subplot(5,4,[1,5,9])			
			%yyaxis left
			b = bar(1,BAR1);
			ylabel('Time [s]'),xlabel('Average'),grid on;
			for pb = 1:2:floor(length(b)/2)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData/(dim*obj.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			yyaxis right
			b2 = bar(1,BAR2);
			yl = ylabel('\DeltaT_{crit} [°C]');grid on;
			%yl.Position(2) = 0;
			%BAR = [(sum(gr1)/obj.Nc/(dim)*100)' (sum(gr2)/obj.Nc/(dim)*100)'];
			%bar(BAR,'FaceColor', "#0072BD");
			%ylabel('[%]');

			subplot(5,4,[2:4, 6:8, 10:12])
			if (nargin < 3) || isempty(x2)
				BAR1 = [obj.Ts*gr1 zeros(obj.Nc,1)];
				BAR2 = [zeros(obj.Nc,1) gt1];
			else
				BAR1 = [obj.Ts*gr1 obj.Ts*gr2 zeros(obj.Nc,1) zeros(obj.Nc,1)];
				BAR2 = [zeros(obj.Nc,1) zeros(obj.Nc,1) gt1 gt2];
			end
			%yyaxis left
			b = bar(1:obj.Nc, BAR1);
			ylabel('Time [s]'),xlabel('Core'),grid on;
			%{
			for pb = 1:2:floor(length(b)/2)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData/(dim*obj.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			%}
			yyaxis right
			b2 = bar(1:obj.Nc, BAR2);
			ylabel('\DeltaT_{crit} [°C]'),grid on;
			%BAR = [(gr1/(dim)*100)' (gr2/(dim)*100)'];
			%bar(BAR,'FaceColor', "#0072BD");
			%ylabel('[%]');
			%total_temp_violation = (obj.Ts*sum(gr))
			
			subplot(5,4,[13:20])
			if (nargin < 3) || isempty(x2)
				BAR1 = [obj.C(1:obj.Nc,:)*max(x1,[],1)' - 273.15, zeros(obj.Nc,1)];
				BAR2 = [zeros(obj.Nc,1), obj.C(1:obj.Nc,:)*mean(x1)' - 273.15];
				ymlj = min(x1, [], "all");
			else
				BAR = [obj.C(1:obj.Nc,:)*max(x1,[],1)' obj.C(1:obj.Nc,:)*max(x2,[],1)'] - 273.15;
				ymlj = min(min(x1, [], "all"),  min(x2, [], "all"));
			end
			
			b = bar(1:obj.Nc, BAR1);
			yl = ylim;
			ylim([min(min(obj.x_init), ymlj)-273.15, yl(2)]);
			ylabel('T_{Max} [°C]'),xlabel('Core'),grid on;
			hold on;
			plot(xlim, [obj.core_crit_temp obj.core_crit_temp]-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			
			ylp = ylim;
			
			yyaxis right
			b2 = bar(1:obj.Nc, BAR2);
			ylabel('T_{Average} [°C]'),grid on;
			ylim(ylp);

		end
		function [fig] = fvplot(obj, f, v)
			
			fig = figure('Name', 'Frequency-Voltage-Domains Graphs');

			ts = obj.tsim /length(f);
			t = [0:1:length(f)-1]*ts;
			
			% TODO: NOT FULLY TESTED!
			rp = 2;
			cp = 2;
			taken = 1;
			stop = 0;
			for i=1:obj.vd
				rp = 1+i;
				if (mod(rp,4)==0) || (mod(rp,3)==0) || (rp == 2)
					taken = min([ceil(rp/3) ceil(rp/2) ceil(rp/4)]);
					for j=2:2:rp			
						if (rp-taken)*j >= obj.vd
							stop = 1;
							cp = j;
							break;
						end			
					end		
				end
				if (stop)
					break;
				end
			end
			
			posif = [];
			posiv = [];
			for i=1:taken
				ind = (i-1)*cp+1;
				posif = [posif, ind:(ind-1)+cp/2];
				posiv = [posiv ind+cp/2:(ind-1)+cp];
			end
			
			ax1 = subplot(rp,cp,posif);
			plot(t, f);
			grid on,xlabel('Time [s]'), ylabel('Applied Freq [GHz]');
			xlim([0 obj.tsim]);
			
			ax2 = subplot(rp,cp,posiv);
			plot(t, v);
			grid on,xlabel('Time [s]'), ylabel('Applied Voltage [V]');
			xlim([0 obj.tsim]);
			
			axi = [];
			dim = sum(obj.VDom);
			for d=1:obj.vd
				axi = [axi subplot(rp,cp,taken*cp+d)];
				hold on;
				idx = obj.VDom(:,d).*[1:1:obj.Nc]';
				idx = idx(idx>0);
				plot(t, f(:,idx));
				idx = sum((v(:,d) * ones(1,obj.FV_levels)) > obj.FV_table(:,1)',2) + 1;
				plot(t, obj.FV_table(idx,3), '--', 'LineWidth',1, 'Color', "#0000FF");
				grid on,xlabel('Time [s]'), ylabel(strcat("Applied Freq [GHz] - Domain: ", int2str(d)));
				xlim([0 obj.tsim]);
			end

			linkaxes([ax1,ax2, axi],'x');
		end
		function [fig] = perfplot(obj, f, w)

			font_size = 6;
			
			fig = figure('Name', 'Performance Comparison');
			movegui(fig, 'northeast');
			
			ax1 = subplot(3,4,[1:4]);
			smref = round((size(f,1)-1)/(size(obj.frtrc,1)-1));
			smf = round((size(obj.frtrc,1)-1)/(size(f,1)-1));
			
			gr = repelem(f(2:1:end,:),max(smf,1),1) - repelem(obj.frtrc(2:end,:),max(smref,1),1);
			plot(obj.tsim*[0:(1/size(gr,1)):1]', [zeros(1,obj.Nc); gr]);
			grid on,xlabel('Time [s]'), ylabel('Reference Difference [GHz]');
			
			subplot(3,4,[5:6]);
			av = sum(obj.frtrc(2:end,:)) / size(obj.frtrc(2:end,:),1);
			b = bar(sum(gr)/size(gr,1)./av*100);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom', 'FontSize', font_size);
			end
			grid on,xlabel('Cores'), ylabel('Mean Difference [%]');
			
			subplot(3,4,7);
			%av = sum(obj.frtrc(2:end,:)) / size(obj.frtrc(2:end,:),1);
			b = bar(sum(gr)*obj.VDom/size(gr,1)./(av*obj.VDom)*100);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			grid on,xlabel('Domains'), ylabel('Mean Difference [%]');
			
			subplot(3,4,8);
			%av = sum(obj.frtrc(2:end,:)) / size(obj.frtrc(2:end,:),1);
			b = bar(sum(sum(gr))/size(gr,1)/sum(av)*100);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			grid on,xlabel('Total'), ylabel('Total Mean Difference [%]');
			
			% % %
			gr = w / (size(obj.wltrc,3)-1) * 100;
			subplot(3,4,[9:10]);
			b = bar(gr);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom', 'FontSize',font_size);
			end
			grid on,xlabel('Cores'), ylabel('Application Completion [%]');
			
			subplot(3,4,11);
			grv = diag(gr)*obj.VDom;
			tt = grv>0;
			min_grv = [];
			for v=1:obj.vd
				min_grv = [min_grv min(grv(tt(:,v),v))];
			end
			b = bar([obj.VDom'*gr ./ sum(obj.VDom)' min_grv']);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom', 'FontSize',font_size);
			end
			if obj.vd == 1
				b.FaceColor = 'flat';
				b.CData(2,:) = [0.8500 0.3250 0.0980];
				xticklabels({'Mean', 'Min'});
			end
			grid on,xlabel('Domains'), ylabel('Application Completion (Mean - Min) [%]');
			
			subplot(3,4,12);
			b = bar([sum(gr)/obj.Nc min(min_grv)]);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			b.FaceColor = 'flat';
			b.CData(2,:) = [0.8500 0.3250 0.0980];
			grid on,xlabel(''), ylabel('Application Completion [%]'), xticklabels({'Mean', 'Min'});
			
			%linkaxes([ax1,ax2,ax3],'x');
			
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

