classdef hpc_lab < handle
	%HPC_LAB Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		%tm;
		tasim = 2; % Free Simulation Time in [s]
		tsim = 1;					% Simulation Time in [s]
		t_init;
        temp_amb (1,1) {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 25.0 + 273.15;	% External Ambient Temperature

		
		Ts_target (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 1e-3;			% Commands min Ts, old Ts_input

		

		graph_show = 1;
		graph_save = 0;

		compare_vs_baseline = 0;

		core_pm;

		%x_init;					% Initial Conditions 
		%urplot;						% Input Reference Plot
		frtrc = {0};						% Freq Reference Trace
		%zrplot;					% Power Noise Plot
		wltrc = {0};						% Wokrload Trace

		%% All Controllers
        %TODO: this is uncommented at the moment to fix other thingy
        %tot_pw_budget; 
        %this is its substitute
        toto_pw_budget;
		chip_pw_budget;
		quad_pw_budget;
		%min_pw_red;


	end

	properties(SetAccess=protected, GetAccess=public)	

		% Internal Variables for Simulation() fnc
	 	B_s;

		osunix;

		perf_max_check;
		pmc_need_update;
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

			if isunix
				obj.osunix = 1;
			else
				obj.osunix = 0;
			end

			obj.perf_max_check{1} = 100;
			obj.pmc_need_update = 0;

			obj.get_size(obj);
        end        
        function obj = anteSimCheckLab(obj, chip)
			% Check if the floorplan has been updated
			lastwarn('');
			chip.VDom = chip.VDom;
			[msgstr, ~] = lastwarn;
			if ~isempty(msgstr)
				warning("[LAB] Incorrect dimension of VDom (voltage domains configuration)!")
				prompt = {['Do you want to run default_VDom_config() to create a default domain config? (1=yes, 0=no)' newline 'ATTENTION! This will OVERWRITE previous values.']};
				dlgtitle = '[LAB] HPC Lab Simulation';
				fieldsize = [1 45];
				definput = {'1'};
				usrin = inputdlg(prompt,dlgtitle,fieldsize,definput);
				if usrin{1} == '1'
					chip.VDom = obj.default_VDom_config(chip.Nh, chip.Nv, chip.vd);
				else
					error("[LAB] VDom has not the correct dimension");
				end				
			end
        end
		[cpxplot, cpuplot, cpfplot, cpvplot, wlop] = simulation(obj, ctrl, chip, CM,show)
        %Move This:
		function obj = taas_fix(obj, wl)
			obj.pmc_need_update = zeros(size(obj.pmc_need_update));
            if ~iscell(wl)
                error("[LAB] wl is not a cell");
            else
                obj.perf_max_check = wl;
            end
			
		end
	end
	methods(Static)
        %TODO: not working, need to fix!
        function VDom = default_VDom_config(Nh, Nv, vd)
            %above part need to be tested
			vNh = Nh;
			vNv = Nv;
			cuts = [1 1];
            vddone = (cuts(1)+1)*(cuts(2)+1); %0
			while (vddone < vd)
			
				fr = vNh>vNv;
				cuts = cuts + [fr, ~fr];
				vNh = Nh/cuts(1);
				vNv = Nv/cuts(2);
			
				tt = cuts>[Nh, Nv];
				cuts = tt.*[Nh, Nv] + ~tt.*cuts;
				vNh = vNh*(vNh>1) + (vNh<=1);
				vNv = vNv*(vNv>1) + (vNv<=1);
			
				if (cuts(2)==Nv) && (cuts(1)==Nh)
					break;
				end
			
				vddone = (cuts(1)+1)*(cuts(2)+1);
			end
			
            % this needs to be fixed, problably you should define a
            %   sequence of cuts and then have a cut index
			steps = floor([Nh, Nv] ./ (cuts+1));
			tbfd = vddone > vd;
			
			VDom = zeros(Nh*Nv, vd);
			for i=1:Nh
				for j=1:Nv
					dom = floor((i-1)/steps(1)) + floor((j-1)/steps(2))*cuts(1) + 1;
					%fixing odd vd
					dom = dom + (dom>1)*tbfd*(-1);
					dom = dom*(dom<=vd) + vd*(dom>vd);
					idx = (i-1)*Nv + j; %(i + (j-1)*Nh);
					VDom(idx, dom) = 1;
				end
			end
		end
		function totSize = get_size(class) 
			props = properties(class); 
			totSize = 0;
			
			for ii=1:length(props) 
  			currentProperty = getfield(class, char(props(ii))); 
  			s = whos('currentProperty'); 
  			totSize = totSize + s.bytes; 
			end
			
			disp(strcat("Size: ", num2str(totSize), " bytes")); 
		end
		function g=gaussian_filter(Filter_size, sigma)
			%size=5; %filter size, odd number
			size=Filter_size;
			g=zeros(size,size); %2D filter matrix
			%sigma=2; %standard deviation
			%gaussian filter
			for i=-(size-1)/2:(size-1)/2
    			for j=-(size-1)/2:(size-1)/2
        			x0=(size+1)/2; %center
        			y0=(size+1)/2; %center
        			x=i+x0; %row
        			y=j+y0; %col
        			g(y,x)=exp(-((x-x0)^2+(y-y0)^2)/2/sigma/sigma);
    			end
			end
			%normalize gaussian filter
			sum1=sum(g);
			sum2=sum(sum1);
			g=g/sum2;

			%plot 3D
			%g1=Gaussian_filter(50,2);
			%g2=Gaussian_filter(50,7);
			%g3=Gaussian_filter(50,11);
			%figure(1);
			%subplot(1,3,1);surf(g1);title('filter size = 50, sigma = 2');
			%subplot(1,3,2);surf(g2);title('filter size = 50, sigma = 7');
			%subplot(1,3,3);surf(g3);title('filter size = 50, sigma = 11');
		end
    end

	%% Approximation
	methods(Static)
		[k0, k1, k2] = pws_ls_approx(chip, I, T, Toff, C, alp, alp_I0, using_voltage)
		[lut, F, T, M_var] = pws_ls_offset(chip, ctrl, Vslot, Tslot, Tmin, Tmax, show )
    end

    %% Populate Controllers
    %methods
    %end

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
		function [fig] = xutplot(obj,chip,x,u)
			t = [0:1:length(x)-1]*chip.Ts;
			fig = figure('Name', 'T-P Evolution');
			movegui(fig, 'northwest');
			sbp = 1 + chip.full_model_layers + (chip.add_states>0);
			%TODO: make it parametric. Also make the names of add_state
			%		parametric in thermal_model
			ax1=subplot(sbp,1,1);
			plot(t, x(:,1+end-chip.add_states:end)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Others [°C]'), hold on;
			legend("Heat-sink", "PCB", "Motherboard", "Air");
			%
			ax2=subplot(sbp,1,2);
			plot(t, x(:,2:2:end-chip.add_states)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Heat Spreader (Cu) [°C]'),hold on;
			%
			ax3=subplot(sbp,1,3);
			plot(t, x(:,1:2:end-chip.add_states)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Cores (Si) [°C]'), hold on;
			if max(max(x(:,1:2:end-1))) >= (chip.core_limit_temp-(0.1*(chip.core_limit_temp-273)))
				yline(chip.core_limit_temp-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			end
			%
			ax4=subplot(sbp,1,sbp);
			plot(t, u);
			grid on,xlabel('Time [s]'),ylabel('Cores Power Output [W]'), hold on;
			if max(max(u)) >= (chip.core_max_power-(0.34*chip.core_max_power))
				yline(chip.core_max_power, '--', 'Max Power', 'LineWidth',1,  'Color', 'b');
			end
			if min(min(u)) <= (chip.core_min_power+(0.34*chip.core_max_power))
				yline(chip.core_min_power, '--','min Power', 'LineWidth',1,  'Color', 'b');
			end
			
			linkaxes([ax1,ax2,ax3,ax4],'x');
		end
        function [fig] = powerconstrplot(obj,chip,idx,u1,u2)
						
			fig = figure('Name', 'Power Constraints Compliance');
			movegui(fig, 'southeast');
			title('Power Constraints analysis')
			
			t1 = [0:1:length(u1)-1]*chip.Ts;
			
			ax1 = subplot(3,4,[1:4]);
			gr = sum(u1,2);
			plot(t1,gr, '.');
			grid on,title('Total Power'),ylabel('[W]'),xlabel('Time [s]'),hold on;
			if ((nargin >=5) && (isempty(u2)==0))
				gr = sum(u2,2);
				plot(t1,gr);
			end
			%if max(gr) >= obj.usum(1)*0.2
			plot([0:max(length(obj.chip_pw_budget{idx})-1,1)]*(obj.tsim/max(length(obj.chip_pw_budget{idx})-1,1)),[obj.chip_pw_budget{idx}], '--', 'LineWidth',1,  'Color', '#CC0000');
			ax2 = subplot(3,4,[5:8]);
			p = plot(t1,u1*chip.VDom, '.');
			grid on,title('Domain Power'),ylabel('[W]'),xlabel('Time [s]'),hold on;
			if ((nargin >=5) && (isempty(u2)==0))
				p = plot(t1,u2*chip.VDom);
			end
			for pl=1:chip.vd
				plot([0:max(size(obj.quad_pw_budget{idx},1)-1,1)]*(obj.tsim/max(size(obj.quad_pw_budget{idx},1)-1,1)), [obj.quad_pw_budget{idx}], '--', 'LineWidth',1,  'Color', p(pl).Color);
			end
			
			u1 = u1(2:end,:);	
			dim = size(u1,1);
			%dim2 = dim
			
			pbt = repelem(obj.chip_pw_budget{idx}(min(2, length(obj.chip_pw_budget{idx})):end),max(round(size(u1,1)/max(length(obj.chip_pw_budget{idx})-1,1)),1),1);
			pbq = repelem(obj.quad_pw_budget{idx}(min(2, size(obj.quad_pw_budget{idx},1)):end,:),max(round(size(u1,1)/max(size(obj.quad_pw_budget{idx},1)-1,1)),1),1 );
			gr1 = [sum(sum(u1,2) > pbt)*chip.Ts; ...
					sum( (u1*chip.VDom) > pbq )'*chip.Ts];
			
			ttqt = [sum( (sum(u1,2) > pbt)>0 ); sum( ((u1*chip.VDom) > pbq )>0 )'];
			ttqt(ttqt<1) = 1;
			gp1 = [ sum((sum(u1,2) - pbt).*(sum(u1,2) > pbt)) / ttqt(1); ...
				sum( ((u1*chip.VDom) - pbq).*((u1*chip.VDom) > pbq))' ./  ttqt(2:end)];
			
			if (nargin < 5) || isempty(u2)
				BAR1 = [gr1(1) 0];
				BAR2 = [0 gp1(1)];
			else
				%TODO:
				%u2 = u2(2:end,:);
				%gr2 = [sum(sum(u2,2) > obj.chip_pw_budget{idx})*chip.Ts; sum( (u2*chip.VDom)' > obj.quad_pw_budget{idx}, 2 )*chip.Ts];
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
				labels = string(round(b(pb).YData/(dim*chip.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end			
			
			yyaxis right
			b2 = bar(1, BAR2);
			grid on,ylabel('\DeltaPw [W]');
			%axr = gca;
			%bar(BAR/(dim*chip.Ts)*100, 'FaceColor', "#0072BD");
			%axr.Color = 'none';
			%axr.YTick = axl.YTick/(dim*chip.Ts)*100;
			%ylabel('[%]');

			
			subplot(3,4,[10:12]);
			if (nargin < 5) || isempty(u2)
				BAR1 = [gr1(2:end) zeros(length(gr1)-1,1)];
				BAR2 = [zeros(length(gp1)-1,1) gp1(2:end)];
			else
				%TODO:
				%BAR = [gr1(2:end) gr2(2:end)];
			end
			%yyaxis left
			%axl = gca; % current axes
			b = bar(1:chip.vd, BAR1);
			%axl.Color = 'none';
			grid on,xlabel('Domain'),ylabel('Time [s]');
			
			for pb = 1:2:floor(length(b)/2)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData/(dim*chip.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			
			yyaxis right
			b2 = bar(1:chip.vd, BAR2);
			grid on,ylabel('\DeltaPw [W]');
			
			%axr = gca;
			%bar(BAR/(dim*chip.Ts)*100, 'FaceColor', "#0072BD");
			%axr.Color = 'none';
			%axr.YTick = axl.YTick/(dim*chip.Ts)*100;
			%ylabel('[%]');

		end
		function [fig] = tempconstrplot(obj,chip,x1,x2)
			
			x1 = x1(2:end,:);
			dim = size(x1,1);
			%dim2 = dim
			gr1=chip.C(1:chip.Nc,:)*sum(x1 > chip.core_limit_temp)';
			tt1 = sum(x1 > chip.core_limit_temp);
			tt1(tt1==0) = 1;
			gt1=chip.C(1:chip.Nc,:)*(sum( (x1 - chip.core_limit_temp) .* (x1 > chip.core_limit_temp)) ./ tt1 )';

			if (nargin < 4) || isempty(x2)
				BAR1 = [sum(chip.Ts*gr1)/chip.Nc 0];
				BAR2 = [0 sum(gt1)/sum(tt1>1)];
			else
				x2 = x2(2:end,:);
				gr2=chip.C(1:chip.Nc,:)*sum(x2 > chip.core_limit_temp)';
				tt2 = sum(x2 > chip.core_limit_temp);
				tt2(tt2==0) = 1;
				gt2=chip.C(1:chip.Nc,:)*(sum( (x2 - chip.core_limit_temp) .* (x2 > chip.core_limit_temp)) ./ tt2 )';
				BAR1 = [sum(chip.Ts*gr1)/chip.Nc sum(chip.Ts*gr2)/chip.Nc 0 0 ];
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
				labels = string(round(b(pb).YData/(dim*chip.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			yyaxis right
			b2 = bar(1,BAR2);
			yl = ylabel('\DeltaT_{crit} [°C]');grid on;
			%yl.Position(2) = 0;
			%BAR = [(sum(gr1)/chip.Nc/(dim)*100)' (sum(gr2)/chip.Nc/(dim)*100)'];
			%bar(BAR,'FaceColor', "#0072BD");
			%ylabel('[%]');

			subplot(5,4,[2:4, 6:8, 10:12])
			if (nargin < 4) || isempty(x2)
				BAR1 = [chip.Ts*gr1 zeros(chip.Nc,1)];
				BAR2 = [zeros(chip.Nc,1) gt1];
			else
				BAR1 = [chip.Ts*gr1 chip.Ts*gr2 zeros(chip.Nc,1) zeros(chip.Nc,1)];
				BAR2 = [zeros(chip.Nc,1) zeros(chip.Nc,1) gt1 gt2];
			end
			%yyaxis left
			b = bar(1:chip.Nc, BAR1);
			ylabel('Time [s]'),xlabel('Core'),grid on;
			%{
			for pb = 1:2:floor(length(b)/2)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData/(dim*chip.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			%}
			yyaxis right
			b2 = bar(1:chip.Nc, BAR2);
			ylabel('\DeltaT_{crit} [°C]'),grid on;
			%BAR = [(gr1/(dim)*100)' (gr2/(dim)*100)'];
			%bar(BAR,'FaceColor', "#0072BD");
			%ylabel('[%]');
			%total_temp_violation = (chip.Ts*sum(gr))
			
			subplot(5,4,[13:20])
			if (nargin < 4) || isempty(x2)
				BAR1 = [chip.C(1:chip.Nc,:)*max(x1,[],1)' - 273.15, zeros(chip.Nc,1)];
				BAR2 = [zeros(chip.Nc,1), chip.C(1:chip.Nc,:)*mean(x1)' - 273.15];
				ymlj = min(x1, [], "all");
			else
				BAR = [chip.C(1:chip.Nc,:)*max(x1,[],1)' chip.C(1:chip.Nc,:)*max(x2,[],1)'] - 273.15;
				ymlj = min(min(x1, [], "all"),  min(x2, [], "all"));
			end
			
			b = bar(1:chip.Nc, BAR1);
			yl = ylim;
			ylim([min(min(obj.temp_amb), ymlj)-273.15, yl(2)]);
			ylabel('T_{Max} [°C]'),xlabel('Core'),grid on;
			hold on;
			plot(xlim, [chip.core_limit_temp chip.core_limit_temp]-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			
			ylp = ylim;
			
			yyaxis right
			b2 = bar(1:chip.Nc, BAR2);
			ylabel('T_{Average} [°C]'),grid on;
			ylim(ylp);

		end
		function [fig] = fvplot(obj, chip, f, v)
			
			fig = figure('Name', 'Frequency-Voltage-Domains Graphs');

			ts = obj.tsim /length(f);
			t = [0:1:length(f)-1]*ts;
			
			% TODO: NOT FULLY TESTED!
			rp = 2;
			cp = 2;
			taken = 1;
			stop = 0;
			for i=1:chip.vd
				rp = 1+i;
				if (mod(rp,4)==0) || (mod(rp,3)==0) || (rp == 2)
					taken = min([ceil(rp/3) ceil(rp/2) ceil(rp/4)]);
					for j=2:2:rp			
						if (rp-taken)*j >= chip.vd
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
			dim = sum(chip.VDom);
			for d=1:chip.vd
				axi = [axi subplot(rp,cp,taken*cp+d)];
				hold on;
				idx = chip.VDom(:,d).*[1:1:chip.Nc]';
				idx = idx(idx>0);
				plot(t, f(:,idx));
				idx = sum((v(:,d) * ones(1,chip.FV_levels)) > chip.FV_table(:,1)',2) + 1;
				plot(t, chip.FV_table(idx,3), '--', 'LineWidth',1, 'Color', "#0000FF");
				grid on,xlabel('Time [s]'), ylabel(strcat("Applied Freq [GHz] - Domain: ", int2str(d)));
				xlim([0 obj.tsim]);
			end

			linkaxes([ax1,ax2, axi],'x');
		end
        function [fig] = perfplot(obj, chip, f, w, cmp,idx)

			font_size = 6;
			
			fig = figure('Name', 'Performance Comparison');
			movegui(fig, 'northeast');
			
			ax1 = subplot(3,4,[1:4]);
			smref = round(max((size(f,1)-1),1)/max((size(obj.frtrc{idx},1)-1),1));
			smf = round(max((size(obj.frtrc{idx},1)-1),1)/max((size(f,1)-1),1));
			
			gr = repelem(f(2:1:end,:),max(smf,1),1) - repelem(obj.frtrc{idx}(2:end,:),max(smref,1),1);
			plot(obj.tsim*[0:(1/size(gr,1)):1]', [zeros(1,chip.Nc); gr]);
			grid on,xlabel('Time [s]'), ylabel('Reference Difference [GHz]');
			
			subplot(3,4,[5:6]);
			av = sum(obj.frtrc{idx}(2:end,:)) / size(obj.frtrc{idx}(2:end,:),1);
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
			b = bar(sum(gr)*chip.VDom/size(gr,1)./(av*chip.VDom)*100);
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
			%this if() is for saveall() and other plot
			if w <= 1.001
				gr = w * 100;
			else
				gr = w ./ cmp * 100;
			end
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
			grv = diag(gr)*chip.VDom;
			tt = grv>0;
			min_grv = [];
			for v=1:chip.vd
				min_grv = [min_grv min(grv(tt(:,v),v))];
			end
			b = bar([chip.VDom'*gr ./ sum(chip.VDom)' min_grv']);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom', 'FontSize',font_size);
			end
			if chip.vd == 1
				b.FaceColor = 'flat';
				b.CData(2,:) = [0.8500 0.3250 0.0980];
				xticklabels({'Mean', 'Min'});
			end
			grid on,xlabel('Domains'), ylabel('Application Completion (Mean - Min) [%]');
			
			subplot(3,4,12);
			b = bar([sum(gr)/chip.Nc min(min_grv)]);
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
		function [] = saveall(obj, xop, uop, fop, vop, wop, path, name)
			if nargin < 7 || (strlength(path)<1)
				if isunix
					path = "/tmp/MATLAB-Figures";
				else
					path = "C:\temp\MATLAB-Figures";
				end
			end
			if nargin < 8 || (strlength(name)<1)
				name = string(datetime('now','Format','dd-MMM-yyyy HH_mm_ss_SSS'));
			end

			if ~exist(path, 'dir')
				mkdir(path);
			end
			if isunix
				separator_os = "/";
			else
				separator_os = "\";
			end
			
			fig = [];
			namefig = "";
			res = 1;

			for i=1:5
				switch i
					case 1
						fig = obj.xutplot(xop, uop);
						namefig = "TP";
						res = 1.5;
					case 2
						fig = obj.powerconstrplot(uop);
						namefig = "Pw";
						res = 1;
					case 3
						fig = obj.tempconstrplot(xop);
						namefig = "T";
						res = 1;
					case 4
						fig = obj.perfplot(fop, wop / 100 );
						namefig = "Perf";
						res = 1;
					case 5 
						fig = obj.fvplot(fop,vop);
						namefig = "FV";
						res = 1.5;
				end
				if ~isempty(fig)
					path_name = strcat(path,separator_os,namefig, "-", name);
					obj.savetofile(fig, path_name, res);
				end
			end %for
		end %function
        function [fig] = leakageplot(obj, T_step, T_offset)
            if (nargin < 3) || (isempty(T_offset))
                T_offset = [10 10];
            end
            if (nargin < 2) || (isempty(T_step))
                T_step = 1;
            end

            V = [];
            T = [];
            V = [obj.V_min:obj.V_discretization_step:obj.V_max];
            d_p = 1;
            T = [obj.temp_amb-T_offset(1):T_step:obj.core_crit_temp+T_offset(2)];
            
            kps = [obj.leak_vdd_k, obj.leak_temp_k, obj.leak_process_k, ...
					            obj.leak_exp_vdd_k, obj.leak_exp_t_k, obj.leak_exp_k];
            obj.leak_exp_vdd_k = 0;
            obj.leak_exp_t_k = 0;
            obj.leak_exp_k = 0;
            
            Power_static = obj.ps_compute(V,T(1),d_p,0);
            
            obj.leak_exp_vdd_k = kps(4);
            obj.leak_exp_t_k = kps(5);
            obj.leak_exp_k = kps(6);
            
            res = [];
            plane = [];
            for i=1:length(T)
	            res(:,i) = obj.ps_compute(V,T(i),d_p,0)';%*0.95 + 0.05;
	            plane(:,i) = Power_static;
            end
            
            F = obj.FV_table( [sum(V>obj.FV_table(:,1))+1]',3);
            
            ci = obj.dyn_ceff_k;
            Ceff_low = [1 zeros(1,obj.ipl-1)] * ci';
            Ceff_high = [zeros(1,obj.ipl-1) 1] * ci';
            Power_dyn_low = Ceff_low .* F .* (V' .* V');
            Power_dyn_high = Ceff_high .* F .* (V' .* V');
            
            T = T - 273.15;
            
            %%
            fig = figure()
            
            subplot(3, 3, [1:2, 4:5]);
            
            surf(T, V, res, 'EdgeAlpha', 0.6)
            hold;
            h = gca;
            surf([T(1) T(end)], V, [Power_static(:), Power_static(:)], 'FaceColor', 'm', ... %"#A2142F", 
	            'FaceAlpha', 0.6)
            
            view([-81.6598139306647 13.4415460949037]);
            
            xlim([T(1) T(end)]);
            ylim([V(1) V(end)]);
            
            xlb = xlabel('Temperature [°C]');
            set(xlb,'rotation',49)
            ylabel('Voltage [V]');
            zlabel('Power [W]');
            ax = gca;
            ax.FontSize = 12;
            
            title("Comparison Exponential leakage (Blue) and Linear leakage (Magenta)", 'FontSize', 20);
            
            %view([-80.5895672336695 16.4241547850151])
            
            subplot(3, 3, 7)
            
            surf(T, V, res)
            hold;
            
            %surf(T, V, repelem(Power_dyn_low, 1,length(T)), 'FaceColor', 'y', ... %"#A2142F", 
            %	'FaceAlpha',0.3)
            
            surf(T, V, plane, 'FaceColor', 'm', ... %"#A2142F", 
	            'FaceAlpha',0.6)
            
            fnz= 13;
            
            xlabel('Temperature [°C]', 'FontSize', fnz);
            ylabel('Voltage [V]', 'FontSize', fnz);
            zlabel('Power [W]', 'FontSize', fnz);
            ax = gca;
            ax.FontSize = fnz;
            xlim([T(1) T(end)]);
            ylim([V(1) V(end)]);
            
            title({"Top View:"; "Comparison Exponential and Linear leakage"}, 'FontSize', 15);
            
            view(0,90)
            %view(0,0)
            %view(90,0)
            
            
            %%%figure()
            
            subplot(3, 3, 3);
            
            surf(T, V, res, 'FaceAlpha',1, 'EdgeAlpha', 0.6);%0.6)
            hold;
            
            h = gca;
            
            for i=1:length(T)
	            planel(:,i) = Power_dyn_low;
	            planeh(:,i) = Power_dyn_high;
            end
            
            %plot3(h.XLim(1):(h.XLim(2) - h.XLim(1))/(length(V)-1):h.XLim(2), V, Power_dyn_high, 'r', 'LineWidth', 1.2);
            
            for i=1:length(T)
	            plot3(T(i)*ones(size(V)), V, Power_dyn_low, 'y', 'LineWidth', 1.2);
	            plot3(T(i)*ones(size(V)), V, Power_dyn_high, 'r', 'LineWidth', 1.2);
	            plot3(T(i)*ones(size(V)), V, Power_static, 'm', 'LineWidth', 1.2);
            end
            %surf(T, V, planel, 'FaceColor', "#77AC30")
            %surf(T, V, planeh, 'FaceColor', "#EDB120", 'FaceAlpha',0.2)
            
            xlabel('Temperature [°C]');
            ylabel('Voltage [V]');
            zlabel('Power [W]');
            ax = gca;
            ax.FontSize = fnz;
            xlim([T(1) T(end)]);
            ylim([V(1) V(end)]);
            
            title({"Comparison Exponential leakge"; "with Max and min Dynamic Power"}, 'FontSize', 15);
            
            view(90,0)
            
            subplot(3, 3, [6, 9]);
            
            surf(T, V, res, 'FaceAlpha',1, 'EdgeAlpha', 0.6);%0.6)
            hold on;
            
            h = gca;
            surf([T(1) T(end)], V, [Power_dyn_low(:), Power_dyn_low(:)], 'FaceColor', 'y', 'FaceAlpha', 0.8);
            surf([T(1) T(end)], V, [Power_dyn_high(:), Power_dyn_high(:)], 'FaceColor', 'r', 'FaceAlpha', 0.8);
            xlb = xlabel('Temperature [°C]');
            set(xlb,'rotation',-35)
            ylb = ylabel('Voltage [V]');
            set(ylb,'rotation',40)
            zlabel('Power [W]');
            view([-138.844996600224 39.7336634094982]);
            ax = gca;
            ax.FontSize = fnz;
            xlim([T(1) T(end)]);
            ylim([V(1) V(end)]);
            
            %title({"Comparison Exponential leakge"; "with Max and min Dynamic Power"}, 'FontSize', 15);
            
            %%figure()
            
            lij = 14;
            
            subplot(3, 3, 8);
            surf([T(1) T(end)], V(end-lij: end), [Power_dyn_high(end-lij:end), Power_dyn_high(end-lij:end)], 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.2);
            hold on;
            surf(T, V, res, 'FaceAlpha',1, 'EdgeAlpha', 0.6);%0.6)
            
            
            %surf(T, V, planel, 'FaceColor', "#77AC30")
            %surf(T, V, planeh, 'FaceColor', "#EDB120", 'FaceAlpha',0.2)
            
            h = gca;
            
            
            
            surf([T(1) T(end)], V, [Power_dyn_low(:), Power_dyn_low(:)], 'FaceColor', 'y', 'FaceAlpha', 0.8, 'EdgeAlpha', 0.3);
            
            %{
            for i=1:length(Power_dyn_high)
	            plot3([20 100], [V(i) V(i)], [Power_dyn_low(i) Power_dyn_low(i)], 'y', 'LineWidth', 0.8);
	            plot3([20 100], [V(i) V(i)], [Power_dyn_high(i) Power_dyn_high(i)], 'r', 'LineWidth', 0.8);
            end
            %}
            
            xlabel('Temperature [°C]');
            ylabel('Voltage [V]');
            zlabel('Power [W]');
            ax = gca;
            ax.FontSize = fnz;
            xlim([T(1) T(end)]);
            ylim([V(1) V(end)]);
            
            title({"Comparison Exponential leakge"; "with Max and min Dynamic Power"}, 'FontSize', 15);
            
            view(0,0)
                        

        end
        function [fig] = discrplot(obj)
            fig = figure();
            subplot(1,2,1);
            
            Fvect = obj.FV_table(1,2):obj.F_discretization_step:obj.FV_table(end,3);
            Vvect = obj.FV_table(sum(Fvect > obj.FV_table(:,3))+1,1);
            
            Fvect = Fvect';
            
            Fvect2 = ones(obj.FV_levels,1)*Fvect';
            
            for i=1:obj.FV_levels
	            Fvect2(i,Fvect2(i,:) > obj.FV_table(i,3)) = 0;
            end
            T = ones(length(Fvect),1)*(70+273.15);
	            
            d_i = ones(length(Fvect),1)*[0,0,0.25,0.75,0];
            d_p = ones(length(Fvect),1);
            
            h = [];
            for j=1:obj.FV_levels
	            pup = obj.power_compute(Fvect2(j,:)', obj.FV_table(j,1)*ones(size(Fvect2,2),1), T, d_i,d_p, 0);
	            %pup = obj.power_compute(Fvect, Vvect, T, d_i,d_p, 0);
	            pup(Fvect2(j,:)'<=0)=NaN;
	            hp = plot(Fvect,pup,'LineWidth', 1.5);
	            h = [h hp];
	            hold on;
            end	
            pup = obj.power_compute(Fvect, Vvect, T, d_i,d_p, 0);
            hp = plot(Fvect,pup, 'LineWidth', 3, 'Color', 'b');
            h = [h hp];
            
            legend(h(end), "Max F-V line", 'FontSize',14);
            
            xlim([Fvect(1), Fvect(end)]);
            hold on;
            fnz = 16;
            xlabel('Frequency [GHz]', 'FontSize',14);
            ylabel('Power [W]','FontSize',14);
            ax = gca;
            ax.FontSize = fnz;
            grid on;
            
            title('P-F Relation @Voltage \in [0.5; 1.2]V', 'FontSize',24 );
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            subplot(1,2,2);
            
            Fvect = obj.FV_table(1,2):obj.F_discretization_step*2:obj.FV_table(end,3);
            Vvect = obj.FV_table(sum(Fvect > obj.FV_table(:,3))+1,1);
            
            hold on;
            
            Fvect = Fvect';
            
            %Fvect = Fvect3;
            %Vvect = Vvect3;
            
            T = ones(obj.FV_levels,1)*(70+273.15);
            d_i = ones(obj.FV_levels,1)*[0,0,0.25,0.75,0];
            d_p = ones(obj.FV_levels,1);
            
            Fvect2 = ones(obj.FV_levels,1)*Fvect';
            for i=1:obj.FV_levels
	            Fvect2(i,Fvect2(i,:) > obj.FV_table(i,3)) = 0;
            end
            
            h =[];
            for j=1:size(Fvect2,2)
	            pup = obj.power_compute(Fvect2(:,j), obj.FV_table(:,1), T, d_i,d_p, 0);
	            %pup = obj.power_compute(Fvect, Vvect, T, d_i,d_p, 0);
	            pup(Fvect2(:,j)'<=0)=NaN;
	            hp = plot(obj.FV_table(:,1),pup, 'LineWidth', 1.5);
	            h = [h hp];
	            hold on;
            end	
            
            %Vvect3 = obj.FV_table(:,1);
            %Fvect3 = obj.FV_table(:,3);
            Fvect3 = obj.FV_table(1,2):obj.F_discretization_step:obj.FV_table(end,3);
            Vvect3 = obj.FV_table(sum(Fvect3 > obj.FV_table(:,3))+1,1);
            T = ones(length(Vvect3),1)*(70+273.15);
            d_i = ones(length(Vvect3),1)*[0,0,0.25,0.75,0];
            d_p = ones(length(Vvect3),1);
            
            pup = obj.power_compute(Fvect3', Vvect3, T, d_i,d_p, 0);
            hp = plot(Vvect3,pup, 'LineWidth', 3, 'Color', 'b');
            h = [h hp];
            
            legend(h(end), "Max F-V line", 'FontSize',14)
            
            xlim([Vvect3(1), Vvect3(end)]);
            
            xlabel('Voltage [V]', 'FontSize',14);
            ylabel('Power [W]', 'FontSize',14);
            ax = gca;
            ax.FontSize = fnz;
            title('P-V Relation @Frequency \in [0.4; 3.6]GHz', 'FontSize',24);
            
            grid on;
        end
    end

	%% Result Analysis
	methods
		res = stats_analysis(obj,ctrl,chip,x,u,f,v,w)
		[wlres, wlop] = base_ideal_unr(obj)
	end

	%% Libraries
	methods(Static)
		[p,n] = numSubplots(n)
	end

	%% Dependent Variables
	methods
		function obj = set.compare_vs_baseline(obj, val)
			if val && (~obj.compare_vs_baseline)
				obj.pmc_need_update(:) = 1;
                for i=1:length(obj.wltrc)
				    obj.perf_max_check{i} = size(obj.wltrc{i},3)-1;
                end
			end
			obj.compare_vs_baseline = val;
		end
		function obj = set.frtrc(obj, val)
            %check which index I'm modifing
            oldval = obj.frtrc;
            if ~iscell(val)
                val{1} = val;
            end
            lenold = length(oldval);
            lennew = length(val);
            %changed = zeros(lenold,1);
            for i=1:length(lenold)
                %ll = size(val{i});
                %changed(i) = oldval{i}~=val{i};
                obj.pmc_need_update(i) = 1;
                obj.perf_max_check{i} = size(obj.wltrc{i},3)-1;
            end
			%obj.pmc_need_update = obj.pmc_need_update || changed;
            if lenold ~= lennew
                obj.pmc_need_update(lenold+1:lennew) = 1;
                for i=lenold+1:lennew
                    if (length(obj.wltrc) < lennew)
                        obj.perf_max_check{i} = 100;
                    else
                        obj.perf_max_check{i} = size(obj.wltrc{i},3)-1;
                    end
                end
            end
			
			obj.frtrc = val;
		end
		function obj = set.wltrc(obj, val)
            %check which index I'm modifing
            oldval = obj.wltrc;
            if ~iscell(val)
                val{1} = val;
            end
            lenold = length(oldval);
            lennew = length(val);
            %changed = zeros(lenold,1);
            for i=1:length(lenold)
                %ll = size(val{i});
                %changed(i) = oldval{i}~=val{i};
                obj.pmc_need_update(i) = 1;
                obj.perf_max_check{i} = size(val{i},3)-1;
            end
			%obj.pmc_need_update = obj.pmc_need_update || changed;
            if lenold ~= lennew
                obj.pmc_need_update(lenold+1:lennew) = 1;
                for i=lenold+1:lennew
                    obj.perf_max_check{i} = size(val{i},3)-1;
                end
            end

			obj.wltrc = val;
		end
	end
end

