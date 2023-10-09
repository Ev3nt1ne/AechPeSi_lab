function [] = simulate_aut(obj, ts, Nsample, exp_gamma)

	if (nargin < 2) || isempty(ts)
		ts = obj.tesim;
	end
	if (nargin < 3) || isempty(Nsample)
		Nsample = 5;
	end
	if (nargin < 4) || isempty(exp_gamma)
		exp_gamma = 4.5;
	end
	
	show = 1;
	save = 1;
	
	x_init = obj.x_init;

	N = ceil(ts / obj.Ts);

	Adl_true = obj.Ad_true;
	Bdl_true = obj.Bd_true;	
	d_p = ones(obj.Nc,1);
	
	xim = 0:(1/(Nsample+1)):1;
	yim = xim.^(exp_gamma);

	pos = floor(yim * N);
	
	pos(1) = 1;

	for i=2:length(pos)
		if pos(i) <= pos(i-1)
			pos(i) = pos(i-1)+1;
		end
	end
	
	test_num = 4;
	xfp = zeros(length(pos), obj.Ns, test_num);
	spnum = obj.numSubplots(length(pos));
	figname = [];
	
	for tst=1:test_num
		x = obj.x_init;
		id = 1;
		xp = zeros(N+1, obj.Ns);
		xp(1,:) = x_init;
		pp = zeros(N+1, obj.Nc);
			
		switch tst
			case 1
				F = ones(obj.Nc,1)*obj.F_min;
				wl = ones(obj.Nc,1)*[1 zeros(1,obj.ipl-1)];
				figname = [figname "Temperature Gradient - P min"];
			case 2
				F = ones(obj.Nc,1)*(obj.F_Max*4+obj.F_min)/5;

				wpid = floor(obj.ipl/2);
				wl = zeros(obj.Nc,obj.ipl);
				wl(:, [wpid, wpid+1]) = 0.5;
				figname = [figname "Temperature Gradient - P med"];
			case 3
				F = ones(obj.Nc,1)*obj.F_Max;
				wl = ones(obj.Nc,1)*[zeros(1,obj.ipl-1) 1];
				figname = [figname "Temperature Gradient - P Max"];
			case 4
				midcore = ceil(obj.Nh/2)*obj.Nv + ceil(obj.Nv/2);
				F = ones(obj.Nc,1)*obj.F_min;
				wl = ones(obj.Nc,1)*[1 zeros(1,obj.ipl-1)];
				F(midcore) = obj.F_Max;
				wl(midcore, :) = [zeros(1,obj.ipl-1) 1];
				figname = [figname "Temperature Dispersion"];
		end
		
		%here it is just a test to check the thermal model behavior, so,
		%	INDEPENDENTLY from the domain, I will make 1 Vdd per core
		V = obj.FV_table(sum(F > (obj.FV_table(:,3)*ones(1,obj.Nc))',2)+1,1);
	
		for s=1:N			
			u = obj.power_compute(F,V,obj.C(1:obj.Nc,:)*x,wl,d_p);
			x = Adl_true*x + Bdl_true*[u;obj.temp_amb*1000];
			xp(s+1,:) = x;
			pp(s+1,:) = u;
			if sum(s==pos)
				xfp(id,:,tst) = x-273.15;
				id = id + 1;
			end
		end
		
		pp(1,:) = pp(2,:);
		if show
			fig = obj.xutplot(xp, pp);
			%fig.Position = [1         865        1920        1080];
		end
		
		if save && ((tst==2) || (tst==4))
			obj.savetofile(fig, strcat("C:\temp\MATLAB-Figures\", figname(tst)), 1);
		end
	end
	
	%% show
	
	if show	
		time = ts/N;
		%figname
		%xfp

		%TODO parametrize 25 and 130
		hclrlim = [ max(min(xfp,[],'All'), 25) min(max(xfp,[],'All'), 130)];	

		for tst=1:length(figname)		
			fig = figure('Name', figname(tst)); %, 'Position', get(0, 'Screensize'));
			tcl = tiledlayout(fig,spnum(1),spnum(2));
			n = size(xfp,1);
			h = gobjects(n,1); 

			for f=1:n
				ax = nexttile(tcl);
				cc = reshape(obj.C(1:obj.Nc,:)*squeeze(xfp(f,:,tst))', obj.Nv, obj.Nh)';
				h(f) = heatmap(cc);

				h(f).Title = strcat('T = ', sprintf('%.5f',(pos(f)*time)));
				h(f).Colormap = obj.customColormap; %turbo; %parula;
				h(f).ColorLimits = hclrlim;
				h(f).ColorbarVisible = 'off';
			end	

			% Create global colorbar that uses the global color limits
			ax = axes(tcl,'visible','off','Colormap',h(1).Colormap,'CLim',hclrlim);
			cb = colorbar(ax);
			cb.Layout.Tile = 'East';

			%fig.Position = [1921    220.2    2560    1440];

			if save && ((tst==2) || (tst==4))
				obj.savetofile(fig, strcat("C:\temp\MATLAB-Figures\", figname(tst), "_map"), 1);
			end

		end
	end
	

	
	%{	
	axpc = axes(fig);
	c = colorbar(axpc,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
	%c.Colormap = obj.customColormap;
	%caxis(axpc,hclrlim);             % set colorbar limits
	
	%%

	hcr = axes(fig);
	cbh = colorbar; 
	% Reposition to figure's left edge, centered vertically
	cbh.Position(1) = .95-cbh.Position(3);
	cbh.Position(2) = 0.5-cbh.Position(4)/2;
	% decrease horizontal extent of subplots to 92% of their current width
	set(hcr, {'Position'}, mat2cell(vertcat(hcr.Position) .* [1 1 .92, 1], ones(size(hcr(:))),4))
	%}

	
	
	
	
	
	
	
	
	
	
end