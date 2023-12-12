%%

WL_N = 3;
old_perf = 1;

p1 = [   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   94.2399
   ];

p2 = [   68.7267
   10.9303
   10.9303
   94.2399
   10.9303
   10.9303
   10.9303
   94.2399
   10.9303
   68.7267
   68.7267
   94.2399
   10.9303
   68.7267
   10.9303
   68.7267
   94.2399
   10.9303
   10.9303
   94.2399
   10.9303
   94.2399
   10.9303
   68.7267
   10.9303
   10.9303
   94.2399
   10.9303
   68.7267
   10.9303
   94.2399
   68.7267
   10.9303
   68.7267
   10.9303
   94.2399];

%p2 = p2(coreid12);

p3 = [   91.2498
   88.9647
   87.6472
   90.7673
   88.1172
   89.4897
   88.7772
   89.0797
   88.9622
   89.5697
   90.2273
   88.7547
   88.6847
   88.2572
   90.0673
   88.8922
   88.5272
   88.0022
   89.2972
   87.3722
   88.4622
   89.3922
   88.1097
   89.5147
   88.3922
   88.8197
   87.2572
   88.7122
   90.2573
   87.3047
   89.3597
   86.9047
   87.7922
   89.0322
   88.6722
   90.6123
   ];


perf_max_check{1} = p1;
perf_max_check{2} = p2;
perf_max_check{3} = p3;


%%

%tres{di, mdli, 2, wli} = hpc.stats_analysis(xop, uop, fop, vop, wlop);

needs2do = 1;
if needs2do
wlres_new =[];
for di=1:4
	for mdli=1:3
		for wli=1:WL_N
			bwl = 1;
			switch wli
				case 1
					%full vector			
					hpc.wrplot = zeros(hpc.Nc, hpc.ipl, ll);
					hpc.wrplot(:,5,:) = 1;
					% Max Freq all time:
					hpc.frplot = 3.45 * ones(min(ceil(hpc.tsim / hpc.Ts_input)+1,(hpc.tsim/hpc.Ts_ctrl+1)), hpc.Nc);
				case 2
					%full idle but 9 cores
					ll = ceil(hpc.tsim*1e6/hpc.quantum_us);
					hpc.wrplot = zeros(hpc.Nc, hpc.ipl, ll);
					hpc.wrplot(:,bwl,:) = 1;

					% Med Freq all time:
					hpc.frplot = hpc.F_min * ones(min(ceil(hpc.tsim / hpc.Ts_input)+1,(hpc.tsim/hpc.Ts_ctrl+1)), hpc.Nc);

					coreid1 = [4 8 12 17 20 22 27 31 36];			
					hpc.wrplot(coreid1,bwl,:) = 0;
					hpc.wrplot(coreid1,5,:) = 1;
					% Max Freq all time:
					hpc.frplot(:, coreid1) = 3.45;
				
					coreid2 = [1 10 11 14 16 24 29 32 34];			
					hpc.wrplot(coreid2,bwl,:) = 0;
					hpc.wrplot(coreid2,3,:) = 0.70;
					hpc.wrplot(coreid2,2,:) = 0.30;
					hpc.frplot(:, coreid2) = 2.7;
					
					%
					%
					%
					coreid12 = [ 1 4 8 10 11 12 14 16 17 20 22 24 27 29 31 32 34 36];
					
					for i=1:3
						f = fres{di, mdli, i, wli};
						fcp = f(2:end,:);					

						smref = round((size(fcp,1))/(size(hpc.frplot(2:end,:),1)));
						smf = round((size(hpc.frplot(2:end,:),1))/(size(fcp,1)));

						fd = repelem(hpc.frplot(2:end,:),max(smref,1),1) - repelem(fcp,max(smf,1),1);
						fd = fd(:, coreid12);

						n2 = reshape(fd,[],1);
						perf_fd2norm(di, mdli, i) = norm(n2);

						w = wlres{di, mdli, i, wli}./ perf_max_check{wli} * 100;

						w = w(coreid12);

						perf_wlMax(di, mdli, i) = max(w);
						perf_wlAv(di, mdli, i) = mean(w);
						perf_wlmin(di, mdli, i) = min(w);
					end				
					
				case 3
					%hpc.wrplot = wlwl3;
					% Max Freq all time:
					hpc.frplot = 3.45 * ones(min(ceil(hpc.tsim / hpc.Ts_input)+1,(hpc.tsim/hpc.Ts_ctrl+1)), hpc.Nc);
				case 4
					%hpc.wrplot = wlwl3;
					% Max Freq all time:
					hpc.frplot = 3.45 * ones(min(ceil(hpc.tsim / hpc.Ts_input)+1,(hpc.tsim/hpc.Ts_ctrl+1)), hpc.Nc);
			end

			for i=1:3
				%[di, mdli, i, wli]
				wlres_new{di, mdli, i, wli} = wlres{di, mdli, i, wli} ./ perf_max_check{wli} * 100;
				tres{di, mdli, i, wli} = hpc.stats_analysis(xres{di, mdli, i, wli}, ures{di, mdli, i, wli}, ...
					fres{di, mdli, i, wli}, vres{di, mdli, i, wli}, wlres_new{di, mdli, i, wli});
			end
		end
	end
end

end %needs2do

%% TEMPERATURE

clc;
%WL = 2;
%MODEL = 2;
%DOM = 4;

%rowj = 0;
BARj = [];

for WL=1:WL_N
	rowj = 0;
	for MODEL=1:3
		for DOM=1:4
			clc;
			rowj = rowj + 1;
			lj = 1;
			%
			disp("Ex Max");
			strj = "";
			for i=1:3
				if isempty(tres{DOM, MODEL, i, WL}.temp.exMaxCr)
					BARj( WL, lj, rowj , i ) = 0;
				else
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.exMaxCr;
				end
				
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.exMaxCr, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			%

			%
			disp("Ex 95");
			strj = "";
			for i=1:3
				%BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.ex95pCr;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.ex95pCr, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			%lj = lj+1;				
			disp(strj);
			%

			%
			disp("Ex Av");
			strj = "";
			for i=1:3	
				BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.exAvCr;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.exAvCr, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			%

			%
			disp("Time Ex tot");
			strj = "";
			for i=1:3
				%BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.exTotTimeCr;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.exTotTimeCr, '%.3f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			%lj = lj+1;
			disp(strj);
			%

			%
			disp("Time Ex Av %");
			strj = "";
			for i=1:3	
				BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.exAvTimeCr/0.75*100;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.exAvTimeCr/0.75*100, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			%

			%
			disp("T Max %");
			strj = "";
			for i=1:3	
				%BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.Max;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.Max, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			%lj = lj+1;
			disp(strj);
			%

			%
			disp("T Av %");
			strj = "";
			for i=1:3	
				BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.Av;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.Av, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			
		end
	end
end

figure();
tt = size(BARj,2);

axgrid = [tt,WL_N];  % [#rows, #cols]
titles = {'Max Exceeded Temperature [°C]', 'Average Exceeded Temperature [°C]', 'Average Exceeded Time [%]', 'Average Temperature [°C]'}; 
xlabelsj = {'W:1D', 'W:4D', 'W:9D', 'W:AD', 'A:1D', 'A:4D', 'A:9D', 'A:AD', 'R:1D', 'R:4D', 'R:9D', 'R:AD' };
wltitles = {'MAX-WL', 'MULTI-WL', 'CLOUD-WL'};
tclMain = tiledlayout(axgrid(1),1); 
tclMain.TileSpacing = 'compact';
tclMain.Padding = 'compact';
tcl = gobjects(1,axgrid(1));
ax = gobjects(axgrid); 
for i=1:tt
	tcl(i) = tiledlayout(tclMain,1,axgrid(2));
	tcl(i).Layout.Tile = i; 
	for WL=1:WL_N
		index = (i-1)*WL_N + WL;
		ax(i,WL) = nexttile(tcl(i));		
		bar(squeeze(BARj(WL, i,:,:)));
		xticklabels(xlabelsj );
		xtickangle(20);
		title(wltitles(WL));				
	end
	title(tcl(i),titles{i});
end
%

%% POWER
clc;
%WL = 1;
%MODEL = 1;
%DOM = 1;

%rowj = 0;
BARj = [];

for WL=1:WL_N
	rowj = 0;
	for MODEL=1:3
		for DOM=1:4
			clc;
			rowj = rowj + 1;
			lj = 1;

			%
			disp("Ex Max W");
			strj = "";
			for i=1:3	
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.power.exMaxW, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			disp(strj);
			%

			%
			disp("Ex Max %");
			strj = "";
			for i=1:3
				if isempty(tres{DOM, MODEL, i, WL}.power.exMaxP)
					BARj( WL, lj, rowj , i ) = 0;
				else
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.power.exMaxP*100;
				end				
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.power.exMaxP*100, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			%

			%
			disp("Ex 95 W");
			strj = "";
			for i=1:3	
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.power.ex95W, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			disp(strj);
			%

			%
			disp("Ex 95 P");
			strj = "";
			for i=1:3
				BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.power.ex95P*100;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.power.ex95P*100, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			%

			%
			disp("Ex Av W");
			strj = "";
			for i=1:3	
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.power.exAvW, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			disp(strj);
			%

			%
			disp("Ex Av P");
			strj = "";
			for i=1:3	
				BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.power.exAvP*100;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.power.exAvP*100, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			%

			%
			disp("Ex Time Tot");
			strj = "";
			for i=1:3	
				BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.power.exTime/1.5*100;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.power.exTime/1.5*100, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			%
		end
	end
end

figure();
tt = size(BARj,2);

axgrid = [tt,WL_N];  % [#rows, #cols]
titles = {'Max Exceeded Power [%]', '95-p Exceeded Power [%]', 'Average Exceeded Power [%]', 'Total Exceeded Time [%]'}; 
xlabelsj = {'W:1D', 'W:4D', 'W:9D', 'W:AD', 'A:1D', 'A:4D', 'A:9D', 'A:AD', 'R:1D', 'R:4D', 'R:9D', 'R:AD' };
wltitles = {'MAX-WL', 'MULTI-WL', 'CLOUD-WL'};
tclMain = tiledlayout(axgrid(1),1); 
tclMain.TileSpacing = 'compact';
tclMain.Padding = 'compact';
tcl = gobjects(1,axgrid(1));
ax = gobjects(axgrid); 
for i=1:tt
	tcl(i) = tiledlayout(tclMain,1,axgrid(2));
	tcl(i).Layout.Tile = i; 
	for WL=1:WL_N
		index = (i-1)*WL_N + WL;
		ax(i,WL) = nexttile(tcl(i));		
		bar(squeeze(BARj(WL, i,:,:)));
		xticklabels(xlabelsj );
		xtickangle(20);
		title(wltitles(WL));				
	end
	title(tcl(i),titles{i});
end


BARj(:,:,:,2) = 0;
figure();
tclMain = tiledlayout(axgrid(1),1); 
tclMain.TileSpacing = 'compact';
tclMain.Padding = 'compact';
tcl = gobjects(1,axgrid(1));
ax = gobjects(axgrid);  
for i=1:tt
	tcl(i) = tiledlayout(tclMain,1,axgrid(2));
	tcl(i).Layout.Tile = i; 
	for WL=1:WL_N
		index = (i-1)*WL_N + WL;
		ax(i,WL) = nexttile(tcl(i));		
		bar(squeeze(BARj(WL, i,:,:)));
		xticklabels(xlabelsj );
		xtickangle(20);
		title(wltitles(WL));				
	end
	title(tcl(i),titles{i});
end

%% Performance
clc;
%WL = 1;
%MODEL = 1;
%DOM = 1;

%rowj = 0;
BARj = [];

for WL=1:WL_N
	rowj = 0;
	for MODEL=1:3
		for DOM=1:4
			clc;
			rowj = rowj + 1;
			lj = 1;

			%
			disp("2-norm");
			strj = "";
			if (WL == 2) && (old_perf)
				for i=1:3
					BARj( WL, lj, rowj , i ) = perf_fd2norm(DOM, MODEL, i) / sqrt(18) / sqrt(4000);
					if i==3
						BARj( WL, lj, rowj , i ) = BARj( WL, lj, rowj , i ) / sqrt(2);
					end
					strj = strcat( strj, num2str(BARj( WL, lj, rowj , i ), '%.2f'));	
					if i<3
						strj = strcat( strj, " | ");
					end
				end
			else
				for i=1:3	
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.perf.fd2norm / sqrt(hpc.Nc) / sqrt(4000);
					if i==3
						BARj( WL, lj, rowj , i ) = BARj( WL, lj, rowj , i ) / sqrt(2);
					end
					strj = strcat( strj, num2str(BARj( WL, lj, rowj , i ), '%.2f'));	
					if i<3
						strj = strcat( strj, " | ");
					end
				end
			end
			
			lj = lj+1;
			disp(strj);
			%
				
			%
			disp("Max perf");
			strj = "";
			if (WL == 2) && (old_perf)
				for i=1:3
					BARj( WL, lj, rowj , i ) = perf_wlMax(DOM, MODEL, i);
					strj = strcat( strj, num2str(BARj( WL, lj, rowj , i ), '%.2f'));	
					if i<3
						strj = strcat( strj, " | ");
					end
				end
			else
				for i=1:3	
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.perf.wlMax;
					strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.perf.wlMax, '%.2f'));	
					if i<3
						strj = strcat( strj, " | ");
					end
				end
			end
			lj = lj+1;
			disp(strj);
			%
			
			%
			disp("average perf");
			strj = "";
			if (WL == 2) && (old_perf)
				for i=1:3
					BARj( WL, lj, rowj , i ) = perf_wlAv(DOM, MODEL, i);
					strj = strcat( strj, num2str(perf_wlAv(DOM, MODEL, i), '%.2f'));	
					if i<3
						strj = strcat( strj, " | ");
					end
				end
			else
				for i=1:3	
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.perf.wlAv;
					strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.perf.wlAv, '%.2f'));	
					if i<3
						strj = strcat( strj, " | ");
					end
				end
			end
			lj = lj+1;
			disp(strj);
			%
			
			%
			disp("min Perf");
			strj = "";
			if (WL == 2) && (old_perf)
				for i=1:3
					BARj( WL, lj, rowj , i ) = perf_wlmin(DOM, MODEL, i);
					strj = strcat( strj, num2str(perf_wlmin(DOM, MODEL, i), '%.2f'));	
					if i<3
						strj = strcat( strj, " | ");
					end
				end
			else
				for i=1:3	
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.perf.wlmin;
					strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.perf.wlmin, '%.2f'));	
					if i<3
						strj = strcat( strj, " | ");
					end
				end
			end
			lj = lj+1;
			disp(strj);
			%
		end
	end
end

figure();
tt = size(BARj,2);

axgrid = [tt,WL_N];  % [#rows, #cols]
titles = {'2-norm of the \Delta F_{pe}', 'Max Executed Workload (MAX-Wl_p) [%]', 'Average Executed Workload (AV-Wl_p) [%]', 'min Executed Workload (MIN-Wl_p) [%]'}; 
xlabelsj = {'W:1D', 'W:4D', 'W:9D', 'W:AD', 'A:1D', 'A:4D', 'A:9D', 'A:AD', 'R:1D', 'R:4D', 'R:9D', 'R:AD' };
wltitles = {'MAX-WL', 'MULTI-WL', 'CLOUD-WL'};
tclMain = tiledlayout(axgrid(1),1); 
tclMain.TileSpacing = 'compact';
tclMain.Padding = 'compact';
tcl = gobjects(1,axgrid(1));
ax = gobjects(axgrid); 
for i=1:tt
	tcl(i) = tiledlayout(tclMain,1,axgrid(2));
	tcl(i).Layout.Tile = i; 
	for WL=1:WL_N
		index = (i-1)*WL_N + WL;
		ax(i,WL) = nexttile(tcl(i));		
		bar(squeeze(BARj(WL, i,:,:)));
		xticklabels(xlabelsj );
		xtickangle(20);
		title(wltitles(WL));	
		if i > 1
			ylim([min(BARj(:, [2:end],:,:),[],"all")-5 min( 100, max(BARj(:, [2:end],:,:),[],"all")+5)]);
		end
	end
	title(tcl(i),titles{i});
end


%%
alg = 1;
qt_comp = 0.5;
maxj = 0;
minj = 10000;
acc1 = 0;
acc2 = 0;
acc3 = 0;
countx = 0;
for d=1:3
	for m=1:3
		countx = countx +1;
		%acc1 = acc1 + tres{d, m, 1, 3}.perf.wlAv;
		%acc2 = acc2 + tres{d, m, 2, 3}.perf.wlAv;
		%acc3 = acc3 + tres{d, m, 3, 3}.perf.wlAv;
		%acc1 = acc1 + perf_wlAv(d, m, 1);
		%acc2 = acc2 + perf_wlAv(d, m, 2);
		%acc3 = acc3 + perf_wlAv(d, m, 3);
		acc1 = acc1 + perft2{d, m, 1, 2}.wlAv;
		acc2 = acc2 + perft2{d, m, 2, 2}.wlAv;
		acc3 = acc3 + perft2{d, m, 3, 2}.wlAv;
		
		%acc1 = acc1 + tres{d, m, 1, 3}.perf.fd2norm / sqrt(hpc.Nc) / sqrt(4000);
		%acc2 = acc2 + tres{d, m, 2, 3}.perf.fd2norm / sqrt(hpc.Nc) / sqrt(4000);
		%acc3 = acc3 + tres{d, m, 3, 3}.perf.fd2norm / sqrt(hpc.Nc) / sqrt(4000) / sqrt(2);
		%acc1 = acc1 + perf_fd2norm(d, m, 1) / sqrt(hpc.Nc) / sqrt(4000);
		%acc2 = acc2 + perf_fd2norm(d, m, 2) / sqrt(hpc.Nc) / sqrt(4000);
		%acc3 = acc3 + perf_fd2norm(d, m, 3) / sqrt(hpc.Nc) / sqrt(4000) / sqrt(2);
		
		%acc = acc + tres{d, m, 3, 2}.power.exTime/1.5*100;
		for w=1:3
			
			%
			quant = tres{d, m, alg, w}.temp.exAvTimeCr/0.75*100;
			%
			if quant > maxj
				maxj = quant;
			end
			if quant < minj
				minj = quant;
			end
			if quant > qt_comp
				%dispstrcat("mod: ", int2str(m), " - wl: ", int2str(w), " - dom: ",int2str(d))
				%disp(quant)
			end
		end
	end
end

disp(strcat("Max: ", num2str(maxj)));
disp(strcat("min: ", num2str(minj)));

acc1 / countx
acc2 / countx
acc3 / countx


%%
clc;

perf_avg = squeeze(BARj( 3, 3, : , : ))'

diff_F_V = perf_avg(1, :) - perf_avg(3,:)
sumd_av = [];
sumd2_av = [];
for i=1:3
	for j=1:3
		inx = (i-1)*4+1;
		sumd_av(j,i) = std(perf_avg(j, [inx:inx+3 ]));
	end
end
for j=1:3
	sumd2_av(j) = std(perf_avg(j, :));
end
sumd_av
sumd2_av = sumd2_av'

perf_min = squeeze(BARj( 1, 4, : , : ))'
sumd_min = [];
sumd2_min = [];
normavmin = [];
for i=1:3
	for j=1:3
		inx = (i-1)*4+1;
		sumd_min(j,i) = std(perf_min(j, [inx:inx+3 ]));
	end
end
for j=1:3
	sumd2_min(j) = std(perf_min(j, :));
end
sumd_min
sumd2_min = sumd2_min'

for j=1:3
	normavmin(j) = norm(perf_avg(j, :) - perf_min(j, :));
end

normavmin = normavmin'

perf_avg - perf_min



%% Performance 2

perft2 = [];
wlres_new = [];
for di=1:4
	for mdli=1:3
		wli=2;
		bwl = 1;

		%full idle but 9 cores
		ll = ceil(hpc.tsim*1e6/hpc.quantum_us);
		hpc.wrplot = zeros(hpc.Nc, hpc.ipl, ll);
		hpc.wrplot(:,bwl,:) = 1;

		% Med Freq all time:
		hpc.frplot = hpc.F_min * ones(min(ceil(hpc.tsim / hpc.Ts_input)+1,(hpc.tsim/hpc.Ts_ctrl+1)), hpc.Nc);

		coreid1 = [4 8 12 17 20 22 27 31 36];			
		hpc.wrplot(coreid1,bwl,:) = 0;
		hpc.wrplot(coreid1,5,:) = 1;
		% Max Freq all time:
		hpc.frplot(:, coreid1) = 3.45;

		coreid2 = [1 10 11 14 16 24 29 32 34];			
		hpc.wrplot(coreid2,bwl,:) = 0;
		hpc.wrplot(coreid2,3,:) = 0.70;
		hpc.wrplot(coreid2,2,:) = 0.30;
		hpc.frplot(:, coreid2) = 2.7;

		%
		%
		%
		coreid12 = [ 1 4 8 10 11 12 14 16 17 20 22 24 27 29 31 32 34 36];
		coreid3 = [2 3 5 6 7 9 13 15 18 19 21 23 25 26 28 30 33 35];
		

		for i=1:3
			for cc=1:2
				fcp = fres{di, mdli, i, wli};
				fcp = fcp(2:end,:);
				wlres_new{di, mdli, i, wli} = wlres{di, mdli, i, wli} ./ perf_max_check{wli} * 100;
				w = wlres_new{di, mdli, i, wli};

				smref = round((size(fcp,1))/(size(hpc.frplot(2:end,:),1)));
				smf = round((size(hpc.frplot(2:end,:),1))/(size(fcp,1)));

				fd = repelem(hpc.frplot(2:end,:),max(smref,1),1) - repelem(fcp,max(smf,1),1);
				
				switch cc
					case 1
						fd = fd(:,coreid1);
						w = w(coreid1);
					case 2
						fd = fd(:,coreid2);
						w = w(coreid2);
					case 3
						fd = fd(:,coreid3);
						w = w(coreid3);
				end

				n2 = reshape(fd,[],1);
				perft2{di, mdli, i, cc}.fd2norm = norm(n2);
				%perf.fdAv2norm = mean(norm
				perft2{di, mdli, i, cc}.fdMax = max(fd(:));
				perft2{di, mdli, i, cc}.fdmin = min(fd(:));
				perft2{di, mdli, i, cc}.fdAv = mean(fd(:));
				perft2{di, mdli, i, cc}.fdSum = sum(fd(:));
				perft2{di, mdli, i, cc}.fdAvSD = std(mean(fd));
				perft2{di, mdli, i, cc}.fdSDMean = mean(std(fd));


				perft2{di, mdli, i, cc}.wlMax = max(w);
				perft2{di, mdli, i, cc}.wlAv = mean(w);
				perft2{di, mdli, i, cc}.wlmin = min(w);
				perft2{di, mdli, i, cc}.wlSD = std(w);
			end
		end
				
	end
end

clc;

BARj = [];

for cc=1:2
	rowj = 0;
	for MODEL=1:3
		for DOM=1:4
			clc;
			rowj = rowj + 1;
			lj = 1;

			%
			disp("2-norm");
			strj = "";
			disp(strcat(num2str(di), " - ", num2str(mdli), " - ", num2str(i), " - ", num2str(cc)));
			for i=1:3	
				BARj( cc, lj, rowj , i ) = perft2{DOM, MODEL, i, cc}.fd2norm / sqrt(hpc.Nc) / sqrt(4000);
				if i==3
					BARj( cc, lj, rowj , i ) = BARj( cc, lj, rowj , i ) / sqrt(2);
				end
				strj = strcat( strj, num2str(BARj( cc, lj, rowj , i ), '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end			
			lj = lj+1;
			disp(strj);
			%
				
			%
			disp("Max perf");
			strj = "";
			for i=1:3	
				BARj( cc, lj, rowj , i ) = perft2{DOM, MODEL, i, cc}.wlMax;
				strj = strcat( strj, num2str(perft2{DOM, MODEL, i, cc}.wlMax, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			%
			
			%
			disp("average perf");
			strj = "";
			for i=1:3	
				BARj( cc, lj, rowj , i ) = perft2{DOM, MODEL, i, cc}.wlAv;
				strj = strcat( strj, num2str(perft2{DOM, MODEL, i, cc}.wlAv, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			%
			
			%
			disp("min Perf");
			strj = "";
			for i=1:3	
				BARj( cc, lj, rowj , i ) = perft2{DOM, MODEL, i, cc}.wlmin;
				strj = strcat( strj, num2str(perft2{DOM, MODEL, i, cc}.wlmin, '%.2f'));	
				if i<3
					strj = strcat( strj, " | ");
				end
			end
			lj = lj+1;
			disp(strj);
			%
		end
	end
end

%%
figure();
tt = size(BARj,2);

axgrid = [tt,2];  % [#rows, #cols]
titles = {'2-norm of the \Delta F_{pe}', 'Max Executed Workload [%]', 'Average Executed Workload [%]', 'min Executed Workload [%]'}; 
xlabelsj = {'W:1D', 'W:4D', 'W:9D', 'W:AD', 'A:1D', 'A:4D', 'A:9D', 'A:AD', 'R:1D', 'R:4D', 'R:9D', 'R:AD' };
wltitles = {'VECT WL @3.45GHz', 'INT/FLOAT WL @2.7GHz', 'CLOUD-WL'};
tclMain = tiledlayout(axgrid(1),1); 
%tclMain.TileSpacing = 'compact';
%tclMain.Padding = 'compact';
tcl = gobjects(1,axgrid(1));
ax = gobjects(axgrid); 
for i=1:tt
	tcl(i) = tiledlayout(tclMain,1,axgrid(2));
	tcl(i).Layout.Tile = i; 
	for cc=1:2
		index = (i-1)*2 + cc;
		ax(i,cc) = nexttile(tcl(i));		
		bar(squeeze(BARj(cc, i,:,:)));
		xticklabels(xlabelsj );
		xtickangle(20);
		title(wltitles(cc));	
		if i > 1
			ylim([min(BARj(:, [2:end],:,:),[],"all")-5 100]);
		end
	end
	title(tcl(i),titles{i});
end


%%



