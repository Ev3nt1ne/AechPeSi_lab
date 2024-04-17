%%

WL_N = 3;
old_perf = 1;


%%

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
				if isempty(tres{DOM, MODEL, i, WL}.temp.exMn.Max)
					BARj( WL, lj, rowj , i ) = 0;
				else
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.exMn.Max;
				end
				
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.exMn.Max, '%.2f'));	
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
				%BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.exMn.Mean95;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.exMn.Mean95, '%.2f'));	
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
				BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.exMn.MeanAv;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.exMn.MeanAv, '%.2f'));	
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
				%BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.exMn.TotTime;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.exMn.TotTime, '%.3f'));	
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
				BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.temp.exMn.MeanTime/0.75*100;
				strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.temp.exMn.MeanTime/0.75*100, '%.2f'));	
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
tp = size(BARj,2);

axgrid = [tp,WL_N];  % [#rows, #cols]
titles = {'Max Exceeded Temperature [°C]', 'Average Exceeded Temperature [°C]', 'Average Exceeded Time [%]', 'Average Temperature [°C]'}; 
xlabelsj = {'W:1D', 'W:4D', 'W:9D', 'W:AD', 'A:1D', 'A:4D', 'A:9D', 'A:AD', 'R:1D', 'R:4D', 'R:9D', 'R:AD' };
wltitles = {'MAX-WL', 'MULTI-WL', 'CLOUD-WL'};
tclMain = tiledlayout(axgrid(1),1); 
tclMain.TileSpacing = 'compact';
tclMain.Padding = 'compact';
tcl = gobjects(1,axgrid(1));
ax = gobjects(axgrid); 
for i=1:tp
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
tp = size(BARj,2);

axgrid = [tp,WL_N];  % [#rows, #cols]
titles = {'Max Exceeded Power [%]', '95-p Exceeded Power [%]', 'Average Exceeded Power [%]', 'Total Exceeded Time [%]'}; 
xlabelsj = {'W:1D', 'W:4D', 'W:9D', 'W:AD', 'A:1D', 'A:4D', 'A:9D', 'A:AD', 'R:1D', 'R:4D', 'R:9D', 'R:AD' };
wltitles = {'MAX-WL', 'MULTI-WL', 'CLOUD-WL'};
tclMain = tiledlayout(axgrid(1),1); 
tclMain.TileSpacing = 'compact';
tclMain.Padding = 'compact';
tcl = gobjects(1,axgrid(1));
ax = gobjects(axgrid); 
for i=1:tp
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
for i=1:tp
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
				for i=1:3	
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.perf.fd.l2norm;
					strj = strcat( strj, num2str(BARj( WL, lj, rowj , i ), '%.2f'));	
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
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.perf.wl.Max;
					strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.perf.wl.Max, '%.2f'));	
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
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.perf.wl.Av;
					strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.perf.wl.Av, '%.2f'));	
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
					BARj( WL, lj, rowj , i ) = tres{DOM, MODEL, i, WL}.perf.wl.min;
					strj = strcat( strj, num2str(tres{DOM, MODEL, i, WL}.perf.wl.min, '%.2f'));	
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
tp = size(BARj,2);

axgrid = [tp,WL_N];  % [#rows, #cols]
titles = {'2-norm of the \Delta F_{pe}', 'Max Executed Workload (MAX-Wl_p) [%]', 'Average Executed Workload (AV-Wl_p) [%]', 'min Executed Workload (MIN-Wl_p) [%]'}; 
xlabelsj = {'W:1D', 'W:4D', 'W:9D', 'W:AD', 'A:1D', 'A:4D', 'A:9D', 'A:AD', 'R:1D', 'R:4D', 'R:9D', 'R:AD' };
wltitles = {'MAX-WL', 'MULTI-WL', 'CLOUD-WL'};
tclMain = tiledlayout(axgrid(1),1); 
tclMain.TileSpacing = 'compact';
tclMain.Padding = 'compact';
tcl = gobjects(1,axgrid(1));
ax = gobjects(axgrid); 
for i=1:tp
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
		acc1 = acc1 + tres{d, m, 1, 2}.perf.wl.Av;
		acc2 = acc2 + tres{d, m, 2, 2}.perf.wl.Av;
		acc3 = acc3 + tres{d, m, 3, 2}.perf.wl.Av;
		
		%acc1 = acc1 + tres{d, m, 1, 3}.perf.fd2norm / sqrt(hpc.Nc) / sqrt(4000);
		%acc2 = acc2 + tres{d, m, 2, 3}.perf.fd2norm / sqrt(hpc.Nc) / sqrt(4000);
		%acc3 = acc3 + tres{d, m, 3, 3}.perf.fd2norm / sqrt(hpc.Nc) / sqrt(4000) / sqrt(2);
		%acc1 = acc1 + perf_fd2norm(d, m, 1) / sqrt(hpc.Nc) / sqrt(4000);
		%acc2 = acc2 + perf_fd2norm(d, m, 2) / sqrt(hpc.Nc) / sqrt(4000);
		%acc3 = acc3 + perf_fd2norm(d, m, 3) / sqrt(hpc.Nc) / sqrt(4000) / sqrt(2);
		
		%acc = acc + tres{d, m, 3, 2}.power.exTime/1.5*100;
		for w=1:3
			
			%
			quant = tres{d, m, alg, w}.temp.exMn.MeanTime/0.75*100;
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
hpc.frtrc;


%% Performance 2

perft2 = [];
wlres_new = [];
for di=1:4
	for mdli=1:3
		wli=2;
		bwl = 1;

		%full idle but 9 cores
		ll = ceil(hpc.tsim*1e6/hpc.quantum_us);
		hpc.wltrc = zeros(hpc.Nc, hpc.ipl, ll);
		hpc.wltrc(:,bwl,:) = 1;
			
		% Med Freq all time:
		tt = min(ceil(hpc.tsim / hpc.Ts_target)+1, (hpc.tsim/1e-4+1));
		hpc.frtrc = hpc.F_min * ones(tt,hpc.Nc);
			
		hpc.wltrc(coreid1,bwl,:) = 0;
		hpc.wltrc(coreid1,5,:) = 1;
		% Max Freq all time:
		hpc.frtrc(:, coreid1) = 3.45;
							
		hpc.wltrc(coreid2,bwl,:) = 0;
		hpc.wltrc(coreid2,3,:) = 0.70;
		hpc.wltrc(coreid2,2,:) = 0.30;
		hpc.frtrc(:, coreid2) = 2.7;

		%
		%
		%
		coreid12 = [ 1 4 8 10 11 12 14 16 17 20 22 24 27 29 31 32 34 36];
		coreid3 = [2 3 5 6 7 9 13 15 18 19 21 23 25 26 28 30 33 35];		

		for i=1:3
			for cc=1:2
				fcp = fres{di, mdli, i, wli};
				fcp = fcp(2:end,:);
				%wlres_new{di, mdli, i, wli} = wlres{di, mdli, i, wli} ./ perf_max_check{wli} * 100;
				w = wlres{di, mdli, i, wli};

				smref = round((size(fcp,1))/(size(hpc.frtrc(2:end,:),1)));
				smf = round((size(hpc.frtrc(2:end,:),1))/(size(fcp,1)));
				
				fd = repelem(hpc.frtrc(2:end,:),max(smref,1),1) - repelem(fcp,max(smf,1),1);
				
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
				BARj( cc, lj, rowj , i ) = perft2{DOM, MODEL, i, cc}.fd2norm;
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
tp = size(BARj,2);

axgrid = [tp,2];  % [#rows, #cols]
titles = {'2-norm of the \Delta F_{pe}', 'Max Executed Workload [%]', 'Average Executed Workload [%]', 'min Executed Workload [%]'}; 
xlabelsj = {'W:1D', 'W:4D', 'W:9D', 'W:AD', 'A:1D', 'A:4D', 'A:9D', 'A:AD', 'R:1D', 'R:4D', 'R:9D', 'R:AD' };
wltitles = {'VECT WL @3.45GHz', 'INT/FLOAT WL @2.7GHz', 'CLOUD-WL'};
tclMain = tiledlayout(axgrid(1),1); 
%tclMain.TileSpacing = 'compact';
%tclMain.Padding = 'compact';
tcl = gobjects(1,axgrid(1));
ax = gobjects(axgrid); 
for i=1:tp
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

%{
main = fieldnames(tres{1,1,1,1,1})
for i=1:length(main)
	
end


for ctl=1:1
	for WL=1:WL_N
		for MODEL=1:3
			for DOM=1:4
				%tres{2, DOM, MODEL, ctl, WL}.power
				%tres{2, DOM, MODEL, ctl, WL}.temp
				%tres{2, DOM, MODEL, ctl, WL}.freq
				%tres{2, DOM, MODEL, ctl, WL}.vdd

				%
				tres{2, DOM, MODEL, ctl, WL}.perf.fd.Max - tres{1, DOM, MODEL, ctl, WL}.perf.fd.Max
				tres{2, DOM, MODEL, ctl, WL}.perf.fd.min - tres{1, DOM, MODEL, ctl, WL}.perf.fd.min
				tres{2, DOM, MODEL, ctl, WL}.perf.fd.Av - tres{1, DOM, MODEL, ctl, WL}.perf.fd.Av
				tres{2, DOM, MODEL, ctl, WL}.perf.fd.StdAv - tres{1, DOM, MODEL, ctl, WL}.perf.fd.StdAv
			end
		end
	end
end
%}
