
%% 1: Instantiate the class
hpc = hpc_lab;

path = "/home/eventine/Documents/MATLAB/Results";
%path = "C:\Users\giovanni.bambini2\OneDrive - Alma Mater Studiorum Università di Bologna\___PhD\Papers\TAAS-2023\Results\Automated";

if isunix
	separator_os = "/";
else
	separator_os = "\";
end
				
%% 2: Define System Parameters: Init + Thermal

% Cores:
hpc.Nc = 36; %72;
hpc.Nh = 6; %9;			% Number of rows
hpc.Nv = 6; %8;			% Number of cols

% Version of the Thermal Model
hpc.model_ver = 0;
% Exponential Leakage
hpc.leak_exp = 1;
% Noise and Variation
hpc.pw_dev_per = 1;
hpc.param_dev_per = 1;
hpc = hpc.create_model_deviation();

hpc.sensor_noise = 1;

% External Ambient Temperature
hpc.t_outside = 25.0 + 273.15;


%% 3: Simulation setup

% SIMULATION TIME (in sec)
hpc.tsim = 2.0;
% Simulation Time-Step
hpc.Ts = 50e-6;
% Workload Quantum-step
hpc.quantum_us = 50;

% Input Time-Step
hpc.Ts_target = 1e-3;

% Create Thermal Model:
hpc.default_floorplan_config();
hpc.create_model_deviation();
hpc.model_init();

% Probability for each WL	
hpc.wl_prob = [0 2 3 2 2.5];	
% minimal execution in us
hpc.wl_min_exec_us = [50e3/8 5e3 5e3 5e3 100e3];
% mean execution in us
hpc.wl_mean_exec_us = [150e3/8 10e3 10e3 10e3 200e3];	

hpc.sensor_noise_amplitude = [1.0; 1.0; 1.0]; %og was 0.5
hpc.F_discretization_step = 0.05;

hpc.pw_gmean = 0.02;
hpc.pw_gvar = 1;
hpc.pw_glim = [-0.02 0.05];
hpc.pw_3sigma_on3 = 0.05;

hpc.create_core_pw_noise();

%% 4 Controller

% Max Temperature
hpc.core_limit_temp = 85 + 273.15;
hpc.core_crit_temp = 95 + 273.15;
% Max Power
tt = min(ceil(hpc.tsim / hpc.Ts_target)+1, (hpc.tsim/1e-4+1));
%Target Power Budget
hpc.tot_pw_budget = 450/36*hpc.Nc*ones(tt,1);
tu = ceil(tt/4);
hpc.tot_pw_budget(tu+1:2*tu) = 2*hpc.Nc;
hpc.tot_pw_budget(2*tu+1:3*tu) = 5*hpc.Nc;
hpc.tot_pw_budget(3*tu+1:3*tu+ceil(tu/2)) = 3*hpc.Nc;
hpc.tot_pw_budget(3*tu+ceil(tu/2)+1:end) = 8*hpc.Nc;

% Delays
hpc.delay_F_mean = 1e-5;		% Mean Frequency application Delay
hpc.delay_V_mean = 1e-5;		% Mean Voltage application Delay

%{
%hpc.wltrc = zeros(size(hpc.wltrc));
%hpc.wltrc(:,hpc.ipl,:) = 1;
%hpc.wltrc(:,2,:) = 1;
hpc.alpha_wl = 0.08; %0.4;

hpc.kp = 1.9518;				% PID Kp
hpc.ki = 73.931;				% PID Ki
%increasing this, smoothen the F, but violate the Temp cap
hpc.aw_up = 2; %0.05;				% PID Anti-Windup Upper Saturation
%this doesn't change if I increase the absolute value
%but decreasing this change things a lot!
%putting it to 0.1 is smoother (only in exp_leakage=0 case)
%putting it to 0.01 nothing works!
%this is the cap on how much the integrator can reduce it. 
% -1 = -pu (+ -proportional!!)
hpc.aw_down_c = -0.5; %-0.75;			% PID Anti-Windup Down Sat. Coefficient
%increasing this >0 does nothing
%decreasing this will more noise and a bit more ""smoothe"", but also increase performance
hpc.sat_up = -0.5; %0;				% PID Upper Saturation

%These 3 actually don't change much the result
hpc.pid_e_down = 2.5;			
hpc.pid_e_up = 1;
hpc.pid_e_band_coeff = 0.5;

hpc.aw_up = 0.05;
hpc.aw_down_c = -0.75;
hpc.sat_up = 0;
%}

%% 5: Init & Check response:

%hpc = hpc.create_thermal_model();
%hpc.thermal_model_ver = 1;
%hpc.exp_leakage = 0;
hpc.t_init = (25+273.15)*ones(hpc.Ns,1); %temp_amb*ones(hpc.Ns,1);
%hpc.x_init(end-1:end,1) = temp_amb+20;
%hpc.tesim = 1.5;
%hpc.simulate_aut();

%hpc.exp_leakage = 1;
%hpc.simulate_aut();

%Others, TODO
hpc.min_pw_red = 0.6;

%% ITERATE ON WORKLOAD
wl_times = 3;
th_models = 3; %2;
ndom = 4;
robust = 0;
show = 0;
savetofile = 0;
hpc.compare_vs_baseline = 1;
hpc.sensor_noise = 0;

test_iter = 10;%20;

tres = [];
tres{test_iter, ndom, th_models, 3, wl_times} = [];

xres = [];
%xres{ndom, th_models, 3, wl_times} = [];
ures = [];
%ures{ndom, th_models, 3, wl_times} = [];
fres = [];
%fres{ndom, th_models, 3, wl_times} = [];
vres = [];
%vres{ndom, th_models, 3, wl_times} = [];
wlres = [];
wlres{test_iter, ndom, th_models, 3, wl_times} = [];

if hpc.Nc == 36
	coreid1 = [4 8 12 17 20 22 27 31 36];
	coreid2 = [1 10 11 14 16 24 29 32 34];
	coreid12 = [ 1 4 8 10 11 12 14 16 17 20 22 24 27 29 31 32 34 36];
else
	error("didn't set coreid up!");
end

%{
p1 = 94.2399*ones(hpc.Nc,1);
p2 = 10.9303*ones(hpc.Nc,1);
p2(coreid1) = 94.2399;
p2(coreid2) = 68.7267;

p3 = [   91.2498; 88.9647; 87.6472; 90.7673; 88.1172; 89.4897; 88.7772; 89.0797; 88.9622; 89.5697;
		 90.2273; 88.7547; 88.6847; 88.2572; 90.0673; 88.8922; 88.5272; 88.0022; 89.2972; 87.3722;
		 88.4622; 89.3922; 88.1097; 89.5147; 88.3922; 88.8197; 87.2572; 88.7122; 90.2573; 87.3047;
		 89.3597; 86.9047; 87.7922; 89.0322; 88.6722; 90.6123
   ];

perf_max_check{1} = p1;
perf_max_check{2} = p2;
perf_max_check{3} = p3;
%}


%% Last modifications
% ?????
%{
hpc.measure_noise = 0; 

hpc.kp = 1.9518*0.3;				% PID Kp
hpc.ki = 73.931/1.7;				% PID Ki
hpc.aw_up = 0.05; %0.05;				% PID Anti-Windup Upper Saturation
hpc.pid_e_down = 5;			
hpc.pid_e_up = 0.5;
hpc.pid_e_band_coeff = 0.3;

hpc.ctrl_MA = 1;
hpc.iterative_fv = ~hpc.ctrl_MA;
%}

max_amb_T = (45+273.15);
max_init_T = (85 + 273.15) - 25;
max_elem_T = 55+273.15;

test_T_step = (max_init_T - hpc.t_outside) /  (test_iter-1);
err_ampl = 0.5;
rand_T_init = rand(hpc.Ns,test_iter)*err_ampl*2 - err_ampl*ones(hpc.Ns,test_iter);
init_cond = ones(hpc.Ns,1) * (hpc.t_outside:test_T_step:max_init_T);

for i=-1:0
    tr = init_cond(end+i,:) > max_amb_T;
    init_cond(end+i,tr) = max_amb_T;
end

for i=-3:-2
    tr = init_cond(end+i,:) > max_elem_T;
    init_cond(end+i,tr) = max_elem_T;
end

init_cond = init_cond + rand_T_init;

%% 
% Import the file
cpth = pwd;
fileToRead1 = strcat(cpth, separator_os, 'wl_journal_good.mat');
newData1 = load('-mat', fileToRead1);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
	assignin('base', vars{i}, newData1.(vars{i}));
end

% Import the file
fileToRead1 = strcat(cpth, separator_os, 'wl4.mat');
newData1 = load('-mat', fileToRead1);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
	assignin('base', vars{i}, newData1.(vars{i}));
end

% Import the file
fileToRead1 = strcat(cpth, separator_os, 'TAAS_Init_cond.mat');
newData1 = load('-mat', fileToRead1);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
	assignin('base', vars{i}, newData1.(vars{i}));
end


figid = 1;

addpath('Controllers/');

%%
for tit=1:test_iter

	hpc.t_init = init_cond(:,tit);

for wli=1:wl_times
	
	tic
		
	bwl = 1;
    tt = min(ceil(hpc.tsim / hpc.Ts_target)+1, (hpc.tsim/1e-4+1));
	switch wli
		case 1
			%full vector
			ll = ceil(hpc.tsim*1e6/hpc.quantum_us);
			hpc.wltrc = zeros(hpc.Nc, hpc.ipl, ll);
			hpc.wltrc(:,5,:) = 1;
			% Max Freq all time:
			hpc.frtrc = 3.45 * ones(tt,hpc.Nc);
		case 2
			%full idle but 9 cores
			ll = ceil(hpc.tsim*1e6/hpc.quantum_us);
			hpc.wltrc = zeros(hpc.Nc, hpc.ipl, ll);
			hpc.wltrc(:,bwl,:) = 1;
			
			% Med Freq all time:
			hpc.frtrc = hpc.F_min * ones(tt,hpc.Nc);
			
			hpc.wltrc(coreid1,bwl,:) = 0;
			hpc.wltrc(coreid1,5,:) = 1;
			% Max Freq all time:
			hpc.frtrc(:, coreid1) = 3.45;
							
			hpc.wltrc(coreid2,bwl,:) = 0;
			hpc.wltrc(coreid2,3,:) = 0.70;
			hpc.wltrc(coreid2,2,:) = 0.30;
			hpc.frtrc(:, coreid2) = 2.7;
		case 3
			hpc.wltrc = wlwl3;
			% Max Freq all time:
			hpc.frtrc = 3.45 * ones(tt,hpc.Nc);
		case 4
			hpc.wltrc = wlwl4;
			% Max Freq all time:
			hpc.frtrc = 3.45 * ones(tt,hpc.Nc);
	end

	if hpc.compare_vs_baseline
		% Compute ideal&unrestricted baseline
		%since atm it is domain independent:
			hpc.vd = hpc.Nc;
			hpc.VDom = eye(hpc.Nc);
		[wmaxc, perfmaxc] = hpc.base_ideal_unr();
		awl = hpc.perf_max_check;
		perf_max_check{wli} = awl;
	end
	
	%% ITERATE ON Thermal Models
	for mdli=1:th_models
	
		switch mdli
			case 1
				hpc.model_ver = 0;
				hpc.model_init();
			case 2
				hpc.model_ver = 1;
				hpc.model_init();
			case 3
				hpc.model_ver = 2;
				hpc.model_init();
		end
		
		%% ITERATE ON DOMAINS
		for di=1:ndom
			% cases:
			% 1: Single Domain
			% 2: 4 Domains
			% 3: 9 Domains
			% 4: All Domains (1 per core)

			switch di
				case 1
					hpc.vd = 1;
					% Structure of the Voltage Domains Nc x vd
					hpc.VDom = ones(hpc.Nc, 1);
				case 2
					hpc.vd = 4;
					% Structure of the Voltage Domains Nc x vd
					hpc.VDom = zeros(hpc.Nc, hpc.vd);
					if hpc.Nc == 36
						hpc.VDom([1:3 7:9 13:15], 1) = 1;
						hpc.VDom([4:6 10:12 16:18], 2) = 1;
						hpc.VDom([19:21 25:27 31:33], 3) = 1;
						hpc.VDom([22:24 28:30 34:36], 4) = 1;
					elseif hpc.Nc == 72
						hpc.VDom([1:6 13:18 25:30], 1) = 1;
						hpc.VDom([7:12 19:24 31:36], 2) = 1;
						hpc.VDom([37:42 49:54 61:66], 3) = 1;
						hpc.VDom([43:48 55:60 67:72], 4) = 1;
					else
						error("didn't set VDom for this Nc!");
					end
				case 3
					hpc.vd = 9;
					% Structure of the Voltage Domains Nc x vd
					hpc.VDom = zeros(hpc.Nc, hpc.vd);
					if hpc.Nc == 36
						hpc.VDom([1:2 7:8], 1) = 1;
						hpc.VDom([3:4 9:10], 2) = 1;
						hpc.VDom([5:6 11:12], 3) = 1;
						hpc.VDom([13:14 19:20], 4) = 1;
						hpc.VDom([15:16 21:22], 5) = 1;
						hpc.VDom([17:18 23:24], 6) = 1;
						hpc.VDom([25:26 31:32], 7) = 1;
						hpc.VDom([27:28 33:34], 8) = 1;
						hpc.VDom([29:30 35:36], 9) = 1;
					elseif hpc.Nc == 72
						hpc.VDom([1:4 9:12], 1) = 1;
						hpc.VDom([5:8 13:16], 2) = 1;
						hpc.VDom([17:20 25:28], 3) = 1;
						hpc.VDom([21:24 29:32], 4) = 1;
						hpc.VDom([33:36 41:44], 5) = 1;
						hpc.VDom([37:40 45:48], 6) = 1;
						hpc.VDom([49:52 57:60], 7) = 1;
						hpc.VDom([53:56 61:64], 8) = 1;
						hpc.VDom([65:72], 9) = 1;
					else
						error("didn't set VDom for this Nc!");
					end
				case 4
					hpc.vd = hpc.Nc;
					% Structure of the Voltage Domains Nc x vd
					hpc.VDom = eye(hpc.Nc);
				otherwise
					fprintf("Error! Wrong case selection!\n\r");
					hpc.vd = 1;
					% Structure of the Voltage Domains Nc x vd
					hpc.VDom = ones(hpc.Nc, 1);	
			end

			hpc.quad_pw_budget = 450/36*hpc.Nc*ones(2,hpc.vd);
			itname = strcat('Model: ', int2str(mdli), ' - Domains: ', int2str(di), ' - WL: ', int2str(wli), ' - test: ', int2str(tit));

			%%

			%hpc.dummy_pw = 1;

			% New ControlPulp
			ctrl = Fuzzy;
			ctrl.C = hpc.Cc;
			ctrl.Ts_ctrl = 500e-6;
			
			[xop, uop, fop, vop, wlop] = hpc.simulation(ctrl, show);
			%xres{di, mdli, 1, wli} = xop;
			%ures{di, mdli, 1, wli} = uop; 
			%fres{di, mdli, 1, wli} = fop;
			%vres{di, mdli, 1, wli} = vop;
			wlres{tit, di, mdli, 1, wli} = wlop;
			if wli == 2
				frtrc_hold = hpc.frtrc;
				fop_hold = fop;
				fop = fop(:,coreid12);
				wlop = wlop(coreid12);
				hpc.frtrc = hpc.frtrc(:,coreid12);
				hpc.taas_fix(awl);
			end
			tres{tit, di, mdli, 1, wli} = hpc.stats_analysis(ctrl, xop, uop, fop, vop, wlop);
			if wli == 2
				hpc.frtrc = frtrc_hold;
				fop = fop_hold;
				wlop = wlres{tit,di, mdli, 1, wli};
				hpc.taas_fix(awl);
			end
			
			if show
				figarray = flip(findobj('Type','figure'));
				%TODO Improve: https://it.mathworks.com/help/matlab/ref/arrayfun.html
				for fn=figid:length(figarray)
					figarray(fn).Name = strcat(figarray(fn).Name, ' Alg: NCP - ', itname);
				end
				figid = length(figarray)+1;
			end
			if savetofile
				hpc.saveall(xop, uop, fop, vop, wlop, path, regexprep(strcat('Alg: NCP - ', itname), ':', '_'));
			else
				hpc.savetofile(hpc.xutplot(xop, uop),strcat(path,separator_os,regexprep(strcat('Alg: NCP - ', itname), "-", "TP")),1);
			end
			
			%hpc.paperTPplot(hpc.Ts*[0:4000*10]', xop, uop);
			
			% Old ControlPulp
			ctrl = CP;
			ctrl.C = hpc.Cc;
			ctrl.Ts_ctrl = 500e-6;
			[xop, uop, fop, vop, wlop] = hpc.simulation(ctrl, show);
			%xres{di, mdli, 2, wli} = xop;
			%ures{di, mdli, 2, wli} = uop; 
			%fres{di, mdli, 2, wli} = fop;
			%vres{di, mdli, 2, wli} = vop;
			wlres{tit, di, mdli, 2, wli} = wlop;			
			if wli == 2
				frtrc_hold = hpc.frtrc;
				fop_hold = fop;
				fop = fop(:,coreid12);
				wlop = wlop(coreid12);
				hpc.frtrc = hpc.frtrc(:,coreid12);
				hpc.taas_fix(awl);
			end
			tres{tit, di, mdli, 2, wli} = hpc.stats_analysis(ctrl, xop, uop, fop, vop, wlop);
			if wli == 2
				hpc.frtrc = frtrc_hold;
				fop = fop_hold;
				wlop = wlres{tit,di, mdli, 2, wli};
				hpc.taas_fix(awl);
			end
			
			if show
				figarray = flip(findobj('Type','figure'));
				%TODO Improve: https://it.mathworks.com/help/matlab/ref/arrayfun.html

				for fn=figid:length(figarray)
					figarray(fn).Name = strcat(figarray(fn).Name, ' Alg: CP - ', itname);
				end
				figid = length(figarray)+1;
			end
			
			if savetofile
				hpc.saveall(xop, uop, fop, vop, wlop, path, regexprep(strcat('Alg: CP - ', itname), ':', '_'));
			end

			
			% IBM
			ctrl = IBM_OCC;
			ctrl.C = hpc.Cc;
			ctrl.Ts_ctrl = 250e-6;
			[xop, uop, fop, vop, wlop] = hpc.simulation(ctrl, show);
			%xres{di, mdli, 3, wli} = xop;
			%ures{di, mdli, 3, wli} = uop; 
			%fres{di, mdli, 3, wli} = fop;
			%vres{di, mdli, 3, wli} = vop;
			wlres{tit, di, mdli, 3, wli} = wlop;
			if wli == 2
				frtrc_hold = hpc.frtrc;
				fop_hold = fop;
				fop = fop(:,coreid12);
				wlop = wlop(coreid12);
				hpc.frtrc = hpc.frtrc(:,coreid12);
				hpc.taas_fix(awl);
			end
			tres{tit, di, mdli, 3, wli} = hpc.stats_analysis(ctrl, xop, uop, fop, vop, wlop);
			if wli == 2
				hpc.frtrc = frtrc_hold;
				fop = fop_hold;
				wlop = wlres{tit,di, mdli, 3, wli};
				hpc.taas_fix(awl);
			end

			if show
				figarray = flip(findobj('Type','figure'));
				%TODO Improve: https://it.mathworks.com/help/matlab/ref/arrayfun.html

				for fn=figid:length(figarray)
					figarray(fn).Name = strcat(figarray(fn).Name, ' Alg: IBM OCC - ', itname);
				end
				figid = length(figarray)+1;
			end
			if savetofile
				hpc.saveall(xop, uop, fop, vop, wlop, path, regexprep(strcat('Alg: OCC - ', itname), {':', ' '}, {'_', ''}));
			end
			

		end % for di domains
	end % for mdli thermal model
	toc
	
end % for wl_times

end %test_iter



