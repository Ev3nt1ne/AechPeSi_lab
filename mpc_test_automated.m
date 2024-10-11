
%% 1: Instantiate the class
hpc = hpc_lab;
chip1 = hpc_chiplet;

path = "/home/eventine/Documents/MATLAB/Results";
%path = "C:\Users\giovanni.bambini2\OneDrive - Alma Mater Studiorum Università di Bologna\___PhD\Papers\TAAS-2023\Results\Automated";

if isunix
	separator_os = "/";
else
	separator_os = "\";
end
				
%% 2: Define System Parameters: Init + Thermal

% Cores:
chip1.Nc = 16; %72;
chip1.Nh = 4; %9;			% Number of rows
chip1.Nv = 4; %8;			% Number of cols

% Version of the Thermal Model
chip1.model_ver = 0;
% Exponential Leakage
chip1.leak_exp = 1;
% Noise and Variation
chip1.pw_dev_per = 1;
chip1.param_dev_per = 1;
chip1 = chip1.create_model_deviation();

chip1.sensor_noise = 1;

% External Ambient Temperature
hpc.temp_amb = 25.0 + 273.15;


%% 3: Simulation setup

% SIMULATION TIME (in sec)
hpc.tsim = 2.0;
% Simulation Time-Step
chip1.Ts = 50e-6;
% Workload Quantum-step
chip1.quantum_us = 50;

% Input Time-Step
hpc.Ts_target = 1e-3;

% Create Thermal Model:
chip1.default_floorplan_config();
chip1.create_model_deviation();
chip1.model_init();

% Probability for each WL	
chip1.wl_prob = [0 2 3 2 2.5];	
% minimal execution in us
chip1.wl_min_exec_us = [50e3/8 5e3 5e3 5e3 100e3];
% mean execution in us
chip1.wl_mean_exec_us = [150e3/8 10e3 10e3 10e3 200e3];	

chip1.sensor_noise_amplitude = [1.0; 1.0; 1.0]; %og was 0.5
chip1.F_discretization_step = 0.05;

chip1.pw_gmean = 0.02;
chip1.pw_gvar = 1;
chip1.pw_glim = [-0.02 0.05];
chip1.pw_3sigma_on3 = 0.05;

chip1.create_core_pw_noise();

%% 4 Controller

% Max Temperature
chip1.core_limit_temp = 85 + 273.15;
chip1.core_crit_temp = 95 + 273.15;
% Max Power
tt = min(ceil(hpc.tsim / hpc.Ts_target)+1, (hpc.tsim/1e-4+1));
%Target Power Budget
hpc.chip_pw_budget{1} = 450/36*chip1.Nc*ones(tt,1);
tu = ceil(tt/4);
hpc.chip_pw_budget{1}(tu+1:2*tu) = 2*chip1.Nc;
hpc.chip_pw_budget{1}(2*tu+1:3*tu) = 5*chip1.Nc;
hpc.chip_pw_budget{1}(3*tu+1:3*tu+ceil(tu/2)) = 3*chip1.Nc;
hpc.chip_pw_budget{1}(3*tu+ceil(tu/2)+1:end) = 8*chip1.Nc;

hpc.toto_pw_budget = hpc.chip_pw_budget{1};

% Delays
chip1.delay_F_mean = 1e-5;		% Mean Frequency application Delay
chip1.delay_V_mean = 1e-5;		% Mean Voltage application Delay

%{
%hpc.wltrc{1} = zeros(size(hpc.wltrc{1}));
%hpc.wltrc{1}(:,hpc.ipl,:) = 1;
%hpc.wltrc{1}(:,2,:) = 1;
chip1.alpha_wl = 0.08; %0.4;

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
hpc.t_init{1} = (25+273.15)*ones(chip1.Ns,1); %temp_amb*ones(chip1.Ns,1);
%hpc.x_init(end-1:end,1) = temp_amb+20;
%hpc.tesim = 1.5;
%hpc.simulate_aut();

%hpc.exp_leakage = 1;
%hpc.simulate_aut();

%Others, TODO
%chip1.min_pw_red = 0.6;

%% ITERATE ON WORKLOAD
wl_times = 3;
th_models = 3; %2;
ndom = 4;
robust = 0;
show = 0;
savetofile = 0;
hpc.compare_vs_baseline = 1;
chip1.sensor_noise = 0;

test_iter = 2; %10;%20;

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

test_T_step = (max_init_T - hpc.temp_amb) /  (test_iter-1);
err_ampl = 0.5;
rand_T_init = rand(chip1.Ns,test_iter)*err_ampl*2 - err_ampl*ones(chip1.Ns,test_iter);
init_cond = ones(chip1.Ns,1) * (hpc.temp_amb:test_T_step:max_init_T);

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

	hpc.t_init{1} = init_cond(:,tit);

for wli=2:wl_times
	
	tic
		
	bwl = 1;
    tt = min(ceil(hpc.tsim / hpc.Ts_target)+1, (hpc.tsim/1e-4+1));
	switch wli
		case 1
			%full vector
			ll = ceil(hpc.tsim*1e6/chip1.quantum_us);
			hpc.wltrc{1} = zeros(chip1.Nc, chip1.ipl, ll);
			hpc.wltrc{1}(:,5,:) = 1;
			% Max Freq all time:
			hpc.frtrc{1} = 3.45 * ones(tt,chip1.Nc);
		case 2
			%full idle but 9 cores

		case 3
			hpc.wltrc{1} = wlwl3;
			% Max Freq all time:
			hpc.frtrc{1} = 3.45 * ones(tt,chip1.Nc);
		case 4
			hpc.wltrc{1} = wlwl4;
			% Max Freq all time:
			hpc.frtrc{1} = 3.45 * ones(tt,chip1.Nc);
	end

	if hpc.compare_vs_baseline
		% Compute ideal&unrestricted baseline
		%since atm it is domain independent:
			chip1.vd = chip1.Nc;
			chip1.VDom = eye(chip1.Nc);
		[wmaxc, perfmaxc] = hpc.base_ideal_unr(chip1,1);
		awl = hpc.perf_max_check{1};
		perf_max_check{wli} = awl;
	end
	
	%% ITERATE ON Thermal Models
	for mdli=1:th_models
	
		switch mdli
			case 1
				chip1.model_ver = 0;
				chip1.model_init();
			case 2
				chip1.model_ver = 1;
				chip1.model_init();
			case 3
				chip1.model_ver = 2;
				chip1.model_init();
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
					chip1.vd = 1;
					% Structure of the Voltage Domains Nc x vd
					chip1.VDom = ones(chip1.Nc, 1);
				case 2
					chip1.vd = 2;
					% Structure of the Voltage Domains Nc x vd
					chip1.VDom = zeros(chip1.Nc, chip1.vd);
                    chip1.VDom([1:8], 1) = 1;
					chip1.VDom([9:16], 2) = 1;
				case 3
					chip1.vd = 4;
					% Structure of the Voltage Domains Nc x vd
					chip1.VDom = zeros(chip1.Nc, chip1.vd);
					chip1.VDom([1,2,5,6], 1) = 1;
					chip1.VDom([3,4,7,8], 2) = 1;
					chip1.VDom([9,10,13,14], 3) = 1;
					chip1.VDom([11,12,15,16], 4) = 1;
				case 4
					chip1.vd = chip1.Nc;
					% Structure of the Voltage Domains Nc x vd
					chip1.VDom = eye(chip1.Nc);
				otherwise
					fprintf("Error! Wrong case selection!\n\r");
					chip1.vd = 1;
					% Structure of the Voltage Domains Nc x vd
					chip1.VDom = ones(chip1.Nc, 1);	
			end

			hpc.quad_pw_budget{1} = 450/36*chip1.Nc*ones(2,chip1.vd);
			itname = strcat('Model: ', int2str(mdli), ' - Domains: ', int2str(di), ' - WL: ', int2str(wli), ' - test: ', int2str(tit));

			%%

			%hpc.dummy_pw = 1;
            CM = [1];

			% New ControlPulp
			ctrl = Fuzzy(chip1);
			ctrl.C = chip1.Cc;
			ctrl.Ts_ctrl = 500e-6;
            ctrl.T_target = ones(ctrl.lNc, 1)*chip1.core_limit_temp;
            
            ctrl.comm_mat = 1;
            ctrl.multi_der_fcn = @mff;
            ctrl.multi_learn_rate = 0.1;
            ctrl.Tsensor_amplitude = chip1.sensor_noise*chip1.sensor_noise_amplitude(chip1.PVT_T);
            ctrl.Fdiscretization_step = chip1.F_discretization_step;
			
			simres = hpc.simulation(ctrl, chip1,CM,show);
			%xres{di, mdli, 1, wli} = xop;
			%ures{di, mdli, 1, wli} = uop; 
			%fres{di, mdli, 1, wli} = fop;
			%vres{di, mdli, 1, wli} = vop;
            wlop = simres(1).wlop;
			wlres{tit, di, mdli, 1, wli} = wlop;
			
			if show
				figarray = flip(findobj('Type','figure'));
				%TODO Improve: https://it.mathworks.com/help/matlab/ref/arrayfun.html
				for fn=figid:length(figarray)
					figarray(fn).Name = strcat(figarray(fn).Name, ' Alg: Fuzzy - ', itname);
				end
				figid = length(figarray)+1;
			end
			if savetofile
				hpc.saveall(simres(1).xop, simres(1).uop, simres(1).fop, simres(1).vop, simres(1).wlop, ...
                    path, regexprep(strcat('Alg: Fuzzy - ', itname), ':', '_'));
			else
				hpc.savetofile(hpc.xutplot(chip1, simres(1).xop, simres(1).uop),...
                    strcat(path,separator_os,regexprep(strcat('Alg: Fuzzy - ', itname), "-", "TP")),1);
			end
			
			%hpc.paperTPplot(hpc.Ts*[0:4000*10]', xop, uop);
	%%		
			% Classic MPC
			ctrl = cp_mpc(chip1);
			ctrl.C = eye(chip1.Ns);
			ctrl.Ts_ctrl = 1e-3;
            ctrl.comm_mat = 1;
            ctrl.multi_der_fcn = @mff;
            ctrl.multi_learn_rate = 0.1;
            ctrl.T_target = ones(ctrl.lNc, 1)*chip1.core_limit_temp;
            ctrl.Fdiscretization_step = chip1.F_discretization_step;
            ctrl.save_solver_stats = 0;
			simres = hpc.simulation(ctrl, chip1,CM,show);
			%xres{di, mdli, 2, wli} = xop;
			%ures{di, mdli, 2, wli} = uop; 
			%fres{di, mdli, 2, wli} = fop;
			%vres{di, mdli, 2, wli} = vop;
            wlop = simres(1).wlop;
			wlres{tit, di, mdli, 2, wli} = wlop;			
			
			if show
				figarray = flip(findobj('Type','figure'));
				%TODO Improve: https://it.mathworks.com/help/matlab/ref/arrayfun.html

				for fn=figid:length(figarray)
					figarray(fn).Name = strcat(figarray(fn).Name, ' Alg: CP - ', itname);
				end
				figid = length(figarray)+1;
			end
			
			if savetofile
				hpc.saveall(simres(1).xop, simres(1).uop, simres(1).fop, simres(1).vop, simres(1).wlop, ...
                    path, regexprep(strcat('Alg: CP - ', itname), ':', '_'));
            end

		end % for di domains
	end % for mdli thermal model
	toc
	
end % for wl_times

end %test_iter



function [x] = mff(v1,id)
    x = 0;
end