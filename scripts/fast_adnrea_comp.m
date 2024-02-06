clc;
clear variables;

%% INIT

% Instantiate main Class:
hpc = hpc_lab;

% Core number
hpc.Nc = 36;

% Rows and Columns
hpc.Nv = 6;
hpc.Nh = 6;

% Number of Voltage Domains
hpc.vd = 4;

% How the cores are distributed per domain
hpc.VDom = zeros(hpc.Nc, hpc.vd);
hpc.VDom([1:3 7:9 13:15], 1) = 1;
hpc.VDom([4:6 10:12 16:18], 2) = 1;
hpc.VDom([19:21 25:27 31:33], 3) = 1;
hpc.VDom([22:24 28:30 34:36], 4) = 1;

% Default Floorplan and Parameter deviation vector
hpc.default_floorplan_config();
hpc.create_model_deviation();

%Decide simulation frequency (discrete timing)
hpc.Ts = 5e-5;

%Decide total simulation time
hpc.tsim = 2;

%Presence/Absence of sensor noise:
hpc.sensor_noise = 1;

%Thermal model version
hpc.model_ver = 0;
hpc.model_init();
hpc.create_core_pw_noise();

% Initial Temperature condition
hpc.t_init = hpc.temp_amb*ones(hpc.Ns,1);

%Input maximum frequency
hpc.Ts_target;

tt = min(ceil(hpc.tsim / hpc.Ts_target)+1, (hpc.tsim/1e-4+1));

% Target Frequency Trajectory
hpc.frtrc = 3.45 * ones(tt, hpc.Nc);

%Target Power Budget
hpc.tot_pw_budget = 450/36*hpc.Nc*ones(tt,1);
tu = ceil(tt/4);
hpc.tot_pw_budget(tu+1:2*tu) = 2*hpc.Nc;
hpc.tot_pw_budget(2*tu+1:3*tu) = 5*hpc.Nc;
hpc.tot_pw_budget(3*tu+1:3*tu+ceil(tu/2)) = 3*hpc.Nc;
hpc.tot_pw_budget(3*tu+ceil(tu/2)+1:end) = 8*hpc.Nc;
hpc.quad_pw_budget = 450/36*hpc.Nc*ones(2,hpc.vd);

%Others, TODO
hpc.min_pw_red = 0.6;


hpc.anteSimCheckTM();
hpc.anteSimCheckLab();
hpc.anteSimCheckPM();


%% Generate Workload trace:
%Probability of each workload type [idle, int, float, mem, vec]
hpc.wl_prob = [0.0345 0.3448 0.3448 0.2414 0.0345];

%Minimum execution time (us) of each workload type
hpc.wl_min_exec_us = [510 100 100 480 520];	

%Mean execution time (us) of each workload type
hpc.wl_mean_exec_us = [520e2 210e1 180e1 560e2 1.2e6/1e1];

%Generate it
hpc.wltrc = hpc.generate_wl_trace(hpc.Nc, hpc.tsim, 0);


%% Fuzzy Controller
ctrl = Fuzzy;
ctrl.C = hpc.C;

%Controller frequency (s)
ctrl.Ts_ctrl = 5e-4;


%Launch Simulation
tic;
hpc.simulation(ctrl,1);
toc;

%% CP
ctrl = CP;
ctrl.C = hpc.C;

%Controller frequency (s)
ctrl.Ts_ctrl = 5e-4;


%Launch Simulation
tic;
hpc.simulation(ctrl,1);
toc;