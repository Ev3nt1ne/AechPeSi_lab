
%{ 
% Missing:

quantum_us = 50;
customColormap;
observable;

properties(SetAccess=protected, GetAccess=public)		
	Controller;					% Persistent variable to optimize controller
	polyFV_opt;
	
	% wl:
	% How much? coefficient
	dur_c = 5/6;
	min_secondary_space = 0.075;
	
	% model
	 wl_index;
	 qt_storage;
	 V_s;
	 F_s;
	 A_s;
	 B_s;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In the constructor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%obj.x_init = ones(obj.Ns,1)*obj.temp_amb;
% TODO: remove this after leakage is fixed:
leak_store = obj.exp_leakage;
obj.exp_leakage = 0;
%obj.min_pw_red = obj.power_compute(ones(obj.Nc,1)*obj.F_min,ones(obj.Nc,1)*obj.V_min,(obj.core_crit_temp)*ones(obj.Nc,1),[zeros(obj.Nc,obj.ipl-1) ones(obj.Nc,1)],ones(obj.Nc,1));
obj.min_pw_red = obj.power_compute(ones(obj.Nc,1)*obj.F_min,ones(obj.Nc,1)*obj.V_Max,(obj.core_crit_temp)*ones(obj.Nc,1),[zeros(obj.Nc,obj.ipl-1) ones(obj.Nc,1)],ones(obj.Nc,1));
obj.exp_leakage = leak_store;

%}

%%

a = thermal_model;

%a.VDom = [1, NaN, inf, 'l'];

%a.Ac_nom = [1 1;1 1];

a.C = ones(a.Nout, a.Ns);

%% Things to put in the constructor
Ac_nom
Bc_nom
Ac_true
Bc_true
C
D
param_dev_per

% things to ALSO fix in set
RC_floorplan_length
CPw_floorplan_length
R_floorplan_material
C_floorplan_material
sensor_active

%
add_states
t_comp
wid_comp
len_comp


%%
hpc = hpc_lab;

hpc.t_init = hpc.temp_amb*ones(hpc.Ns,1);
hpc.sim_tm_autonomous()

%%
tm = thermal_model();

tm.C

tm.sensors_active = [1 1 1 0];

tm = tm.model_init();
tm.C
%%
tnA = hpc.Ac_nom;

hpc2 = hpc_system();

% Cores:
hpc2.Nc = 9; %72;
hpc2.Nh = 3; %9;			% Number of rows
hpc2.Nv = 3; %8;			% Number of cols

% Version of the Thermal Model
hpc2.thermal_model_ver = 0;
% Exponential Leakage
hpc2.exp_leakage = 1;
% Noise and Variation
hpc2.model_variation = 1;
hpc2 = hpc2.create_core_pw_noise();
hpc2.measure_noise = 1;

% External Ambient Temperature
temp_amb = 25.0 + 273.15;

% Create Thermal Model:
hpc2 = hpc2.create_thermal_model();

toA = hpc2.Ac_nom;

dif = tnA - toA


difB = hpc.Bc_nom - hpc2.Bc_nom



