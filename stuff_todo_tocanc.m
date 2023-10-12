
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

hpc.tm.Nc

