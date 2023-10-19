
%{ 
% Missing:

quantum_us = 50;

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

freq_fact = 8;				% Times per seconds expected Target Freq Changes




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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHODS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dxdt, power] = nl_model_dyn(obj,A,B,x,u,d_i,d_p, ot) %, pw, ceff, pw_levels)

			Ni_f = obj.Nc;
			Ni_v = obj.vd;
			
			dxdt = zeros(obj.Ns,1);

			F = u(1:Ni_f);
			V = obj.VDom*u(Ni_f+1:Ni_f+Ni_v);
			
			power = obj.power_compute(F,V,obj.C(1:obj.Nc,:)*x,d_i,d_p);
			dxdt = A*x + B*[power; ot];
		end

%%%%%%



%}


%%
hpc = hpc_lab;

hpc.t_init = hpc.temp_amb*ones(hpc.Ns,1);
hpc.sim_tm_autonomous()

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

%%

pw = power_model;
pw.pw_dev_per = [1; 1.5; 0.5; 10; 0.1];

%{
pw.power_compute(1,1,300,[0 0.5 0.5 0 0], 1, 0)
pw.power_compute([1;1;1],1,[300;300;300],[0 0.5 0.5 0 0], 1, 0)
pw.power_compute(1,[1;1;1],[300;300;300],[0 0.5 0.5 0 0], 1, 0)
pw.power_compute(1,1,[300;300;300],[0 0.5 0.5 0 0], [1;1;1], 0)
pw.power_compute(1,1,[300;300;300],ones(3,1)*[0 0.5 0.5 0 0], 1, 0)
pw.power_compute([1;1;1],1,[300;300;300],[0 0.5 0.5 0 0], 1, 1)
pw.power_compute(1,[1;1;1],[300;300;300],[0 0.5 0.5 0 0], 1, 1)
pw.power_compute(1,1,[300;300;300],[0 0.5 0.5 0 0], [1;1;1], 1)
pw.power_compute(1,1,[300;300;300],ones(3,1)*[0 0.5 0.5 0 0], 1, 1)

pw.power_compute([1;1;1],[1;1;1],[300;300;300],[0 0.5 0.5 0 0], 1, 0)
pw.power_compute(1,[1;1;1],[300;300;300],ones(3,1)*[0 0.5 0.5 0 0], 1, 0)

pw.power_compute([1;1;1],[1;1;1],[300;300;300],[0 0.5 0.5 0 0], 1, 1)
pw.power_compute(1,[1;1;1],[300;300;300],ones(3,1)*[0 0.5 0.5 0 0], 1, 1)
%}




