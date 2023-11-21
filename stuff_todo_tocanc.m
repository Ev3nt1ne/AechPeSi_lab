
%{ 
% Missing:



properties(SetAccess=protected, GetAccess=public)		
	polyFV_opt;
end

freq_fact = 8;				% Times per seconds expected Target Freq Changes

		Ni_nl;						% Number of inputs for the non-linear case
		Ni_c;						% Number of controllable inputs
		Ni_nc;						% Number of non-controllable inputs





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
addpath Controllers/

hpc.t_init = hpc.temp_amb*ones(hpc.Ns,1);
hpc.wltrc = hpc.generate_wl_trace(hpc.Nc, hpc.tsim, 0);
%hpc.sim_tm_autonomous()
wl_bkp = hpc.wltrc;

%%
ctrl = Fuzzy;
ctrl.C = hpc.C;
%%
ctrl = black_wolf;
ctrl.Ts_ctrl = 5e-3;
ctrl.C = eye(hpc.Ns);
%TODO
ctrl.Cty = zeros(ctrl.Nhzn, hpc.Nout);
ctrl.Ctu = zeros(ctrl.Nhzn, hpc.Nc);
ctrl.R = eye(hpc.Nc);
ctrl.R2 = zeros(hpc.Nc);
ctrl.Q = zeros(hpc.Ns);
%%
hpc.x_init = hpc.temp_amb * ones(hpc.Ns,1);
hpc.frtrc = 3.45 * ones(min(ceil(hpc.tsim / hpc.Ts_target)+1,(hpc.tsim/ctrl.Ts_ctrl+1)), hpc.Nc);
ts = ceil(hpc.tsim / hpc.Ts_target)+1;
hpc.tot_pw_budget = 450/36*hpc.Nc*ones(ts,1);
ts = ceil(ts/4);
hpc.tot_pw_budget(ts+1:2*ts) = 2*hpc.Nc;
hpc.tot_pw_budget(2*ts+1:3*ts) = 5*hpc.Nc;
hpc.tot_pw_budget(3*ts+1:3*ts+ceil(ts/2)) = 3*hpc.Nc;
hpc.tot_pw_budget(3*ts+ceil(ts/2)+1:end) = 8*hpc.Nc;
hpc.quad_pw_budget = 450/36*hpc.Nc*ones(2,hpc.vd);

hpc.min_pw_red = 0.6;
%%
hpc.simulation(ctrl,1);

%%
ll = ceil(hpc.tsim*1e6/hpc.quantum_us);
hpc.wltrc = zeros(hpc.Nc, hpc.ipl, ll);
hpc.wltrc(:,5,:) = 1;

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


%% MUL and DIV simulations
clc;

Ts_target = 5.99;
Ts_ctrl = 1;
N = 10060;

% Inputs
% initial offset:
target_ts_offset = 0;
% frequency
ts_div = Ts_target / Ts_ctrl;
target_mul = floor(ts_div);
ctrl_mul = floor(1/ts_div) - 1;
ctrl_mul = (ctrl_mul<0)*0 + (ctrl_mul>=0)*ctrl_mul;
target_counter = target_mul + target_ts_offset;
target_index = 0;
% Manage non-integer part:
if target_mul>ctrl_mul
	ntd = (ts_div-target_mul)/target_mul/target_mul;
	nip_sign = -1;
else
	ntd = (1/ts_div)-(ctrl_mul+1);
	nip_sign = +1;
end
nip_mul = ceil(1/ntd);
nip_counter = nip_mul + target_ts_offset;

for s=1:N

	% Input step managing
	target_counter = (target_counter-1)*(target_counter>0) + (target_counter<=0)*(target_mul-1);
	nip_counter = (nip_counter-1)*(nip_counter>0) + (nip_counter<=0)*(nip_mul-1);
	target_index = target_index + (target_counter==target_mul-1) + nip_sign*(nip_counter==nip_mul-1) + ctrl_mul;

	%disp(strcat("s: ",int2str(s),", cnt: ",int2str(target_counter)," / ", int2str(nip_counter), ", idx: ",int2str(target_index)));

end
disp(strcat("s: ",int2str(s),", cnt: ",int2str(target_counter)," / ", int2str(nip_counter), ", idx: ",int2str(target_index)));

disp(num2str(N/ts_div))

disp(strcat("ts_div: ",num2str(ts_div), ", target_mul: ",int2str(target_mul), ", ctrl_mul: ", int2str(ctrl_mul)))
disp(strcat("nip_mul: ", num2str(nip_mul), ", ntd: ", num2str(ntd)))


%%
%Taking a big sample for proper testing
A = rand(2500,500,100);
v = randi(size(A,3),size(A,1),1);

fun1 = @() forloop(A,v);
fun2 = @() vectorization(A,v);

tic
z1 = fun1();
toc

z2 = fun2();


%Check for equality
isequal(z1,z2)

%%

function F = forloop(A,v)
[m,n,l] = size(A);
F = zeros(m,n);
for k=1:m
    F(k,:) = A(k,:,v(k));
end

end
function F = vectorization(A,v)
[m,n,l] = size(A);
tic
idx = sub2ind([m,n,l],repelem(1:m,1,n),repmat(1:n,1,m),repelem(v(:).',1,n));
toc
tic
F = (A(idx));
F = reshape(F,[],m).';
toc
end
