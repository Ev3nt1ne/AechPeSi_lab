
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
hpc.Ts = 1e-5;
hpc.tsim = 2;
addpath Controllers/

hpc.t_init = hpc.temp_amb*ones(hpc.Ns,1);
%%
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
hpc.sensor_noise = 1;
%TODO
ctrl.Cty = zeros(ctrl.Nhzn, hpc.Nc);
ctrl.Ctu = zeros(ctrl.Nhzn, hpc.Nc);
R_coeff = 10; %10;
R2_coeff = R_coeff/1; %100; %30; %10

ctrl.R = R_coeff*eye(hpc.Nc);
ctrl.R2 = zeros(hpc.Nc);

for v=1:hpc.vd
	%Here I could create an accumulation thing that need to be
	%optimized (e.g. reduced) that contains the deltaF among quadrant 
	ctrl.R2 = ctrl.R2 + hpc.VDom(:,v)*hpc.VDom(:,v)' / sum(hpc.VDom(:,v))^2 * R2_coeff;
end
%TODO: Assuming that the diagonal is full
ctrl.R2(~eye(size(ctrl.R2))) = ctrl.R2(~eye(size(ctrl.R2))) * (-1);
cores_per_dom = sum(hpc.VDom);
ctrl.R2(logical(eye(size(ctrl.R2)))) = ctrl.R2(~~eye(size(ctrl.R2))) .* hpc.VDom*cores_per_dom'; %here I need ~~ to converto to logical values

%hpc.R2 = hpc.R2 + (R_coeff/10)*eye(size(hpc.R2));
ctrl.Q = zeros(hpc.Ns);
%%
ctrl = cp_mpc;
ctrl.Ts_ctrl = 5e-3;

ctrl.C = eye(hpc.Ns);
hpc.sensor_noise = 1;
%TODO
ctrl.Cty = zeros(ctrl.Nhzn, hpc.Nc);
ctrl.Ctu = zeros(ctrl.Nhzn, hpc.Nc);
R_coeff = 10; %10;

ctrl.R = R_coeff*eye(hpc.Nc);
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
tic;
hpc.simulation(ctrl,1);
toc;
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

%{
obj.psoff_lut = [
	1.0924    0.7147    0.1889   -0.4509   -1.1900   -1.9638   -2.5989
	0.8312    0.4534   -0.0710   -0.7070   -1.4400   -2.1963   -2.7762
	0.5706    0.1929   -0.3301   -0.9616   -1.6873   -2.4235   -2.9392
	0.3110   -0.0668   -0.5880   -1.2144   -1.9317   -2.6445   -3.0856
	0.0523   -0.3254   -0.8446   -1.4651   -2.1727   -2.8582   -3.2127
		-0.2052   -0.5829   -1.0997   -1.7133   -2.4096   -3.0635   -3.3175
		-0.4613   -0.8391   -1.3531   -1.9586   -2.6418   -3.2590   -3.3963
		-0.7159   -1.0936   -1.6044   -2.2007   -2.8686   -3.4431   -3.4449
		-0.9686   -1.3464   -1.8535   -2.4388   -3.0891   -3.6141   -3.4584
		-1.2192   -1.5970   -2.0997   -2.6726   -3.3022   -3.7698   -3.4313
		-1.4674   -1.8451   -2.3429   -2.9011   -3.5069   -3.9077   -3.3569
		-1.7127   -2.0904   -2.5823   -3.1236   -3.7017   -4.0249   -3.2276
		-1.9546   -2.3324   -2.8175   -3.3391   -3.8850   -4.1182   -3.0346
		-2.1927   -2.5704   -3.0478   -3.5464   -4.0550   -4.1836   -2.7676
		-2.4263   -2.8041   -3.2723   -3.7444   -4.2095   -4.2167   -2.4146
];

Fv = [0.4000    0.9000    1.4000    1.9000    2.4000    2.9000    3.4000]';
Tv = [20    25    30    35    40    45    50    55    60    65    70    75    80    85    90]';
obj.F0v = ones(hpc_class.Nc, length(Fv))*diag(Fv);
obj.T0v = ones(hpc_class.Nc, length(Tv))*diag(Tv);
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

hpc.pws_ls_offset(4, 11, 10)


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
