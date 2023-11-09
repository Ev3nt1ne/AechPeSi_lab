%
% This Script extracts the MPC model by extracting the yalmip objective and constraints
% for different problem sizes and saves QP vectors and matrices in .mat files.
%

% TODO: change these path if you are on a different system
clear all;
addpath('../');                      % add git repo base to search path
basename = '../../src';   % directory were .mat files are saved to

experiments = ...
[ ...
  6  6  3; ...
  7  7  3; ...
  8  8  3; ...
  9  9  3; ...
  10  10  3; ...
  11  11  3; ...
  12  12  3; ...
  6  6  2; ...
  7  7  2; ...
  8  8  2; ...
  9  9  2; ...
  10  10  2; ...
  11  11  2; ...
  12  12  2; ...
  6 7 7;
  ];

for values = experiments.'
    
    %% 1: only instantiate
    Nh = values(1); 			% Number of Processing Elements in each row
    Nv = values(2); 			% Number of PE's in each column
    Nc = Nh*Nv;         % Total Number of controlled PE's in thermal model
    hpc = hpc_system(Nc,Nh,Nv);
    hpc.Nhzn = values(3);           % Control Horizon
    
    %% 2: Define System Parameters: Init + Thermal
    % Cores:
    display(char(compose('Generating linear thermal model with Nh=%d, Nv=%d and horizon=%d',hpc.Nh, hpc.Nv, hpc.Nhzn)));
    hpc % verbose output:
    
    % Version of the Thermal Model
    % hpc.thermal_model_ver = 0;
    % Exponential Leakage
    % hpc.exp_leakage = 1;
    % Noise and Variation
    %hpc.model_variation = 1;
    
    % hpc.pw_gmean = 0.02;
    % hpc.pw_gvar = 1;
    % hpc.pw_glim = [-0.02 0.05];
    % hpc.pw_3sigma_on3 = 0.05;
    % hpc = hpc.create_core_pw_noise();
    % hpc.measure_noise = 1;
    
    % External Ambient Temperature
    temp_amb = 25.0 + 273.15;
    
    % Create Thermal Model:
    hpc = hpc.create_thermal_model();
    
    hpc.x_init = temp_amb*ones(hpc.Ns,1);
    
    %% 3: Simulation setup
    
    % SIMULATION TIME (in sec)
    hpc.tsim = 2.0;
    % Simulation Time-Step
    hpc.Ts = 50e-6;
    % Workload Quantum-step
    hpc.quantum_us = 50;
    
    % Controller Time-Step
    hpc.Ts_ctrl = 500e-6;
    % Input Time-Step
    hpc.Ts_input = 1e-3;
    
    % Probability for each WL	
    hpc.wl_prob = [0 2 3 2 2.5];	
    % minimal execution in us
    hpc.wl_min_exec_us = [50e3/8 5e3 5e3 5e3 100e3];
    % mean execution in us
    hpc.wl_mean_exec_us = [150e3/8 10e3 10e3 10e3 200e3];
    
    % hpc = hpc.generate_reference();
    
    hpc.T_noise_max = 0.5;
    hpc.F_discretization_step = 0.05;
    
    
    %% 4 Controller
    
    % Max Temperature
    hpc.core_crit_temp = 85 + 273.15;
    % Max Power
    ts = ceil(hpc.tsim / hpc.Ts_input)+1;
    hpc.tot_pw_budget = 450*ones(ts,1);
    ts = ceil(ts/4);
    hpc.tot_pw_budget(ts+1:2*ts) = 4*hpc.Nc;
    hpc.tot_pw_budget(2*ts+1:3*ts) = 2*hpc.Nc;
    hpc.tot_pw_budget(3*ts+1:3*ts+ceil(ts/2)) = 6*hpc.Nc;
    hpc.tot_pw_budget(3*ts+ceil(ts/2)+1:end) = 8*hpc.Nc;
    
    % TODO: ??? why
    % TODO: any constraint change is only allowed to change the l, u vector 
    %       not the A matrix -> otherwise numerical refactorization has to be done.
    % TODO: but changes in l, u might infringe on recursive feasibility!
    %       Actually total power budget should be a soft constraint to avoid problems in that regard.
    % TODO: insentive on setpoint should be given by the reference!
    hpc.tot_pw_budget = hpc.tot_pw_budget/36*hpc.Nc;
    % Delays
    hpc.delay_F_mean = 1e-5;		% Mean Frequency application Delay
    hpc.delay_V_mean = 1e-5;		% Mean Voltage application Delay
    
    %hpc.frplot = 3.45 * ones(min(ceil(hpc.tsim / hpc.Ts_input)+1,(hpc.tsim/hpc.Ts_ctrl+1)), hpc.Nc);
    
    %% MPC init:
    
    % Robust Controls
    robust = 1;
    
    R_coeff = 10; %10;
    R2_coeff = 100; %30; %10
    
    hpc.R = R_coeff*eye(hpc.Nc);
    hpc.R2 = zeros(hpc.Nc);
    
    for v=1:hpc.vd
    	%Here I could create an accumulation thing that need to be
    	%optimized (e.g. reduced) that contains the deltaF among quadrant 
    	hpc.R2 = hpc.R2 + hpc.VDom(:,v)*hpc.VDom(:,v)' / sum(hpc.VDom(:,v))^2 * R2_coeff;
    end
    %TODO: Assuming that the diagonal is full
    hpc.R2(~eye(size(hpc.R2))) = hpc.R2(~eye(size(hpc.R2))) * (-1);
    cores_per_dom = sum(hpc.VDom);
    hpc.R2(logical(eye(size(hpc.R2)))) = hpc.R2(~~eye(size(hpc.R2))) .* hpc.VDom*cores_per_dom'; %here I need ~~ to converto to logical values
    
    hpc.yMax = ones(hpc.Nout,1)*hpc.core_crit_temp;
    % TODO: check that core power constraint is not unbounded
    hpc.umin = hpc.core_min_power;
    hpc.uMax = hpc.core_Max_power;
    hpc.usum = [hpc.tot_pw_budget(1) ...
    					hpc.quad_pw_budget(1,:)];
    
    hpc.mpc_robustness = 0.7;
    
    %{
    if robust == 1
    	hpc.Cty = hpc.find_unc_set(hpc.mpc_robustness*100)*hpc.C';
    	%hpc.Ctu = hpc.max_uncertainty(??,hpc.mpc_robustness*100)*ones(hpc.Nhzn, hpc.Ni_c);
    	hpc.Ctu = repelem(hpc.mpc_robustness* ...
    		abs(hpc.power_compute(hpc.F_Max*ones(hpc.Nc,1),hpc.V_Max*ones(hpc.Nc,1),ones(hpc.Nc,1)*hpc.core_crit_temp,[1 zeros(1,hpc.ipl-1)],1) ...
    		- hpc.power_compute(hpc.F_Max*ones(hpc.Nc,1),hpc.V_Max*ones(hpc.Nc,1),ones(hpc.Nc,1)*hpc.core_crit_temp,[zeros(1,hpc.ipl-1) 1],1) ), ...
    		1, hpc.Nhzn)';
    %}
    if robust == 1
    	hpc.Cty = 5 * ones(hpc.Nhzn, hpc.Nout);
    	hpc.Ctu = 0.3 * ones(hpc.Nhzn, hpc.Ni_c);
    else
    	hpc.Cty = zeros(hpc.Nhzn, hpc.Nout);
    	hpc.Ctu = zeros(hpc.Nhzn, hpc.Ni_c);
    end
    
    %%
    hpc.R2 = zeros(hpc.Nc);
    
    hpc.exp_leakage = 0;
    hpc.ctrl_MA = 0;
    hpc.iterative_fv = ~hpc.ctrl_MA;
    
    %% linear mpc
    hpc = hpc.lin_mpc_setup();
    
    
    % dump Quadratic Program matrix data for further optimization, code generation
    filename = char(compose('%s/_HPC_%dx%d_H%d.mat',basename, hpc.Nh, hpc.Nv, hpc.Nhzn));
    display(strcat('Saving model to file:  ',filename));
    ops = sdpsettings('verbose',1,'solver','osqp', 'usex0',0);
    hpc.ylmp_constraints;
    yalmip2mat(filename,hpc.ylmp_constraints,hpc.ylmp_objective,ops,hpc.ylmp_opt_variables,hpc.ylmp_opt_output);
    %yalmip_model = export(constraints,hpcective,ops,opt_variables,opt_output);
    %a = yalmip2osqp(yalmip_model);
    %ops = sdpsettings('verbose',1,'solver','osqp', 'savesolverinput',1);
    %qp = optimize(constraints,hpcective,ops);
 end
