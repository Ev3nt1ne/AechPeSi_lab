%
% This Script extracts the MPC model by extracting the yalmip objective and constraints
% for different problem sizes and saves QP vectors and matrices in .mat files.
%

% TODO: change these path if you are on a different system
clear variables;
addpath('../');                      % add git repo base to search path
addpath('../Controllers');                      % add git repo base to search path
basename = '../../src';   % directory were .mat files are saved to

experiments = ...
[ ...
  %3  3  3; ...
  %4  4  3; ...
  %5  5  3; ...
  %6  6  3; ...
  %7  7  3; ...
  %8  8  3; ...
  %9  9  3; ...
  %10  10  3; ...
  %11  11  3; ...
  %12  12  3; ...
  3  3  2; ...
  4  4  2; ...
  5  5  2; ...
  6  6  2; ...
  7  7  2; ...
  8  8  2; ...
  9  9  2; ...
  10  10  2; ...
  11  11  2; ...
  12  12  2; ...
  3  3  3; ...
  6  6  3; ...
  9  9  3; ...
  12  12  3; ...
  3  3  4; ...
  6  6  4; ...
  9  9  4; ...
  12  12  4; ...
  %6 7 7;
  ];

tic;
for values = experiments.'
    %% Create HPC problem:
    hpc = hpc_lab;

    % Rows and Columns
    hpc.Nv = values(2);     % Number of Processing Elements in each row
    hpc.Nh = values(1);     % Number of PE's in each column

    % Core number
    hpc.Nc = hpc.Nh*hpc.Nv;         % Total Number of controlled PE's in thermal model
    
    % Number of Voltage Domains
    hpc.vd = hpc.Nc;

    % How the cores are distributed per domain
    %hpc.VDom = ... ;
    % Alternatively
    hpc.default_VDom_config();

    % Default Floorplan and Parameter deviation vector
    hpc.default_floorplan_config();
    hpc.create_model_deviation();

    %Decide simulation frequency (discrete timing)
    hpc.Ts = 5e-5;

    %Decide total simulation time
    hpc.tsim = 0;

    %Presence/Absence of sensor noise:
    hpc.sensor_noise = 1;
    
    %Thermal model version
    hpc.model_ver = 0;
    
    hpc.model_init();
    %hpc.create_core_pw_noise();
    
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
    
    %% MPC
    ctrl = cp_mpc();
    ctrl.C = eye(hpc.Ns);

    % options
    ctrl.ops.warm_start = true;

    % Control Horizon
    ctrl.Nhzn = values(3);
    
    %Controller frequency (s)
    ctrl.Ts_ctrl = 5e-3;
    
    %Robust Margins for T and P
    ctrl.Cty = zeros(ctrl.Nhzn, hpc.Nc);
    ctrl.Ctu = zeros(ctrl.Nhzn, hpc.Nc);
    
    %Reference Tracking Objective Matrix
    R_coeff = 10;
    ctrl.R = R_coeff*eye(hpc.Nc);
    
    %Others
    ctrl.R2 = zeros(hpc.Nc);
    ctrl.Q = zeros(hpc.Ns);
    
    %% 4 verbose
    % check
    %hpc.anteSimCheckTM();
    %hpc.anteSimCheckLab();
    %hpc.anteSimCheckPM();
    
    %yalmip problem
    hpc.init_compute_model(hpc.Ad_true, hpc.Bd_true);
    Nsim = 1; %Nsim = ceil(obj.tsim / ctrl.Ts_ctrl);
    ctrl = ctrl.init_fnc(hpc, Nsim);
    hpc
    ctrl

    %% 5 dump Quadratic Program matrix data for further optimization, code generation
    filename = char(compose('%s/_HPC_%dx%d_H%d.mat',basename, hpc.Nh, hpc.Nv, ctrl.Nhzn));
    display(strcat('Saving model to file:  ',filename));
    ops = sdpsettings('verbose',1,'solver','osqp', 'usex0',0);
    ctrl.ylmp_constraints;
    yalmip2mat(filename,ctrl.ylmp_constraints,ctrl.ylmp_objective,ops,ctrl.ylmp_opt_variables,ctrl.ylmp_opt_output);
    %yalmip_model = export(constraints,hpcective,ops,opt_variables,opt_output);
    %a = yalmip2osqp(yalmip_model);
    %ops = sdpsettings('verbose',1,'solver','osqp', 'savesolverinput',1);
    %qp = optimize(constraints,hpcective,ops);
end
toc;
