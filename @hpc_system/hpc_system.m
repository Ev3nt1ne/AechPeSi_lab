
% Requires MPT Toolbox: https://www.mpt3.org/
% Requires Ipopt: https://github.com/ebertolazzi/mexIPOPT

classdef hpc_system
	
	properties		
		%% System Structure
		Nc = 36;					% Number of Cores
		Nh = 6;						% Number of rows
		Nv = 6;						% Number of cols
		vd = 4;						% Voltage/Power Domains
		VDom;						% Structure of the Voltage Domains Nc x vd
				
		Ts_ctrl = 5e-4;				% Controller Ts
		Ts_obs = 1e-3;				% Observer Ts
		Ts_mpc = 5e-3;				% MPC Ts
		Ts_input = 1e-3;			% Commands min Ts
		Obs_poles = [0.8 0.1];		% Poles of the Luemberg Observer
		
		%% Matlab/Simulation
		thermal_model_ver = 0;		% Version of the Thermal Model
		model_variation = 1;		% If the true system should be != from nominal
		exp_leakage = 1;
		measure_noise = 1;
		T_noise_max = 1.0;
		
		tsim = 1;					% Simulation Time in [s]
		tesim = 2;					% Free Simulation Time in [s]
		Ts = 250e-6;				% Discretization Time
		x_init;						% Initial Conditions 
		urplot;						% Input Reference Plot
		frplot;						% Freq Reference Plot
		zrplot;						% Power Noise Plot
		wrplot;						% Wokrload Plot
		freq_fact = 8;				% Times per seconds expected Target Freq Changes
		
		quantum_us = 50;
		
		% Matlab Matrices	
		Ac_nom;
		Bc_nom;
		Ac_true;
		Bc_true;
		
		C;
		D;
		
		observable;
		
		customColormap;
		
		%% All Controllers
		tot_pw_budget;
		quad_pw_budget;
		
		%% MPC Controller
		separate_omega;
		xref;
		uref;
		yref;
		usum;
		Q;
		R;
		R2;
		Nhzn = 3;
		Ctx;
		Ctu;
		Cty;
		
		mpc_robustness = 1;
		
		umin;
		uMax;
		xmin;
		xMax;
		ymin;
		yMax;
        ylmp_constraints;                % save yalmip constraints 
        ylmp_objective;                   % save yalmip objective
        ylmp_opt_variables;              % save the yalmip optimization variables used
        ylmp_opt_output;                 % save the yalmip variable that is extracted after optimization
        
		
		
		mpc_h0_T;
		mpc_h0_F;
		mpc_h0_app;
		
		
		%% SoA/PID Controller
		kp = 1.9518;				% PID Kp
		ki = 73.931;				% PID Ki
		aw_up = 0.05;				% PID Anti-Windup Upper Saturation
		aw_down_c = -0.75;			% PID Anti-Windup Down Sat. Coefficient
		T_margin = 5;				% PID Temperature Margin
		voltage_rule = 100;			% Percentile in the Choice of the Voltage in the domain

		sat_up = 0;					% PID Upper Saturation
		
		pid_e_down = 2.5;			% PID banding, down constraint (absolute value)
		pid_e_up = 1;				% PID banding, up constraint (absolute value)
		pid_e_band_coeff = 0.7;		% PID banding coefficient
		
		ctrl_MA = 1;				% If CP should use Moving Average Voltage Selection
		ctrl_fixedv = 0;			% if CP should use Fixed Voltage = V_Max
		ctrl_fixedv_value = 1.2;	% the value of the fixed Voltage
		iterative_fv = 0;			% iterative fv solution. Generally: = ~ctrl_MA
		pw_inverse = 0;				% use Inverse pw reduction wrt temp
		
		alpha_wl = 0.4;				% Moving Average filter parameter for Workload
		
		dummy_pw = 0;
		
		%% System Parameters
		temp_amb = 25.0 + 273.15;	% External Ambient Temperature

		leak_vdd = 324.29;
		leak_temp = 1.0;
		leak_process = 528.04;
		
		ceff_pw_coeff = [306.741, 694.866, 1235.452, 1651.217, 1600.92];
		%[466.741, 694.866, 1235.452, 1651.217, 1600.92];
		exp_leak_coeff = [4507.748, 29.903, 6033.035];
		
		min_pw_red;

		core_crit_temp = 358.15;

		FV_table = [0.50, 0.40, 1.3500;
				 0.55, 0.40, 1.6000;
				 0.60, 1.35, 1.8000;
				 0.65, 1.35, 2.0000;
				 0.70, 1.35, 2.2000;
				 0.75, 1.60, 2.4000;
				 0.80, 1.80, 2.6000;
				 0.85, 1.80, 2.7500;
				 0.90, 2.00, 2.9000;
				 0.95, 2.00, 3.0500;
				 1.00, 2.00, 3.2000;
				 1.05, 2.60, 3.3500;
				 1.10, 2.60, 3.4500;
				 1.15, 2.60, 3.5500;
				 1.20, 2.60, 3.6600];
			 
		% Power Noise
		%pn_mu = 0.1;
		%pn_sigma = 2.0; %0.5
		%pn_alpha = 0.2; %0.8
		pw_gmean = 0.02;
		pw_gvar = 1;
		pw_glim = [-0.02 0.05];
		pw_3sigma_on3 = 0.05;
		
		core_pw_noise_char;
		ThMoNoise;
		
		%
		delay_F_mean = 1e-5;		% Mean Frequency application Delay
		delay_V_mean = 1e-5;		% Mean Voltage application Delay
		delay_F_max = 5e-5;			% Max Frequency application Delay
		delay_V_max = 2e-5;			% Max Voltage application Delay
		%
		F_discretization_step = 0.05;
		
		%% WL Parameters
		
		wl_prob = [0.0345 0.3448 0.3448 0.2414 0.0345];			% Probability for each WL	
		wl_min_exec_us = [510 100 100 480 520];					% minimal execution in us
		wl_mean_exec_us = [520e2 210e1 180e1 560e2 1.2e6/1e1];	% mean execution in us

		wl_uniq_min = [0.9 0.3 0.4 0.6 0.9];					% minimal percentage of presence when Primary Wl
		wl_comb_array = [0 1 0 0 0;								% Combination for each wl with each other:
						 1 0 1 1 1;								%	1 = they can be present together, 0 = they cannot
						 0 1 0 1 1;								%	Matrix has to be symmetric, with 0 on the diagonal.
						 0 1 1 0 1;
						 0 1 1 1 0];
		max_num_wl = 3;											% max number of contemporary wls (primary [1] + secondary). It may not be always respected.
		mem_wl = [0.95 0.15 0.05 0.20 0];							% memory boundness
		
	end
	
	properties(Dependent)
		
		%% System Structure
		ipl;						% Instruction power levels
		FV_levels;
		
		%% Matlab/Simulation
		Ns;							% Number of States
		Ni;							% Number of inputs
		Ni_nl;						% Number of inputs for the non-linear case
		Ni_c;						% Number of controllable inputs
		Ni_nc;						% Number of non-controllable inputs
		Nout;						% Number of "Outputs"
		
		quantum_instr;
		
		% Matlab Matrices
		Ad_nom;
		Bd_nom;
		Ad_true;
		Bd_true;
		Ad_mpc;
		Bd_mpc;
		Ad_ctrl;
		Bd_ctrl;
		Ad_obs;
		Bd_obs;
		
		%% System Parameters
		static_pw_coeff;
		polyFV;
		
		V_Max;
		V_min;
		F_Max;
		F_min;
		core_Max_power;
		core_min_power;		
	end
	
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

	
	%% Models
	methods
		function R = spreading_r_computation(obj, source_area, plate_area, plate_thickness, k_source, R_plate)
			
			src_r = (source_area / pi)^(1/2);
			plt_r = (plate_area / pi)^(1/2);
			tau = plate_thickness / plt_r;
			
			epsi = src_r/plt_r;
			h = 1 / (R_plate * plate_area);
			Biot = h * plt_r / k_source;
			
			delta_c = pi + 1/(pi^(1/2)*epsi);
			phi_c = ( tanh(delta_c*tau) + delta_c/Biot ) / ( 1 + delta_c/Biot*tanh(delta_c*tau) );			
			
			psi_avg = (1-epsi)^(3/2)*phi_c / 2;
			psi_max = (1-epsi)*phi_c/(pi^(1/2));
			
			R = psi_avg / k_source / (source_area^(1/2));
			
		end
		function [A, B] = lin_model_create(obj, d, th_model_ver, temp) %, pw, ceff, pw_levels)
			
			lNc = obj.Nc;
			lNh = obj.Nh;
			lNv = obj.Nv;
			
			debug = 0;
			
			air_pos = 0;
			mb_pos = 1;
			pcb_pos = 2;
			al_pos = 3;			
			
			% create Distance/position matrix
			% This matrix should be given in input

			core_cols = lNv;
			core_rows = lNh;

			extl_cols = 2;
			extr_cols = 2;

			extt_rows = 1;
			extb_rows = 1;

			cols = core_cols + extl_cols + extr_cols;
			rows = core_rows + extt_rows + extb_rows;

			%TODO instead of initial values, add air
			dRNmat = 0.05*ones(rows, cols); %zeros(rows, cols);
			dRWmat = 0.05*ones(rows, cols); %zeros(rows, cols);
			dRSmat = 0.05*ones(rows, cols); %zeros(rows, cols);
			dREmat = 0.05*ones(rows, cols); %zeros(rows, cols);
			
			lCNmat = 0.05*ones(rows, cols); %zeros(rows, cols);
			lCSmat = 0.05*ones(rows, cols); %zeros(rows, cols);
			lCWmat = 0.05*ones(rows, cols); %zeros(rows, cols);
			lCEmat = 0.05*ones(rows, cols); %zeros(rows, cols);
			
			lCcoreNmat = zeros(rows, cols); %zeros(rows, cols);
			lCcoreSmat = zeros(rows, cols); %zeros(rows, cols);
			lCcoreWmat = zeros(rows, cols); %zeros(rows, cols);
			lCcoreEmat = zeros(rows, cols); %zeros(rows, cols);
			
			k1mat = ones(rows, cols)*1e-16; %zeros(rows, cols);
			k2mat = ones(rows, cols)*0.025; %zeros(rows, cols);
			c1mat = zeros(rows, cols);
			c2mat = zeros(rows, cols);

			alpha_k_si = -4e-3;
			alpha_k_cu = -1e-4;
			k_si = 127; %148;
			k_cu = 398.5; %400;
			k_air = 0.025;
			k_al = 225.94;
			k_pcb = 3.096;
			k_mb = 0.167;
			alpha_c_si = 1e-3;
			alpha_c_cu = 3e-4;
			c_si = 1.7243e+06; %1.66e6;
			c_cu = 3.4794e+06; %3.44e6;
			c_air = 1004*1.29*1.5; %SHC*air density*"air compression factor" % SHC=Specific Heat Capacity
			c_al = 921*2698; %SHC*density
			c_pcb = 753*2900; %SHC*density
			c_mb = 795*1900; %SHC*density

			
			% Doing cores:
			for r = (1+extt_rows):(rows - extb_rows)
				for c = (1+extl_cols):(cols - extr_cols)
					dRNmat(r,c) = 2.30e-3;
					dRSmat(r,c) = 1e-3;

					dRWmat(r,c) = 0.75e-3;
					dREmat(r,c) = 0.75e-3;
					
					lCNmat(r,c) = 2.30e-3;
					lCSmat(r,c) = 1e-3;
					lCWmat(r,c) = 0.75e-3;
					lCEmat(r,c) = 0.75e-3;
					
					lCcoreNmat(r,c) = 1e-3;
					lCcoreSmat(r,c) = 1e-3;
					lCcoreWmat(r,c) = 0.75e-3;
					lCcoreEmat(r,c) = 0.75e-3;

					% Internal Mesh connection
					if (mod(c,2) == 1) && (c~=(core_cols+extl_cols))
						dREmat(r,c) = dREmat(r,c) + 0.25e-3;
					end
					
					% Eternal Mesh connection: row
					if r == (1+extt_rows)
						dRNmat(r,c) = dRNmat(r,c) + 0.3e-3;
					end
					if r == (rows - extb_rows)
						dRSmat(r,c) = dRSmat(r,c) + 0.3e-3;
					end
					% Eternal Mesh connection: cols
					if c == (1+extl_cols)
						dRWmat(r,c) = dRWmat(r,c) + 0.3e-3;
					end
					if c == (cols - extr_cols)
						dREmat(r,c) = dREmat(r,c) + 0.3e-3;
					end
					
					k1mat(r,c) = k_si;
					k2mat(r,c) = k_cu;
					c1mat(r,c) = c_si;
					c2mat(r,c) = c_cu;
				end	
			end
			
			% thickness
			t_si = 5e-4; %3e-4;
			t_cu = 1e-3; %1.1e-3;
			%TODO
			t_al = 1.5e-2;%5e-3;
			t_air = 10e-2;%5e-3;
			t_pcb = 1e-3;
			t_mb = 2e-3;
			air_factor = 100;
			air_time_factor = 10; %20			
			
			% Additional vertical resistance
			R_TIM1 = 1;%0.25; %6 4.5
			R_TIM2 = 2.5; %6; %4; 0.75;
			case_fan_dis = 0.25; %0.85 %0.1 - 1
			heatsink_fan_dis = 10;
			fins_coeff = 3;
			fan_nom_speed = 1000; 
			
			pw2therm_coeff = 1.05; %0.90 %0.41; 0.7;
			
			%TODO: this does not consider multi-chiplet and stuff
			chiplet_width = max(sum(dREmat(1+extt_rows:end-extb_rows,1+extl_cols:end-extr_cols),2) + sum(dRWmat(1+extt_rows:end-extb_rows,1+extl_cols:end-extr_cols),2));
			chiplet_length = max(sum(dRNmat(1+extt_rows:end-extb_rows,1+extl_cols:end-extr_cols),2) + sum(dRSmat(1+extt_rows:end-extb_rows,1+extl_cols:end-extr_cols),2));
	
			%TODO parametrize 1cm and 10cm additional dimension
			al_li = (chiplet_length+2*1e-2);
			al_wi = (chiplet_width+2*1e-2);
			air_li = 20e-2;
			air_wi = 20e-2;
			A_al = al_li*al_wi;
			A_air = air_li*air_wi;
			A_pcb = (chiplet_length + 1.5e-2)*(chiplet_width + 1.5e-2);
			A_mb = 10e-2*10e-2;
			
			C_pcb = c_pcb * A_pcb*t_pcb;
			C_mb = c_mb * A_mb*t_mb;
			C_air = c_air * A_air*t_air;
			C_al = c_al * A_al*t_al;
			
			%partial
			Ri_air_tot_v = t_air/A_air/k_air/air_factor;
			Ri_air_al_v = (t_air-t_al-t_cu-t_si-t_pcb)/A_air/k_air/air_factor;
			Ri_air_AL1 = (air_li-al_li)/2 /(air_wi*t_air) /k_air/air_factor;
			Ri_air_AL2 = (air_wi-al_wi)/2 /(air_li*t_air) /k_air/air_factor;
			%
			Ri_al_v = t_al/A_al/k_al;
			Ri_alAL1_v = al_li /(al_wi*t_al) /k_al;
			Ri_alAL2_v = al_wi /(al_li*t_al) /k_al;
			%
			Ri_pcb_v = t_pcb/A_pcb/k_pcb;
			Ri_mb_v = t_mb/(A_mb - A_al)/k_mb;
			%
			RaL = 1.948E+09 * (t_air-t_al-t_cu-t_si-t_pcb)^3;
			h_alair_cond = RaL^(1/4) * 0.54 * k_air / (t_air-t_al-t_cu-t_si-t_pcb);
			R_alair_cond = 1 / (h_alair_cond * A_al * fins_coeff);
			
			RaL_AL12 = 1.948E+09 * (t_air-t_al-t_cu-t_si-t_pcb + t_al/2)^3;
			h_alairAL12_cond = k_air/(t_air-t_al-t_cu-t_si-t_pcb + t_al/2) * ...
				(0.825 + (0.387 * RaL_AL12^(1/6))/(1+(0.492/7.039E-01)^(9/16))^(8/27))^2;
			R_alairAL1_cond = 1 / (h_alairAL12_cond * (al_wi*t_al));
			R_alairAL2_cond = 1 / (h_alairAL12_cond * (al_li*t_al));
			
			RaL = 1.948E+09 * (t_air)^3;
			h_mbair_cond = RaL^(1/4) * 0.54 * k_air / (t_air);
			R_mbair_cond = 1 / (h_mbair_cond * (A_mb - A_al));
			
			%source_area, plate_area, plate_thickness, k_source, R_plate
			R_spr_alair = obj.spreading_r_computation(A_al, A_air, t_air, k_al, Ri_air_al_v);
			R_alair_v = Ri_al_v + R_alair_cond + R_spr_alair;
			%
			R_spr_alair1 = obj.spreading_r_computation((al_wi*t_al), (air_wi*t_air), air_li, k_al, Ri_air_AL1);
			R_alairAL1_v = Ri_alAL1_v + R_alairAL1_cond + R_spr_alair1;
			
			R_spr_alair2 = obj.spreading_r_computation((al_li*t_al), (air_li*t_air), air_wi, k_al, Ri_air_AL2);
			R_alairAL2_v = Ri_alAL2_v + R_alairAL2_cond + R_spr_alair2;
			%
			R_spr_pcbmb = obj.spreading_r_computation(A_pcb, A_mb, t_mb, k_pcb, Ri_mb_v);
			R_pcbmb_v = Ri_pcb_v + Ri_mb_v + R_spr_pcbmb;
			R_spr_mbair = obj.spreading_r_computation(A_mb, A_air, t_air, k_mb, Ri_air_tot_v);
			R_mbair_v = Ri_mb_v + R_mbair_cond + R_spr_mbair;
			
			%TODO NOWWW
			%R_air_v = R_air_v * obj.ThMoNoise(fix(i/2)+1,7);
			%C_air = Ci_air * obj.ThMoNoise(fix(i/2)+1,7);
			si_pcb_fact = 0.1;
			pcb_mb_fact = 15;
			mb_air_fact = 10;

			A=zeros(obj.Ns);
			
			for i=1:2*lNc
				if (mod(i,2)>0)				
					% Core Position
					ci = fix((i-1)/2)+1;
					irow = fix((ci-1)/core_cols)+1 + extt_rows;
					icol = mod((ci-1),core_cols)+1 + extl_cols;

					% Computing Core R
					li = dRNmat(irow,icol) + dRSmat(irow,icol);
					wi = dREmat(irow,icol) + dRWmat(irow,icol);
					Ri_si_hn = li/(wi*t_si)/k1mat(irow,icol);
					Ri_si_hs = Ri_si_hn;
					Ri_si_he = wi/(li*t_si)/k1mat(irow,icol);
					Ri_si_hw = Ri_si_he;

					Ri_si_v = t_si/(wi*li)/k1mat(irow,icol);
					Ri_cu_v = t_cu/(wi*li)/k2mat(irow,icol);
					R_spr_cual = obj.spreading_r_computation((wi*li), A_al, t_al, k_cu, Ri_al_v);
					R_spr_sipcb = obj.spreading_r_computation((wi*li), A_pcb, t_pcb, k_si, Ri_pcb_v);

					Ri_cu_hn = li/(wi*t_cu)/k2mat(irow,icol);
					Ri_cu_hs = Ri_cu_hn;
					Ri_cu_he = wi/(li*t_cu)/k2mat(irow,icol);
					Ri_cu_hw = Ri_cu_he;

					% Computing Core C
					li = lCNmat(irow,icol) + lCSmat(irow,icol);
					wi = lCEmat(irow,icol) + lCWmat(irow,icol);
					Ci_si = c1mat(irow, icol) * li*wi*t_si;
					Ci_cu = c2mat(irow, icol) * li*wi*t_cu;
					%
					
					% Computing neighbourhood R
					li = dRNmat(irow-1,icol) + dRSmat(irow-1,icol);
					wi = dREmat(irow-1,icol) + dRWmat(irow-1,icol);
					Rn_si_h = li/(wi*t_si)/k1mat(irow-1,icol);
					Rn_cu_h = li/(wi*t_cu)/k2mat(irow-1,icol);
					li = dRNmat(irow+1,icol) + dRSmat(irow+1,icol);
					wi = dREmat(irow+1,icol) + dRWmat(irow+1,icol);
					Rs_si_h = li/(wi*t_si)/k1mat(irow+1,icol);
					Rs_cu_h = li/(wi*t_cu)/k2mat(irow+1,icol);
					li = dRNmat(irow,icol+1) + dRSmat(irow,icol+1);
					wi = dREmat(irow,icol+1) + dRWmat(irow,icol+1);
					Re_si_h = wi/(li*t_si)/k1mat(irow,icol+1);
					Re_cu_h = wi/(li*t_cu)/k2mat(irow,icol+1);
					li = dRNmat(irow,icol-1) + dRSmat(irow,icol-1);
					wi = dREmat(irow,icol-1) + dRWmat(irow,icol-1);
					Rw_si_h = wi/(li*t_si)/k1mat(irow,icol-1);
					Rw_cu_h = wi/(li*t_cu)/k2mat(irow,icol-1);

					% Series
					R_si_hn = Ri_si_hn + Rn_si_h;
					R_si_hs = Ri_si_hs + Rs_si_h;
					R_si_he = Ri_si_he + Re_si_h;
					R_si_hw = Ri_si_hw + Rw_si_h;

					R_sicu_v = Ri_si_v + Ri_cu_v + R_TIM1;

					R_cu_hn = Ri_cu_hn + Rn_cu_h;
					R_cu_hs = Ri_cu_hs + Rs_cu_h;
					R_cu_he = Ri_cu_he + Re_cu_h;
					R_cu_hw = Ri_cu_hw + Rw_cu_h;					
					
					R_cual_v = Ri_cu_v + R_TIM2 + Ri_al_v + R_spr_cual;
					R_sipcb_v = Ri_si_v + Ri_pcb_v + R_spr_sipcb;
					
					C_si = Ci_si;
					C_cu = Ci_cu;

					% Noise
					if d==1
					R_si_hn = R_si_hn * obj.ThMoNoise(fix(i/2)+1,1);
					R_si_hs = R_si_hs * obj.ThMoNoise(fix(i/2)+1,1);
					R_si_he = R_si_he * obj.ThMoNoise(fix(i/2)+1,1);
					R_si_hw = R_si_hw * obj.ThMoNoise(fix(i/2)+1,1);

					R_sicu_v = R_sicu_v * obj.ThMoNoise(fix(i/2)+1,3);
					R_cual_v = R_cual_v * obj.ThMoNoise(fix(i/2)+1,4);

					R_cu_hn = R_cu_hn * obj.ThMoNoise(fix(i/2)+1,2);
					R_cu_hs = R_cu_hs * obj.ThMoNoise(fix(i/2)+1,2);
					R_cu_he = R_cu_he * obj.ThMoNoise(fix(i/2)+1,2);
					R_cu_hw = R_cu_hw * obj.ThMoNoise(fix(i/2)+1,2);

					C_si = Ci_si * obj.ThMoNoise(fix(i/2)+1,5);
					C_cu = Ci_cu * obj.ThMoNoise(fix(i/2)+1,6);

					end
					
					switch th_model_ver
						case 10
							C_si = 1e-3; %C_si / 6 *1.1;
							C_cu = 1.2e-2; %C_cu / 2;

							R_si_hn = R_si_hn * 1.44;
							R_si_hs = R_si_hs * 1.44;
							R_si_he = R_si_he * 2.67;
							R_si_hw = R_si_hw * 2.67;

							R_sicu_v = 5; %R_si_v * 5;

							R_cu_hn = R_cu_hn * 0.56;
							R_cu_hs = R_cu_hs * 0.56;
							%R_cu_he
							%R_cu_hw

							R_cual_v = 29; %R_cu_v * 5 * 1.3;
						case 1	
							ci = fix((i-1)/2)+1;
							irow = fix((ci-1)/core_cols)+1;
							icol = mod((ci-1),core_cols)+1;
							
							pw2therm_coeff = 0.9;
							
							C_si = C_si*1.1;
							C_cu = C_cu*1.1;

							g1=obj.gaussian_filter(max(obj.Nv, obj.Nh),1.2)*16;
							g1 = g1 + (1-0.025 - g1(1));
							g2=obj.gaussian_filter(max(obj.Nv, obj.Nh),2.5)*16;
							%g2 = g2 + (1-0.025 - g2(1));
							g2 = 1./g2 * 0.75;						
							
							R_si_hn = R_si_hn * g2(irow,icol); % * 2.67;
							R_si_hs = R_si_hs * g2(irow,icol); % * 2.67;
							R_si_he = R_si_he * g2(irow,icol); % * 2.67;
							R_si_hw = R_si_hw * g2(irow,icol); % * 2.67;							

							R_sicu_v = R_sicu_v / 2.2; %1.5;

							R_cu_hn = R_cu_hn * g2(irow,icol) * 2; %* 4;
							R_cu_hs = R_cu_hs * g2(irow,icol) * 2; % * 4;
							R_cu_he = R_cu_he * g2(irow,icol) * 2; % * 4;
							R_cu_hw = R_cu_hw * g2(irow,icol) * 2; % * 4;

							R_cual_v = R_cual_v; %* 1.5;% * 1.3;							
							R_sipcb_v = R_sipcb_v; % * 1.5;
							R_pcbmb_v = R_pcbmb_v;
							
							%R_alair_v
							%R_alairAL1_v
							%R_alairAL2_v							
							
						case 2 
							ci = fix((i-1)/2)+1;
							irow = fix((ci-1)/core_cols)+1;
							icol = mod((ci-1),core_cols)+1;							
							
							C_si = C_si*1.1;
							C_cu = C_cu*1.1;
							
							%pw2therm_coeff = 0.9;

							g1=obj.gaussian_filter(max(obj.Nv, obj.Nh),1.2)*16;
							g1 = g1 + (1-0.025 - g1(1));
							g2=obj.gaussian_filter(max(obj.Nv, obj.Nh),2.5)*16;
							%g2 = g2 + (1-0.025 - g2(1));
							g2 = 1./g2 * 0.7;						
							
							R_si_hn = R_si_hn * g2(irow,icol); % * 2.67;
							R_si_hs = R_si_hs * g2(irow,icol); % * 2.67;
							R_si_he = R_si_he * g2(irow,icol); % * 2.67;
							R_si_hw = R_si_hw * g2(irow,icol); % * 2.67;							

							R_sicu_v = R_sicu_v / 1.2;

							R_cu_hn = R_cu_hn * g2(irow,icol) * 2; %* 4;
							R_cu_hs = R_cu_hs * g2(irow,icol) * 2; % * 4;
							R_cu_he = R_cu_he * g2(irow,icol) * 2; % * 4;
							R_cu_hw = R_cu_hw * g2(irow,icol) * 2; % * 4;

							R_cual_v = R_cual_v / 1.5;% * 1.3;							
							R_sipcb_v = R_sipcb_v / 1.5;
							R_pcbmb_v = R_pcbmb_v / 1.07;
							
							%R_alair_v
							%R_alairAL1_v
							%R_alairAL2_v
							
							air_therm_dev = 0.2; %exp(obj.Ts / (C_al*Ri_al_v) ) - 1; %0.4;
					end
					
					if debug
						clc;
						drow = fix((ci-1)/core_cols)+1;
						dcol = mod((ci-1),core_cols)+1;
						ci
						
						C_si
						%C_core
						dC(drow, dcol) = C_si;
						%dCcore(drow, dcol) = C_core;
						C_cu
						
						R_si_hn
						R_si_hs
						R_si_he
						R_si_hw
						
						dRN(drow, dcol) = R_si_hn;
						dRS(drow, dcol) = R_si_hs;
						dRE(drow, dcol) = R_si_he;
						dRW(drow, dcol) = R_si_hw;

						R_sicu_v
						
						dRV(drow, dcol) = R_sicu_v;

						R_cu_hn
						R_cu_hs
						R_cu_he
						R_cu_hw

						R_cual_v
					end
				end
				
				%for j=1:2*lNc
					
				% ==============
				% diagonal terms
				% ==============

				%Silicon (die temp dyn) 
				if (mod(i,2)>0)
					A(i,i) = -1/(R_sicu_v*C_si) -1/C_si * (1/R_si_he+1/R_si_hn+1/R_si_hw+1/R_si_hs) - 1/(C_si*R_sipcb_v) * si_pcb_fact;
					
					%PCB:
					A(i, end-pcb_pos) = 1/(C_si*R_sipcb_v) * si_pcb_fact;
					A(end-pcb_pos, i) = 1/(C_pcb*R_sipcb_v) * si_pcb_fact;
				% Heat-Spread (Cooper dyn)
				else %if (mod(i,2)==0)
					A(i,i) = -1/(R_sicu_v*C_cu) -1/(R_cual_v*C_cu) -1/C_cu * (1/R_cu_he+1/R_cu_hn+1/R_cu_hw+1/R_cu_hs);
					
					%heat sink:
					A(i, end-al_pos) = 1/(C_cu*R_cual_v);
					A(end-al_pos, i) = 1/(C_al*R_cual_v);
				end

				% ==============
				% Vertical coupling terms btween silicon and copper (flow from die to spreader)
				% ==============
				if (mod(i,2)>0)
					A(i,i+1)=1/(R_sicu_v*C_si);
				else %if (mod(i,2)==0)
					A(i,i-1)=1/(R_sicu_v*C_cu);
				end

				% ==============
				% horizontal coupling of die thermal dyn (silicon)
				% ==============
				%if(lNc~=1) %avoid single core case
				if (mod(i,2)>0)
					% vertici
					if(i==1)
					   A(i,i+2)=1/(R_si_he*C_si); 
					   A(i,i+2*lNv)=1/(R_si_hs*C_si);
					elseif(i==2*lNv-1)
					   A(i,i-2)=1/(R_si_hw*C_si); 
					   A(i,i+2*lNv)=1/(R_si_hs*C_si);
					elseif (i==(lNh-1)*2*lNv+1) 
					   A(i,i+2)=1/(R_si_he*C_si); 
					   A(i,i-2*lNv)=1/(R_si_hn*C_si);
					elseif (i==2*lNc-1) 
					   A(i,i-2)=1/(R_si_hw*C_si); 
					   A(i,i-2*lNv)=1/(R_si_hn*C_si);

					% bordi verticali (west and east)
					elseif ((mod(i,2*lNv)==1)&&(i~=(lNh-1)*2*lNv+1)&&(i~=1)) 
					   A(i,i+2)=1/(R_si_he*C_si); 
					   A(i,i+2*lNv)=1/(R_si_hs*C_si);
					   A(i,i-2*lNv)=1/(R_si_hn*C_si);
					elseif (((mod(i+1,2*lNv)==0) && (i~=2*lNc-1) && (i~=2*lNv-1))) 
					   A(i,i-2)=1/(R_si_hw*C_si); 
					   A(i,i+2*lNv)=1/(R_si_hs*C_si);
					   A(i,i-2*lNv)=1/(R_si_hn*C_si);

					% horizonatal boundaries (upper and lower)
					elseif ((i~=1)&& (i<2*lNv-2))
					   A(i,i-2)=1/(R_si_hw*C_si); 
					   A(i,i+2)=1/(R_si_he*C_si); 
					   A(i,i+2*lNv)=1/(R_si_hs*C_si);
					elseif ((i>(lNh-1)*2*lNv+1) && (i<2*lNc-1))
					   A(i,i-2)=1/(R_si_hw*C_si); 
					   A(i,i+2)=1/(R_si_he*C_si); 
					   A(i,i-2*lNv)=1/(R_si_hn*C_si);

					%internal
					else
					  A(i,i-2)=1/(R_si_hw*C_si); 
					  A(i,i+2)=1/(R_si_he*C_si); 
					  A(i,i-2*lNv)=1/(R_si_hn*C_si);
					  A(i,i+2*lNv)=1/(R_si_hs*C_si);           
					end
				end

				% ==============
				% horizontal coupling of sperader thermal dyn (copper)
				% ==============
				if (mod(i,2)==0)
					% vertexes
					if(i==2)
					   A(i,i+2)=1/(R_cu_he*C_cu); 
					   A(i,i+2*lNv)=1/(R_cu_hs*C_cu);
					elseif(i==2*lNv)
					   A(i,i-2)=1/(R_cu_hw*C_cu); 
					   A(i,i+2*lNv)=1/(R_cu_hs*C_cu);
					elseif (i==(lNh-1)*2*lNv+2) 
					   A(i,i+2)=1/(R_cu_he*C_cu); 
					   A(i,i-2*lNv)=1/(R_cu_hn*C_cu);
					elseif (i==2*lNc) 
					   A(i,i-2)=1/(R_cu_hw*C_cu); 
					   A(i,i-2*lNv)=1/(R_cu_hn*C_cu);

					% vertical boundaries (west and east)
					elseif ((mod(i,2*lNv)==2)&&(i~=(lNh-1)*2*lNv+2) &&(i~=2)) 
					   A(i,i+2)=1/(R_cu_he*C_cu); 
					   A(i,i+2*lNv)=1/(R_cu_hs*C_cu);
					   A(i,i-2*lNv)=1/(R_cu_hn*C_cu);
					elseif (((mod(i,2*lNv)==0) && (i~=2*lNc) && (i~=2*lNv)))
					   A(i,i-2)=1/(R_cu_hw*C_cu); 
					   A(i,i+2*lNv)=1/(R_cu_hs*C_cu);
					   A(i,i-2*lNv)=1/(R_cu_hn*C_cu);

					% horizontal boundaries (upper and lower)
					elseif((i~=2)&& (i<2*lNv))
					   A(i,i-2)=1/(R_cu_hw*C_cu); 
					   A(i,i+2)=1/(R_cu_he*C_cu); 
					   A(i,i+2*lNv)=1/(R_cu_hs*C_cu);
					elseif ((i>(lNh-1)*2*lNv+2) && (i<2*lNc))
					   A(i,i-2)=1/(R_cu_hw*C_cu); 
					   A(i,i+2)=1/(R_cu_he*C_cu); 
					   A(i,i-2*lNv)=1/(R_cu_hn*C_cu);

					%internal
					else
					  A(i,i-2)=1/(R_cu_hw*C_cu); 
					  A(i,i+2)=1/(R_cu_he*C_cu); 
					  A(i,i-2*lNv)=1/(R_cu_hn*C_cu);
					  A(i,i+2*lNv)=1/(R_cu_hs*C_cu);           
					end
				end
				%end %single core case
				
				% Last Columns Poles (airs)
				% HeatSpreader
				if mod(i,2)==0
					% diagonal terms for vertexes
					if (i==2) 
						A(i,end-air_pos)=1/C_cu * (1/R_cu_hn + 1/R_cu_hw);
						A(end-air_pos,i)=1/C_air * (1/R_cu_hn + 1/R_cu_hw);
					elseif (i==2*lNv)
						A(i,end-air_pos)=1/C_cu * (1/R_cu_hn + 1/R_cu_he);
						A(end-air_pos,i)=1/C_air * (1/R_cu_hn + 1/R_cu_he);
					elseif (i==(lNh-1)*2*lNv+2) 
						A(i,end-air_pos)=1/C_cu * (1/R_cu_hs + 1/R_cu_hw);
						A(end-air_pos,i)=1/C_air * (1/R_cu_hs + 1/R_cu_hw);
					elseif (i==2*lNc)
						A(i,end-air_pos)=1/C_cu * (1/R_cu_hs + 1/R_cu_he);
						A(end-air_pos,i)=1/C_air * (1/R_cu_hs + 1/R_cu_he);
						
					% diagonal terms for vertical (west and east) boundary nodes
					elseif ((mod(i,2*lNv)==2)&&(i~=(lNh-1)*2*lNv+2)&&(i~=2))
						A(i,end-air_pos)=1/(C_cu*R_cu_hw);
						A(end-air_pos,i)=1/(C_air*R_cu_hw);
					elseif ((mod(i,2*lNv)==0) && (i~=2*lNc) && (i~=2*lNv))
						A(i,end-air_pos)=1/(C_cu*R_cu_he);
						A(end-air_pos,i)=1/(C_air*R_cu_he);
					% diagonal terms for horizontal (upper and lower) boundary nodes
					elseif ((i~=2)&& (i<2*lNv))
						A(i,end-air_pos)=1/(C_cu*R_cu_hn);
						A(end-air_pos,i)=1/(C_air*R_cu_hn);
					elseif ((i>(lNh-1)*2*lNv+2) && (i<2*lNc))
						A(i,end-air_pos)=1/(C_cu*R_cu_hs);
						A(end-air_pos,i)=1/(C_air*R_cu_hs);
						
					% diagonal terms corresponding to other nodes (internal)
					%else
						%A(i,end)=0;
						%A(end,i)=0;
					end
					
				% Silicium
				else
					% diagonal terms for vertexes
					if (i==1)
						A(i,end-air_pos)=1/C_si * (1/R_si_hn + 1/R_si_hw);
						A(end-air_pos,i)=1/C_air * (1/R_si_hn + 1/R_si_hw);					
					elseif (i==2*lNv-1)
						A(i,end-air_pos)=1/C_si * (1/R_si_hn + 1/R_si_he);
						A(end-air_pos,i)=1/C_air * (1/R_si_hn + 1/R_si_he);	
					elseif (i==(lNh-1)*2*lNv+1)
						A(i,end-air_pos)=1/C_si * (1/R_si_hs + 1/R_si_hw);
						A(end-air_pos,i)=1/C_air * (1/R_si_hs + 1/R_si_hw);						
					elseif (i==2*lNc-1)
						A(i,end-air_pos)=1/C_si * (1/R_si_hs + 1/R_si_he);
						A(end-air_pos,i)=1/C_air * (1/R_si_hs + 1/R_si_he);	

					% diagonal terms for vertical (west and east) boundary nodes
					elseif (mod(i,2*lNv)==1)&&(i~=(lNh-1)*2*lNv+1)&&(i~=1)
						A(i,end-air_pos)=1/(C_si*R_si_hw);
						A(end-air_pos,i)=1/(C_air*R_si_hw);
					elseif (mod(i+1,2*lNv)==0) && (i~=2*lNc-1) && (i~=2*lNv-1)
						A(i,end-air_pos)=1/(C_si*R_si_he);
						A(end-air_pos,i)=1/(C_air*R_si_he);
					% diagonal terms for horizontal (upper and lower) boundary nodes
					elseif (i~=1)&& (i<2*lNv-2)
						A(i,end-air_pos)=1/(C_si*R_si_hn);
						A(end-air_pos,i)=1/(C_air*R_si_hn);
					elseif (i>(lNh-1)*2*lNv+1) && (i<2*lNc-1)
						A(i,end-air_pos)=1/(C_si*R_si_hs);
						A(end-air_pos,i)=1/(C_air*R_si_hs);

					%  diagonal terms corresponding to other nodes (internal)
					%else
						%B(i,end) = 0;
					end
				end
				
				 %%%%%%%%
				% CASE 2 %
				 %%%%%%%%
				% for performance, add it above
				% HeatSpreader part 2
				if mod(i,2)==0
					if (th_model_ver == 2)
						%for j=1:lNc
						j=i/2; %remember
						row = fix((j-1)/core_cols)+1;
						col = mod((j-1),core_rows)+1;

						%1: Add/change the final input
							A(2*j,end-air_pos) = A(2*j,end-air_pos)-1/(R_cual_v*C_cu) + ...
								(1/(R_cual_v*C_cu))*(1-air_therm_dev)^(row-1);
						%2: add previous rows part
						for k=1:row-1
							id2 = j - (core_cols*(row-k));
							A(2*j, 2*id2) = A(2*j, 2*id2) + ...
								air_therm_dev * (1-air_therm_dev)^(row-1-k) * ...
								1/(R_cual_v*C_cu);
						end
						%end
						
						%3: fix self
						%This is not needed beacuse all coefficients stay
						%the same!
						A(2*j, 2*j) = -(sum(A(2*j,:)) - A(2*j, 2*j));
						
					end	%(th_model_ver == 2)
				end %mod(i,2)==0

			end % for i

			%HeatSink
			A(end-al_pos, end-air_pos) = 1/C_al * (1/R_alair_v + 2/R_alairAL1_v + 2/R_alairAL2_v) + heatsink_fan_dis/C_al;
			A(end-air_pos, end-al_pos) = 1/C_air * (1/R_alair_v + 2/R_alairAL1_v + 2/R_alairAL2_v) + heatsink_fan_dis/C_al; %here C_al or C_air?
			A(end-al_pos, end-al_pos) = 0;
			A(end-al_pos, end-al_pos) = -sum(A(end-al_pos,:));
			
			%PCB
			%negligible with air
			A(end-pcb_pos, end-mb_pos) =1/(C_pcb*R_pcbmb_v) * pcb_mb_fact;
			A(end-mb_pos, end-pcb_pos) =1/(C_mb*R_pcbmb_v) * pcb_mb_fact;
			A(end-pcb_pos, end-pcb_pos) = 0;
			A(end-pcb_pos, end-pcb_pos) = -sum(A(end-pcb_pos,:));
			
			%Motherboard
			A(end-mb_pos, end-air_pos) = 1/(C_mb * R_mbair_v) * mb_air_fact;
			A(end-air_pos, end-mb_pos) = 1/(C_air * R_mbair_v) * mb_air_fact;
			A(end-mb_pos, end-mb_pos) = 0;
			A(end-mb_pos, end-mb_pos) = -sum(A(end-mb_pos,:));
			
			%Air
			%A(end, end) = -(lNc+(lNh-2)*2+(lNv-2)*2+4*2)/(R_cu_v*C_cu) -100*lNc/4; % -0.149075;
			A(end-air_pos, end-air_pos)=0;
			A(end-air_pos, end-air_pos) = -sum(A(end-air_pos,:)) - case_fan_dis; % -0.149075;

			% ==============
			% Make it symmetrical
			% ==============
			% CANNOT DO THIS ANYMORE:
			%	because now C_si is different between cores, making
			%	different Rs
			%	NOT doing THIS, means that the Resistance of core i will
			%	be different (due to noise) seen by different cores.
% 			for i=1:2*lNc
% 				for j=1:i-2
% 					if A(i,j)~=A(j,i)
% 						if rand(1) >= 0.5
% 							A(i,i) = A(i,i) + (A(i,j)-A(j,i));
% 							A(i,j) = A(j,i);
% 						else
% 							A(j,j) = A(j,j) + (A(j,i)-A(i,j));
% 							A(j,i) = A(i,j);
% 						end					
% 					end
% 				end
% 				for j=i+2:2*lNc
% 					if A(i,j)~=A(j,i)
% 						if rand(1) >= 0.5
% 							A(i,i) = A(i,i) + (A(i,j)-A(j,i));
% 							A(i,j) = A(j,i);
% 						else
% 							A(j,j) = A(j,j) + (A(j,i)-A(i,j));
% 							A(j,i) = A(i,j);
% 						end					
% 					end
% 				end
% 			end
			
			% ==============
			% input matrix
			% ==============

			% inputs: core power (controlled) ambient temperature (disturbance)
			B=zeros(obj.Ns,obj.Ni);
			k=0;
			for i=1:2*lNc
				if mod(i,2)>0					
					k=k+1;
					
					% Core Position
					ci = fix((i-1)/2)+1;
					irow = fix((ci-1)/core_cols)+1 + extt_rows;
					icol = mod((ci-1),core_cols)+1 + extl_cols;
					
					li = lCcoreNmat(irow,icol) + lCcoreSmat(irow,icol);
					wi = lCcoreEmat(irow,icol) + lCcoreWmat(irow,icol);
					C_core = c1mat(irow, icol) * li*wi*t_si;					
					if d==1
						C_core = C_core * obj.ThMoNoise(fix(i/2)+1,5);
					end
					
					if (th_model_ver == 1) %|| (th_model_ver == 2)
						g2=obj.gaussian_filter(max(obj.Nv, obj.Nh),3)*16;
						g2 = g2 + (1-0.01 - g2(1));
						C_core = C_core / g2(irow-extt_rows,icol-extl_cols);% / 1.5;
					end							
							
					B(i,k)=1/C_core * pw2therm_coeff;
				end
			end

			B(end-air_pos,:) = [zeros(1, lNc) case_fan_dis/fan_nom_speed];

		end %function
		function pu = power_compute(obj, F, V, T, d_i,d_p, noise) 
			
			if (nargin < 7) || isempty(noise)
				noise = obj.model_variation;
			end
			
			P = [obj.static_pw_coeff obj.ceff_pw_coeff] / 1000;

			if (obj.exp_leakage == 1)
				P = [P(1:size(obj.static_pw_coeff,2)) obj.exp_leak_coeff/1000 P(size(obj.static_pw_coeff,2)+1:end)];
			end

			%ot = u(Ni_f+Ni_v+1:end);

			Power_static = V*P(1) + d_p*P(3);
			ci_index = 4;
			if size(P,2) > ((ci_index-1)+obj.ipl)
				maxT = 125+273.15;
				Power_static = Power_static .* exp(V*P(4) + (min(T,ones(length(T),1)*maxT)-273.15)*P(5) - ones(length(V),1)*P(6));
				ci_index = ci_index+3;
			end

			ci = P(ci_index:ci_index+obj.ipl-1);
			
			tt = size(d_i,1) > 1;
			if tt && noise
				npw = obj.core_pw_noise_char;
			else
				npw = 1;
			end

			Ceff = d_i * ci' .* npw;
			Power_dyn = Ceff .* F .* (V .* V);

			pu = Power_static + Power_dyn;		
		end
		function [dxdt, power] = nl_model_dyn(obj,A,B,x,u,d_i,d_p, ot) %, pw, ceff, pw_levels)

			Ni_f = obj.Nc;
			Ni_v = obj.vd;
			
			dxdt = zeros(obj.Ns,1);

			F = u(1:Ni_f);
			V = obj.VDom*u(Ni_f+1:Ni_f+Ni_v);
			
			power = obj.power_compute(F,V,obj.C(1:obj.Nc,:)*x,d_i,d_p);
			dxdt = A*x + B*[power; ot];
		end
		function [sp] = wl_prob_compute(obj, primary_wl, new_wl)
			
			hm_min = obj.wl_uniq_min(primary_wl);
			
			% Compute How Much?
			hmp = rand(1) + hm_min*obj.dur_c;
			hmp = hmp + (hmp>1)*(1-hmp) + (hmp<hm_min)*(hm_min-hmp);
			% Compute prob secondary wls
			%sp = (rand(1,5)-obj.wl_prob).*(obj.wl_comb_array(primary_wl,:)|(new_wl>0));
			sp = rand(1,obj.ipl);
			ttsp = (sp>(1-obj.wl_prob));
			if sum(ttsp) > 0
				sp = ttsp .* sp;
			else
				sp = sp .* obj.wl_prob;
			end	
			sp = sp.*(obj.wl_comb_array(primary_wl,:)|(new_wl>0));
			sp = sp - (sp<=0).*sp + (sp<=0.001)*0.001;
			sp = sp + new_wl;
			% Remove excessive secondary wls
			while sum(sp>0) > obj.max_num_wl-1
				spmin = min(sp(sp>0&sp<1));
				if isempty(spmin)
					break;
				end
				sp = sp - sp.*(sp==spmin);
			end
			% Compute how much? secondary wls
			sp = sp - new_wl;
			% Fix case sp=0 and new_wl != 0
			if sum(sp) == 0
				sp = obj.wl_comb_array(primary_wl,:);
			end
			% Normalize on 1-hmp, and forzing hmp to be <1 in case new_wl != 0
			if sum(new_wl)>0
				max_hmp = max(hm_min, 1 - obj.min_secondary_space*min(obj.max_num_wl,sum(obj.wl_comb_array(primary_wl,:)))); 
				hmp = hmp - (hmp>max_hmp)*(hmp-(max_hmp));
			end
			sp = sp/sum(sp)*(1-hmp);
			
			sp(primary_wl) = sp(primary_wl) + hmp;
			
		end
		function [obj] = init_compute_model(obj, A, B)
			obj.qt_storage = obj.quantum_instr*ones(obj.Nc, 1);
			obj.wl_index = ones(obj.Nc, 1);
			obj.F_s = obj.F_min*ones(obj.Nc,1);
			obj.V_s = obj.V_min*ones(obj.vd,1);
			obj.A_s = A;
			obj.B_s = B;
			% cannot put it here, if I want it to be constant among several
			% runs
			%obj.core_pw_noise_char = ones(obj.Nc,1);
			%if (obj.model_variation)
			%obj = obj.create_core_pw_noise();
			%obj = obj.create_thermal_model_noise();
			%end
		end
		function [uplot, xplot, d_is, pw_ms, obj] = compute_model(obj, N, x, V, F, d_p)
			
			d_is = zeros(1,obj.ipl);
			delay_F_index = zeros(obj.Nc,1);
			delay_V_index = zeros(obj.vd,1);		
			pw_ms = 0;
			
			delay_F_div = ceil(obj.delay_F_mean / obj.Ts);
			delay_V_div = ceil(obj.delay_V_mean / obj.Ts);
			
			uplot = zeros(N, obj.Ni_c);
			xplot = zeros(N, obj.Ns);
			
			for sim=1:N
				%index = (ix0-1)*N + sim;

				% Delay Application
				% Application of V before F
				delay_F_index = delay_F_index + obj.VDom*(V==obj.V_s);
				delay_V_index = delay_V_index + (V~=obj.V_s);
				tv = delay_V_index > (delay_V_div*ones(obj.vd,1));
				obj.V_s = obj.V_s - obj.V_s.*tv + V.*tv;
				%
				tf = delay_F_index > (delay_F_div*ones(obj.Nc,1));
				obj.F_s = obj.F_s - obj.F_s.*tf + F.*tf;
				%noise 

				%wl
				wlp = zeros(obj.Nc, obj.ipl);
				wlpn = zeros(obj.Nc, obj.ipl);
				for c=1:obj.Nc
					wlp(c,:) = obj.wrplot(c,:,max(mod(obj.wl_index(c), size(obj.wrplot,3)),1));
					wlpn(c,:) = obj.wrplot(c,:,max(mod(obj.wl_index(c)+1, size(obj.wrplot,3)),1+1));
				end
				%if (obj.F_s(1) > obj.F_min)
				%	asd = 1;
				%end
				instr = obj.F_s * 1e9 * (obj.Ts);
				mem_instr = obj.F_min * 1e9 * (obj.Ts);
				a_citr = sum(wlp.*obj.mem_wl,2);
				a_citrn = sum(wlpn.*obj.mem_wl,2);
				citr = a_citr.*mem_instr + (1-a_citr).*instr;
				citrn = a_citrn.*mem_instr + (1-a_citrn).*instr;
				pwl = obj.qt_storage./citr;
				ttwl = (pwl>1);
				pwl = pwl + ttwl.*(1-pwl);

				wl = pwl .* wlp + ...
					(1-pwl) .* wlpn;

				obj.qt_storage = obj.quantum_instr*(1-ttwl) - (citr - obj.qt_storage);
				%if sum((obj.qt_storage <0)>0)
				%	stop=1;
				%end
				obj.wl_index = obj.wl_index + (1-ttwl);

				d_is = d_is + wl;

				pu_s = obj.power_compute(obj.F_s,obj.VDom*obj.V_s,obj.C(1:obj.Nc,:)*x,wl,d_p);
				pw_ms = pw_ms + sum(pu_s);
				%x = A_nom*x + B_nom*[u;temp_amb*1000];
				x = obj.A_s*x + obj.B_s*[pu_s;obj.temp_amb*1000];		

				uplot(sim,:) = pu_s;
				xplot(sim,:) = x;
			end	%ssim
			d_is = d_is / N;
			pw_ms = pw_ms / N;
		end
	end %methods
	
	%% Controllers
	methods
		function obj = init_mpc(obj, Nhzn, Q, R, xref, uref, yref)
			
			obj.xref = zeros(obj.Ns, 1);
			obj.uref = zeros(obj.Ni_c, 1);
			obj.yref = zeros(obj.Nout, 1);
			obj.Q = zeros(obj.Ns);
			obj.R = zeros(obj.Ni_c);
			obj.R2 = zeros(obj.Ni_c);
			%obj.Nhzn = 3;
			obj.Ctx = zeros(obj.Nhzn, obj.Ns);
			obj.Ctu = zeros(obj.Nhzn, obj.Ni_c);
			obj.Cty = zeros(obj.Nhzn, obj.Nout);
			
			obj.mpc_robustness = 0.5;
			
			obj.umin = -inf; %obj.core_min_power;
			obj.uMax = inf; %obj.core_Max_power;
			obj.xmin = -inf; %273.15-50;
			obj.xMax = inf; %obj.core_crit_temp;
			obj.ymin = -inf; %273.15-50;
			obj.yMax = inf; %obj.core_crit_temp;
			obj.usum = ones(1,obj.vd+1)*inf;
			
			if nargin >= 2
				obj.Nhzn = Nhzn;				
			end
			if nargin >= 3
				obj.Q = Q;
			end
			if nargin >= 4
				obj.R = R;
			end
			if nargin >= 5
				obj.xref = xref;
			end
			if nargin >= 6
				obj.uref = uref;
			end
			if nargin >= 7
				obj.yref = yref;
			end
			
		end %init_mpc		
		function obj = lin_mpc_setup(obj)
			
			% It's good practice to start by clearing YALMIPs internal database 
			% Every time you call sdpvar etc, an internal database grows larger
			yalmip('clear')

			x = sdpvar(repmat(obj.Ns,1,obj.Nhzn+1),repmat(1,1,obj.Nhzn+1));
			u = sdpvar(repmat(obj.Ni_c,1,obj.Nhzn),repmat(1,1,obj.Nhzn));
			
			sdpvar ot; %TODO: decide if single or with horizon
			%sdpvar ly_xref;
			ly_uref = sdpvar(obj.Ni_c,1);
			%sdpvar ly_yref;
			ly_usum = sdpvar(obj.Nhzn,obj.vd+1);
			%ly_Cty = dpvar(repmat(obj.Nout,1,obj.Nhzn),repmat(1,1,obj.Nhzn));

			constraints = [];
			
			Adl_mpc = obj.Ad_mpc;
			Bdl_mpc = obj.Bd_mpc;
			
			Bu = Bdl_mpc(:,1:obj.Ni_c);
			Bd = Bdl_mpc(:,obj.Ni_c+1:end);

			for k = 1:obj.Nhzn
					constraints = [constraints, x{k+1} == Adl_mpc*x{k}+Bu*u{k}+Bd*ot];
					
					%TODO: Ct I did only for y, and only for max
					if (~isinf(obj.umin))
						constraints = [constraints, u{k} >= obj.umin];
					end
					if (~isinf(obj.uMax))
						constraints = [constraints, u{k} <= obj.uMax - obj.Ctu(k,:)'];
					end
					if (~isinf(obj.xmin))
						constraints = [constraints, x{k+1} >= obj.xmin];
					end
					if (~isinf(obj.xMax))
						constraints = [constraints, x{k+1} <= obj.xMax];
					end
					if (~isinf(obj.ymin))
						constraints = [constraints, obj.C*x{k+1} >= obj.ymin];
					end
					if (~isinf(obj.yMax))
						constraints = [constraints, obj.C*x{k+1} <= obj.yMax - obj.Cty(k,:)'];
					end
					%if (~isinf(obj.usum(1)))
						constraints = [constraints,sum(u{k}) <= ly_usum(k,1) - sum(obj.Ctu(k,:),2)];
					%end
					%if (~isinf(obj.usum(2:end)))
					%
						constraints = [constraints, obj.VDom'*u{k} <= ly_usum(k,2:end)' - (obj.Ctu(k,:)*obj.VDom)'];
					%
					%end
			end
			
			objective = 0;

			% TODO: do Q and x (with xref), and "Q" and yref
			for k = 1:obj.Nhzn
				objective = objective + (u{k}-ly_uref)'*obj.R*(u{k}-ly_uref) + u{k}'*obj.R2*(u{k}) + x{k}'*obj.Q*x{k}; % + (x{k}-ly_xref)'*obj.Q*(x{k}-lx_uref);
			end
			
			ops = sdpsettings('verbose',1,'solver','quadprog', 'usex0',1);
			ops = sdpsettings('verbose',1,'solver','osqp', 'usex0',0); %You have specified an initial point, but the selected solver (OSQP) does not support warm-starts through YALMIP
			ops.quadprog.TolCon = 1e-2;
			ops.quadprog.TolX = 1e-5;
			ops.quadprog.TolFun = 1e-3;
			%ops.quadprog.MaxPCGIter = max(1, ops.quadprog.MaxPCGIter * 3);
			ops.quadprog.MaxIter = 50;

            % save yalmip model for potential extraction to external solver or code generation
            obj.ylmp_opt_variables = {x{1},ot,ly_uref,ly_usum};
            obj.ylmp_opt_output = u{1};
            obj.ylmp_constraints = constraints;
            obj.ylmp_objective = objective;
            % Controller
            obj.Controller = optimizer(constraints,objective,ops,obj.ylmp_opt_variables,obj.ylmp_opt_output);
			
		end %lin_mpc_setup	
		function [k0, k1, k2] = createPsApprox(obj)
			%creating the asdsad
			Fint = [obj.F_min obj.F_Max];
			Tint = [20 90];

			%Fint = [2.8 obj.F_Max];
			%Tint = [60 90];

			a = Tint(1);
			b = Tint(2);
			c = Fint(1);
			d = Fint(2);

			Ke = -obj.exp_leak_coeff(3)/1000;
			KT = obj.exp_leak_coeff(2)/1000;
			KV = obj.exp_leak_coeff(1)/1000;
			Ks = obj.leak_process / 1000;
			Kw = obj.leak_vdd / 1000;

			V0 = 0.9;
			%V0 = 1.2;
			KV1 = exp(KV*V0)*(1-KV*V0);
			KV2 = KV / (1 - KV*V0);

			alp = 0.2995; %0.15;

			c1 = exp(Ke)*KV1*(exp(KT*b) - exp(KT*a)) / (KT*(d-c)*(b-a));
			c2 = (d^(alp+1) - c^(alp+1))*(Ks*KV2 + Kw) / (alp+1);
			c3 = KV2*Kw*(d^(2*alp+1) - c^(2*alp+1)) / (2*alp+1);
			c4 = Ks*(d-c);

			m0 = c1 * (c2 + c3 + c4);

			c1 = 6*exp(Ke)*KV1*( KT*(b-a)*(exp(KT*b) + exp(KT*a)) - 2*(exp(KT*b) - exp(KT*a)) ) / (KT^2*(d-c)*(b-a)^3);
			c2 = (d^(alp+1) - c^(alp+1))*(Ks*KV2 + Kw) / (alp+1);
			c3 = KV2*Kw*(d^(2*alp+1) - c^(2*alp+1)) / (2*alp+1);
			c4 = Ks*(d-c);

			m1 = c1 * (c2 + c3 + c4);

			c1 = 2*exp(Ke)*KV1*(exp(KT*b) - exp(KT*a)) / (KT*(b-a));
			cd1 = (d-c)*(d+c)^2 - 2*(d^(2*alp+2) - c^(2*alp+2))*(d+c)*(alp+1)^(-1) + 4*(d^(4*alp+3) - c^(4*alp+3))*(4*alp+3)^(-1);
			c2 = (d^(alp+1) - c^(alp+1))*(Ks*(d^(alp+1)+c^(alp+1))-(d+c)*(Ks*KV2+Kw))/(alp+1);
			c3 = KV2*Kw*(d^(2*alp+1)-c^(2*alp+1))*((d^(2*alp+1)+c^(2*alp+1)) - (d+c)) / (2*alp+1);
			c4 = 2*(d^(3*alp+2)-c^(3*alp+2))*(Ks*KV2 + Kw) / (3*alp+2);
			c5 = Ks*(d^2 - c^2);

			m2 = c1/cd1 * (c2 + c3 + c4 - c5);

			k1 = m1;
			k2 = m2;
			k0 = m0 - (m1*(b+a) + m2*(d+c))/2;			
		end
		function obj = lin_mpc_setup2(obj)
			
			% It's good practice to start by clearing YALMIPs internal database 
			% Every time you call sdpvar etc, an internal database grows larger
			yalmip('clear')

			x = sdpvar(repmat(obj.Ns,1,obj.Nhzn+1),repmat(1,1,obj.Nhzn+1));
			u = sdpvar(repmat(obj.Ni_c,1,obj.Nhzn),repmat(1,1,obj.Nhzn));
			
			sdpvar ot; %TODO: decide if single or with horizon
			if obj.separate_omega
				w = sdpvar(obj.Ni_c,1); %TODO: decide if single or with horizon
			end
			h0v = sdpvar(obj.Ni_c,1);
			%sdpvar ly_xref;
			ly_uref = sdpvar(obj.Ni_c,1);
			%sdpvar ly_yref;
			ly_usum = sdpvar(obj.Nhzn,obj.vd+1);
			%ly_Cty = dpvar(repmat(obj.Nout,1,obj.Nhzn),repmat(1,1,obj.Nhzn));

			Adl_mpc = obj.Ad_mpc;
			Bdl_mpc = obj.Bd_mpc;
			
			Bu = Bdl_mpc(:,1:obj.Ni_c);
			Bd = Bdl_mpc(:,obj.Ni_c+1:end);
			
			[~, h1, h2] = obj.createPsApprox();
			
			
			C2 = eye(obj.Ns);
			C2([2:2:end]) = 0;
			C2([end-4:end], :) = 0; %TODO parametrize
			
			C3 = obj.C(1:obj.Nc,:)';
			
			constraints = [];

			for k = 1:obj.Nhzn
					if obj.separate_omega
						constraints = [constraints, x{k+1} == (Adl_mpc + C2*h1)*x{k}-(C2*h1*273.15*ones(obj.Ns,1))+(Bu+(obj.C(1:obj.Nc,:)'*h2))*(u{k}.*w)+Bd*ot + C3*h0v];
						%[constraints, x{k+1} == (Adl_mpc + C2*h1)*x{k}-(C2*h1*273.15*ones(obj.Ns,1))+(Bu+h2)*(u{k}.*w)+Bd*ot + h0];

						%TODO: Ct I did only for y, and only for max
						if (~isinf(obj.umin))
							%constraints = [constraints, u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k}-273.15)) + h0 >= obj.umin];
						end
						if (~isinf(obj.uMax))
							%constraints = [constraints, u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k}-273.15)) + h0 <= obj.uMax - obj.Ctu(k,:)'];
						end
						if (~isinf(obj.xmin))
							%constraints = [constraints, x{k+1} >= obj.xmin];
						end
						if (~isinf(obj.xMax))
							%constraints = [constraints, x{k+1} <= obj.xMax];
						end
						if (~isinf(obj.ymin))
							%constraints = [constraints, obj.C)*x{k+1} >= obj.ymin];
						end
						if (~isinf(obj.yMax))
							constraints = [constraints, obj.C*x{k+1} <= obj.yMax - obj.Cty(k,:)'];
						end
						%if (~isinf(obj.usum(1)))
							constraints = [constraints,sum(u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k+1}-273.15)) + h0v) <= ly_usum(k,1) - sum(obj.Ctu(k,:),2)];
							%constraints = [constraints,sum(u{k}*(h2+1) + h1*(273.15) + h0) <= ly_usum(k,1) - sum(obj.Ctu(k,:),2)];
						%end
						%if (~isinf(obj.usum(2:end)))
						%
							%constraints = [constraints, obj.VDom'*(u{k}.*w*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k}-273.15)) + h0) <= ly_usum(k,2:end)' - (obj.Ctu(k,:)*obj.VDom)'];
						%
						%end
						
					else
						constraints = [constraints, x{k+1} == (Adl_mpc + C2*h1)*x{k}-(C2*h1*273.15*ones(obj.Ns,1))+(Bu+(obj.C(1:obj.Nc,:)'*h2))*(u{k})+Bd*ot + C3*h0v];
						%[constraints, x{k+1} == (Adl_mpc + C2*h1)*x{k}-(C2*h1*273.15*ones(obj.Ns,1))+(Bu+h2)*(u{k}.*w)+Bd*ot + h0];

						%TODO: Ct I did only for y, and only for max
						if (~isinf(obj.umin))
							%constraints = [constraints, u{k}*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k}-273.15)) + h0 >= obj.umin];
						end
						if (~isinf(obj.uMax))
							%constraints = [constraints, u{k}*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k}-273.15)) + h0 <= obj.uMax - obj.Ctu(k,:)'];
						end
						if (~isinf(obj.xmin))
							%constraints = [constraints, x{k+1} >= obj.xmin];
						end
						if (~isinf(obj.xMax))
							%constraints = [constraints, x{k+1} <= obj.xMax];
						end
						if (~isinf(obj.ymin))
							%constraints = [constraints, obj.C)*x{k+1} >= obj.ymin];
						end
						if (~isinf(obj.yMax))
							constraints = [constraints, obj.C*x{k+1} <= obj.yMax - obj.Cty(k,:)'];
						end
						%if (~isinf(obj.usum(1)))
							constraints = [constraints,sum(u{k}*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k+1}-273.15)) + h0v) <= ly_usum(k,1) - sum(obj.Ctu(k,:),2)];
							%constraints = [constraints,sum(u{k}*(h2+1) + h1*(273.15) + h0) <= ly_usum(k,1) - sum(obj.Ctu(k,:),2)];
						%end
						%if (~isinf(obj.usum(2:end)))
						%
							%constraints = [constraints, obj.VDom'*(u{k}*(h2+1) + h1*(obj.C(1:obj.Nc,:)*(x{k}-273.15)) + h0) <= ly_usum(k,2:end)' - (obj.Ctu(k,:)*obj.VDom)'];
						%
						%end
					end
			end
			
			objective = 0;

			% TODO: do Q and x (with xref), and "Q" and yref
			for k = 1:obj.Nhzn
				objective = objective + (u{k}-ly_uref)'*obj.R*(u{k}-ly_uref) + u{k}'*obj.R2*(u{k}) + x{k}'*obj.Q*x{k}; % + (x{k}-ly_xref)'*obj.Q*(x{k}-lx_uref);
			end
			
			ops = sdpsettings('verbose',1,'solver','quadprog', 'usex0',1);
			ops = sdpsettings('verbose',1,'solver','osqp', 'usex0',0); %You have specified an initial point, but the selected solver (OSQP) does not support warm-starts through YALMIP
			ops.quadprog.TolCon = 1e-2;
			ops.quadprog.TolX = 1e-5;
			ops.quadprog.TolFun = 1e-3;
			ops.convertconvexquad = 0;
			%ops.quadprog.MaxPCGIter = max(1, ops.quadprog.MaxPCGIter * 3);
			ops.quadprog.MaxIter = 50;
			if obj.separate_omega
				obj.Controller = optimizer(constraints,objective,ops,{x{1},w,ot, h0v, ly_uref,ly_usum},{u{1}, x{2}});
			else
				obj.Controller = optimizer(constraints,objective,ops,{x{1},ot, h0v, ly_uref,ly_usum},{u{1}, x{2}});
			end
			obj.Controller
			
		end %lin_mpc_setup
		function uout = lin_mpc2(obj, currentx, w, ot, h0v, uref, usum)
			
			if nargin < 4
				error("Needs current x, w and d");
			end
			if nargin < 6
				usum = obj.usum;
			end
			if nargin < 5
				uref = obj.uref;
			end
			
			if size(usum,1) < obj.Nhzn
				usum(size(usum,1)+1:obj.Nhzn,:) = repmat(usum(size(usum,1),:), obj.Nhzn - size(usum,1),1);
			end
			
			if obj.separate_omega
				[uout, problem,~,~,optimizer_object] = obj.Controller({currentx, w, ot, h0v, uref, usum});
			else
				[uout, problem,~,~,optimizer_object] = obj.Controller({currentx, ot, h0v, uref, usum});
			end

			% Analyze error flags
			if problem
				warning(yalmiperror(problem))
			end
			
		end %lin_mpc
		function uout = lin_mpc(obj, currentx, ot, uref, usum)
			
			if nargin < 3
				error("Needs current x and d");
			end
			if nargin < 5
				usum = obj.usum;
			end
			if nargin < 4
				uref = obj.uref;
			end
			
			if size(usum,1) < obj.Nhzn
				usum(size(usum,1)+1:obj.Nhzn,:) = repmat(usum(size(usum,1),:), obj.Nhzn - size(usum,1),1);
			end

			[uout, problem,~,~,optimizer_object] = obj.Controller({currentx, ot, uref, usum});

			% Analyze error flags
			if problem
				warning(yalmiperror(problem))
			end
			
		end %lin_mpc		
		function obj = nl_mpc_setup(obj)
			
			% It's good practice to start by clearing YALMIPs internal database 
			% Every time you call sdpvar etc, an internal database grows larger
			yalmip('clear')
			
			x = sdpvar(repmat(obj.Ns,1,obj.Nhzn+1),repmat(1,1,obj.Nhzn+1));
			u = sdpvar(repmat(obj.Ni_nl,1,obj.Nhzn),repmat(1,1,obj.Nhzn));
			
			d_i = sdpvar(obj.Nc,obj.ipl,obj.Nhzn); %(repmat(obj.Nc,obj.ipl,obj.Nhzn),repmat(obj.ipl,1,obj.Nhzn));
			d_p = ones(obj.Nc,1); %TODO
			sdpvar ot;
			pu = sdpvar(repmat(obj.Nc,1,obj.Nhzn+1),repmat(1,1,obj.Nhzn+1));
			ly_usum = sdpvar(obj.Nhzn,obj.vd+1);
			ly_uref = sdpvar(obj.Ni_nl,1);
			
			%temp = binvar(obj.FV_levels,1); %obj.vd,obj.Nhzn);
			%tempFV = sdpvar(obj.FV_levels,1);
			%tempFV1 = double2sdpvar(obj.FV_table(:,1));
			%tempFV3 = double2sdpvar(obj.FV_table(:,3));
			
			obj.polyFV_opt = obj.polyFV;
			
			Adl_mpc = obj.Ad_mpc;
			Bdl_mpc = obj.Bd_mpc;
			
			constraints = [];
			for k = 1:obj.Nhzn
				[mx, pu{k}] = obj.nl_model_dyn(Adl_mpc,Bdl_mpc,x{k},u{k},d_i(:,:,k),d_p, ot);
					constraints = [constraints, x{k+1} == mx];
					
					%TODO: Ct I did only for y, and only for max
					if (~isinf(obj.umin))
						constraints = [constraints, u{k} >= obj.umin];
					end
					if (~isinf(obj.uMax))
						constraints = [constraints, u{k} <= obj.uMax];
					end
					if (~isinf(obj.xmin))
						constraints = [constraints, x{k+1} >= obj.xmin];
					end
					if (~isinf(obj.xMax))
						constraints = [constraints, x{k+1} <= obj.xMax];
					end
					if (~isinf(obj.ymin))
						constraints = [constraints, obj.C*x{k+1} >= obj.ymin];
					end
					if (~isinf(obj.yMax))
						constraints = [constraints, obj.C*x{k+1} <= obj.yMax - obj.Cty(k,:)'];
					end
					%if (~isinf(obj.usum(1)))
						constraints = [constraints,sum(u{k}) <= ly_usum(k,1) - sum(obj.Ctu(k,:),2)];
					%end
					%if (~isinf(obj.usum(2:end)))
						constraints = [constraints, obj.VDom'*u{k} <= ly_usum(k,2:end)' - (obj.Ctu(k,:)*obj.VDom)'];
					%end
					
					% u{k}(1:obj.Nc) > ones(obj.Nc, obj.FV_levels)*diag(obj.FV_table(:,3)) :compare the Frequencies to the FV table
					% then I take the sum(,2) and add 1 to identify the V_max indexes
					% I apply those indexes to obj.FV_table to find V_max
					% I diagonalize it and multiply for VDom to check for domains
					% I take the max of each and transpose for verticality
					% [zeros(obj.Nc,1) inf] to fix the >= problem of F_Max

					%constraints = [constraints, u{k}(obj.Nc+1:end) >= max(diag(obj.FV_table(sum(repmat(u{k}(1:obj.Nc),1, obj.FV_levels) >= ones(obj.Nc, obj.FV_levels)*diag(obj.FV_table(:,3)),2) + ones(obj.Nc,1),1))*obj.VDom)'];
					%[constraints, u{k}(obj.Nc+1:end) >= max(diag(obj.FV_table(,1))*obj.VDom)'];
					
					Model = [];
					
					%Model = [ implies( temp(:,:,k), ones(obj.FV_levels,obj.vd)*diag(max(diag(u{k}(1:obj.Nc,1))*obj.VDom)) >= repmat(obj.FV_table(:,3) + [zeros(obj.FV_levels-1,1); inf],1,obj.vd) ), ...
					%	implies( 1-temp(:,:,k), ones(obj.FV_levels,obj.vd)*diag(max(diag(u{k}(1:obj.Nc,1))*obj.VDom)) <= repmat(obj.FV_table(:,3) + [zeros(obj.FV_levels-1,1); inf], 1,obj.vd) ), ...
						%u{k}(obj.Nc+1:end,1) >= tempFV1(sum(temp(:,:,k))+ones(1, obj.vd)) %, tempFV == obj.FV_levels(:,1)
					 %];
					
					% WORKING! (compiled) but not working in optimizer!
					%for v = 1:obj.vd
					%	Model = [implies(temp, ones(obj.FV_levels,1)*max(diag(u{k}(1:obj.Nc,1))*obj.VDom(:,v)) >= obj.FV_table(:,3) + [zeros(obj.FV_levels-1,1); obj.F_Max*5]), ...
					%		implies(1-temp, ones(obj.FV_levels,1)*max(diag(u{k}(1:obj.Nc,1))*obj.VDom(:,v)) <= obj.FV_table(:,3) + [zeros(obj.FV_levels-1,1); obj.F_Max*5]), ...
					%		u{k}(obj.Nc+v,1) >= tempFV1(sum(temp)+1)];
						
					%	constraints = [constraints, Model ];
					%end					
					% LINEARIZE
					l = max(diag(u{k}(1:obj.Nc,1))*obj.VDom);
					constraints = [constraints, u{k}(obj.Nc+1:end,1) >= (obj.polyFV_opt*[l.^3;l.^2;l;ones(1,obj.vd)])'];
					
					
					%Model = [implies(temp, max(diag(u{k}(1:obj.Nc,1))*obj.VDom) >= (obj.FV_table(:,3) + [zeros(obj.FV_levels-1,1); inf])), 
					%	implies(1-temp, max(diag(u{k}(1:obj.Nc,1))*obj.VDom) <= (obj.FV_table(:,3) + [zeros(obj.FV_levels-1,1); inf])),
					%	u{k}(obj.Nc+1:end,1) >= obj.FV_table(sum(temp)+ones(1, obj.vd),1) ];
					
					%for v = 1:obj.vd
					%	for sc = 1:obj.FV_levels
							% error: sum() is not working
							%constraints = [constraints, implies( sum(max(diag(u{k}(1:obj.Nc))*obj.VDom(:,v)) >= (obj.FV_table(:,3) + [zeros(obj.FV_levels-1,1); inf] )+1 == sc, u{k}(obj.Nc+v) <= obj.FV_table(sc, 1))];
							%constraints = [constraints, implies( sum(max(diag(u{k}(1:obj.Nc))*obj.VDom(:,v)) >= (obj.FV_table(:,3) + [zeros(obj.FV_levels-1,1); inf] )+1) == sc, u{k}(obj.Nc+v) <= obj.FV_table(sc, 1))];
							%constraints = [constraints]; %, implies( max((max(diag(u{k}(1:obj.Nc))*obj.VDom(:,v)) >= (obj.FV_table(:,3) + [zeros(obj.FV_levels-1,1); inf] ) ) .* (1:1:obj.FV_levels)') +1 == sc, u{k}(obj.Nc+v) >= obj.FV_table(sc, 1))];
					%	end
					%end

					
			end
			
			objective = 0;

			% TODO: do Q and x (with xref), and "Q" and yref
			for k = 1:obj.Nhzn
				objective = objective + (u{k}-ly_uref)'*obj.R*(u{k}-ly_uref); % + (x{k}-ly_xref)'*obj.Q*(x{k}-lx_uref);
			end
			
			ops = sdpsettings('verbose',0,'solver','ipopt', 'usex0',1); %, 'allownonconvex',0
			ops.ipopt.tol = 1e-3;
			obj.Controller = optimizer(constraints,objective,ops,{x{1},d_i,ot,ly_uref,ly_usum},u{1});
			
		end %nl_mpc_setup
		function uout = nl_mpc(obj, currentx, d_i, ot, uref, usum)
			
			if nargin < 4
				error("Needs current x and d");
			end
			if nargin < 6
				usum = obj.usum;
			end
			if nargin < 5
				uref = obj.uref;
			end
			
			if size(usum,1) < obj.Nhzn
				usum(size(usum,1)+1:obj.Nhzn,:) = usum(size(usum,1),:);
			end

			[uout, problem,~,~,optimizer_object] = obj.Controller({currentx, d_i, ot, uref, usum});

			% Analyze error flags
			if problem
				warning(yalmiperror(problem))
			end
			
		end %nl_mpc	
		function [operc, omax] = max_uncertainty(obj, x, perc)
			if (nargin < 2)
				perc = 100;
			end
			omax = max(abs(x),[],"all");
			operc = prctile(abs(x),perc,"all");
		end
		function [Ct] = find_unc_set(obj, perc, show)
			
			if (nargin < 2) || isempty(perc)
				perc = 100;
			end
			if (nargin < 3) || isempty(show)
				show = 0;
			end
			
			air_pos = 0;
			mb_pos = 1;
			pcb_pos = 2;
			al_pos = 3;
			
			%w_max = obj.max_uncertainty(obj.zrplot,perc);
			%w_max = umax - umm;	
			
			copper_fix_temp = 10;
			
			% Compute forward reachable sets
			Nlim = obj.Nhzn;
			% max u
			%umax = [obj.core_Max_power*ones(obj.Ni_c,1); obj.temp_amb*1000];
			x = obj.core_crit_temp*ones(obj.Ns,1); 
			% fix heatspreader
			x(2:2:obj.Nc*2) = x(2:2:obj.Nc*2) - copper_fix_temp;
			%fix others:
			x(end-air_pos) = x(end-air_pos) - copper_fix_temp*4;
			x(end-mb_pos) = x(end-mb_pos) - copper_fix_temp*3;
			x(end-pcb_pos) = x(end-pcb_pos) - copper_fix_temp*1.5;
			x(end-al_pos) = x(end-al_pos) - copper_fix_temp*2;			
			
			disc = c2d(ss(obj.Ac_true, obj.Bc_true, obj.C, obj.D), obj.Ts_mpc);
			Adl_true2 = disc.A; %eye(Ns) + Ts_mpc*Ac_true;
			Bdl_true2 = disc.B; %Ts_mpc*Bc_true;
			Ct = zeros(Nlim, obj.Ns);
			Ctn = zeros(Nlim, obj.Ns);
			for s=1:Nlim
				umax = obj.power_compute(obj.F_Max, obj.V_Max, obj.C(1:obj.Nc,:)*x, [zeros(1, obj.ipl-1) 1], 1);
				x = Adl_true2*x+Bdl_true2*[umax; obj.temp_amb*(1000+1)];
				Ct(s,:) = x;
			end
			
			x = obj.core_crit_temp*ones(obj.Ns,1); 
			% fix heatspreader
			x(2:2:obj.Nc*2) = x(2:2:obj.Nc*2) - copper_fix_temp;
			%fix others:
			x(end-air_pos) = x(end-air_pos) - copper_fix_temp*4;
			x(end-mb_pos) = x(end-mb_pos) - copper_fix_temp*3;
			x(end-pcb_pos) = x(end-pcb_pos) - copper_fix_temp*1.5;
			x(end-al_pos) = x(end-al_pos) - copper_fix_temp*2;	
			
			Adl_mpc = obj.Ad_mpc;
			Bdl_mpc = obj.Bd_mpc;
			for s=1:Nlim
				umm = obj.power_compute(obj.F_Max, obj.V_Max, obj.C(1:obj.Nc,:)*x, [1 zeros(1, obj.ipl-1)], 1);
				x = Adl_mpc*x+Bdl_mpc*[umm; obj.temp_amb*(1000)];
				Ctn(s,:) = x;
			end
			
			Ct = Ct-Ctn;

			%show
			if show
				figure("Name", "MPC Robust Constraint");
				ax1 = subplot(1,2,1);
				plot(Ct*obj.C(1:obj.Nc,:)');
				grid on,xlabel('Steps'),ylabel('Cores State Difference [°C]'),hold on;
				ax1 = subplot(1,2,2);
				plot(Ct(:,2:2:end-max([air_pos,al_pos,mb_pos,pcb_pos])));
				grid on,xlabel('Steps'),ylabel('Copper State Difference [°C]'),hold on;
			end
			
			%saturate
			Ct = Ct + (Ct<0).*(-Ct);
			
			%output:
			Ct = Ct * (perc/100);
		end
		function [voltage_choice] = cp_voltage_choice(obj, Ft)
			%TODO issue when dimension of domains are not equal!!!
			for v=1:obj.vd
				extrV = sum( Ft(:,v) > (obj.FV_table(:,3) + [zeros(obj.FV_levels-1,1); inf])', 2);
				vote_cast(:,v) = extrV(nonzeros(obj.VDom(:,v).*[1:obj.Nc]')) + 1;
			end
			if size(vote_cast,1) == 1
				%problem: matrix become array and prctile does not work anymore
				vote_cast(2,:) = vote_cast;
			end
			voltage_choice = obj.FV_table(round(prctile(vote_cast,obj.voltage_rule)),1);
		end
		function [upid, ointeg] = cp_pid(obj, T, pid_target, pid_integ, pu)
			
			kp_l = obj.kp;
			ki_l = obj.ki*obj.Ts_ctrl;

			% PID error
			e = pid_target - T;

			% PID error banding
			eA = ~(obj.pid_e_down>0).*(e>0) + ~(obj.pid_e_up>0).*(e<=0);
			if (obj.pid_e_down>0)
				eA = (e*(1-obj.pid_e_band_coeff)/obj.pid_e_down + obj.pid_e_band_coeff) .* (e>0);
			end
			if (obj.pid_e_up>0)
				eA = eA + ((-e*(1-obj.pid_e_band_coeff)/obj.pid_e_up + obj.pid_e_band_coeff) .* (e<=0));
			end
			%eA = eA*obj.pid_e_band_coeff + ~eA;	
			eA = eA + (eA>1).*(1-eA);
			e = e .* eA;

			% PID Integral
			ointeg = pid_integ + ki_l*e;
			eI = ointeg>obj.aw_up;
			ointeg = eI*obj.aw_up + (~eI).*ointeg;
			eI = ointeg<(pu*obj.aw_down_c);
			ointeg = eI.*(pu*obj.aw_down_c) + (~eI).*ointeg;

			% PID Proporitonal (& Output)
			upid = kp_l*e + ointeg;

			% PID Output Saturation
			eU = upid>obj.sat_up;
			upid = eU*obj.sat_up + (~eU).*upid;
			eU = upid < (-(pu-ones(obj.Nc,1).*obj.min_pw_red*0.7));
			upid = eU.*(-(pu-ones(obj.Nc,1).*obj.min_pw_red*0.7)) + (~eU).*upid;
		end
		function [pu, pw_storage] = cp_pw_dispatcher(obj, T, pid_target, delta_p, ipu)
			
			safety_pw = -1e-3;
			
			if obj.dummy_pw
				perc = delta_p / sum(ipu);
				red = ipu*perc;
				pu = ipu - red;
				
				%check on pw
				%since the difference is negligible (0.7369 vs 0.8503), but also when
				%	leakage, I will just use MAX_CEFF.
				%min_pw_red = computed above
				tt = (pu < obj.min_pw_red*0.7);
				pu = pu-(pu.*tt) + obj.min_pw_red*0.7.*tt;
				pw_storage = delta_p - sum(pu);
			elseif obj.pw_inverse
				%compute alpha
				%check number
				safety = 1;
				tt = (pid_target - T - safety)>=0;
				tt_value = obj.core_crit_temp - pid_target;
				if tt_value < safety
					tt_value = safety1; %0.1; %TODO (considering that a margin of 10°C is a max -otherwise too much performance loss-)
				end
				pu = ipu;
				pw_storage = -delta_p;
				while(pw_storage < safety_pw) %-5/obj.vd)
					%1
					Alpha = (obj.core_crit_temp - T).*tt + tt_value.*(~tt);
					
					%2 
					% this seems worse
					%Alpha = 1./(pid_target - T).*tt + 0.99.*((tt - 1) * (-1));

					%test on numbers .....

					%normalize
					aa = sum(Alpha);
					Alpha = Alpha/aa;	

					%test on numbers .....

					%check on power
					%moved below application
					%delta_power*Alpha > pu-obj.core_min_power
					%alpha = (pu-obj.core_min_power)/delta_p	

					%apply
					pu = pu + Alpha*pw_storage;
					%check on pw
					%since the difference is negligible (0.7369 vs 0.8503), but also when
					%	leakage, I will just use MAX_CEFF.
					%min_pw_red = computed above
					tt2 = (pu < obj.min_pw_red*0.7);
					pu = pu-(pu.*tt2) + obj.min_pw_red*0.7.*tt2;
					pw_storage = (sum(ipu) - sum(pu)) - delta_p;
					if sum(tt2)==length(ipu)
						break;
					end
				end
			else
				%compute alpha
				%check number
				safety = 5;
				tt = (pid_target - T - safety)>=0;
				tt_value = obj.core_crit_temp - pid_target;
				tt2 = zeros(length(ipu),1);
				if tt_value < safety
					tt_value = safety; %0.1; %TODO (considering that a margin of 10°C is a max -otherwise too much performance loss-)
				end
				pu = ipu;
				pw_storage = -delta_p;
				while(pw_storage < safety_pw) %-5/obj.vd)
					%1
					Alpha = 1./(obj.core_crit_temp - T).*tt + 1./tt_value.*(~tt);

					%2 
					% this seems worse
					%Alpha = 1./(pid_target - T).*tt + 0.99.*((tt - 1) * (-1));

					%test on numbers .....

					%normalize
					Alpha = Alpha.*(~tt2);
					aa = sum(Alpha);
					Alpha = Alpha/aa;	

					%test on numbers .....

					%check on power
					%moved below application
					%delta_power*Alpha > pu-obj.core_min_power
					%alpha = (pu-obj.core_min_power)/delta_p	

					%apply
					pu = pu + Alpha*pw_storage;
					%check on pw
					%since the difference is negligible (0.7369 vs 0.8503), but also when
					%	leakage, I will just use MAX_CEFF.
					%min_pw_red = computed above
					tt2 = (pu < obj.min_pw_red*0.7);
					pu = pu-(pu.*tt2) + obj.min_pw_red*0.7.*tt2;
					pw_storage = (sum(ipu) - sum(pu)) - delta_p;
					if sum(tt2)==length(ipu)
						break;
					end
				end
			end		
		
		end
		function [pu, pw_storage] = cp_pw_dispatcher_c(obj, Ceff, delta_p, ipu)
			
				safety_pw = -1e-3;
				
				%perc = delta_p / sum(ipu);
				%red = ipu*perc;
				%pu = ipu - red;
				
				pu = ipu;
				tt = zeros(length(ipu),1);
				%pw_storage = (sum(ipu) - sum(pu)) - delta_p;
				pw_storage = -delta_p;
				%here it should be pw_storage = - delta_p; CHECK!
				while(pw_storage < safety_pw) %-5/obj.vd)
					Alpha = 1./(Ceff+(sum(Ceff)/length(Ceff)));
					Alpha = Alpha.*(~tt);
					%normalize
					aa = sum(Alpha);
					Alpha = Alpha/aa;
					%apply
					red = Alpha*pw_storage;
					red = red + (-red>Alpha.*pu).*(-Alpha.*pu - red);
					%pu = pu + Alpha*pw_storage;
					pu = pu + red;

					%check on pw
					%since the difference is negligible (0.7369 vs 0.8503), but also when
					%	leakage, I will just use MAX_CEFF.
					%min_pw_red = computed above
					tt = (pu <= obj.min_pw_red);
					pu = pu-(pu.*tt) + obj.min_pw_red.*tt;
					pw_storage = (sum(ipu) - sum(pu)) - delta_p;
					if sum(tt)==length(ipu)
						break;
					end
				end				
		end
		function [wl, Ceff] = cp_wl_process(obj, prev_wl, iwl)
			wl = prev_wl*(1-obj.alpha_wl) + iwl*obj.alpha_wl;
			Ceff = wl * (obj.ceff_pw_coeff / 1000)';
		end
		function [pw_adapt] = cp_pw_adapt(obj, prev_pw_adapt, pw_m, pw_est, pbc, a1, a2)
			if nargin < 6
				alpha_up = 0.04;
			else
				alpha_up = a1;
			end
			if nargin < 7
				alpha_down = 0.04; %0.08
			else
				alpha_down = a2;
			end

			delta = pw_m-pw_est;
			
			alpha_pw = (delta>0)*alpha_down + (delta<=0)*alpha_up;
			
			%TODO: Here I need to do the thing that connect the alpha with
			%a poly function when power budget change. atm power budget is
			%constant
			
			pw_adapt = prev_pw_adapt*(1-alpha_pw) + delta*alpha_pw;	
			
			%here? above?
			if pbc
				pw_adapt = delta;
			end
		end
	end %methods
	
	%% Matlab Sim
	
	% ==============
	% Plots
	% ==============
	methods (Static)
		[p,n] = numSubplots(n)
		g=gaussian_filter(Filter_size, sigma)
	end
		
	methods %(Static)		
		function [fig] = xutplot(obj,x,u)
			add_state = obj.Ns-obj.Nc*2; %TODO
			t = [0:1:length(x)-1]*obj.Ts;
			fig = figure('Name', 'Metrics'); %, 'Position', get(0, 'Screensize'));
			movegui(fig, 'northwest');
			ax1=subplot(4,1,1);
			plot(t, x(:,1+end-add_state:end)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Temp Rack [°C]'), hold on;
			legend("Heat-sink", "PCB", "Motherboard", "Air");
			%yline(obj.core_crit_temp-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			ax2=subplot(4,1,2);
			plot(t, x(:,2:2:end-add_state)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Copper [°C]'),hold on;
			if max(max(x(:,2:2:end-add_state))) >= (obj.core_crit_temp-(0.1*(obj.core_crit_temp-273)))
				yline(obj.core_crit_temp-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			end
			ax3=subplot(4,1,3);
			plot(t, x(:,1:2:end-add_state)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Silicon [°C]'), hold on;
			if max(max(x(:,1:2:end-1))) >= (obj.core_crit_temp-(0.1*(obj.core_crit_temp-273)))
				yline(obj.core_crit_temp-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			end
			
			ax4=subplot(4,1,4);
			%fig2 = figure('Name', 'Power', 'Position', get(0, 'Screensize'));
			plot(t, u);
			grid on,xlabel('Time [s]'),ylabel('Input Power [W]'), hold on;
			if max(max(u)) >= (obj.core_Max_power-(0.34*obj.core_Max_power))
				yline(obj.core_Max_power, '--', 'Max Power', 'LineWidth',1,  'Color', 'b');
			end
			if min(min(u)) <= (obj.core_min_power+(0.34*obj.core_Max_power))
				yline(obj.core_min_power, '--','min Power', 'LineWidth',1,  'Color', 'b');
			end
			
			linkaxes([ax1,ax2,ax3,ax4],'x');
		end
		function [fig] = tempplot(obj,t,x)
			fig = figure('Name', 'Temperature'); %, 'Position', get(0, 'Screensize'));
			plot(t, x-273.15);
			grid on,xlabel('Time [s]'),ylabel('State [°C]'),hold on;
			if max(max(x)) >= (obj.core_crit_temp-(0.1*(obj.core_crit_temp-273)))
				yline(obj.core_crit_temp-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			end
		end
		function [fig1, fig2] = powerplot(obj,t1,u)
			fig1 = figure('Name', 'Power');
			ax1=subplot(3,1,1);
			plot(t1, u);
			grid on,xlabel('Time [s]'),ylabel('Input Power [W]'), hold on;
			if max(max(u)) >= (obj.core_Max_power-(0.34*obj.core_Max_power))
				yline(obj.core_Max_power, '--', 'Max Power', 'LineWidth',1,  'Color', 'b');
			end
			if min(min(u)) <= (obj.core_min_power+(0.34*obj.core_Max_power))
				yline(obj.core_min_power, '--','min Power', 'LineWidth',1,  'Color', 'b');
			end
			
			ax2=subplot(3,1,2);
			dimf = t1(end)/obj.Ts_input;
			plot(obj.Ts_input*[0:dimf]', obj.urplot);
			grid on,xlabel('Time [s]'),ylabel('Input Ref [W]'), hold on;
			if max(max(u)) >= (obj.core_Max_power-(0.34*obj.core_Max_power))
				yline(obj.core_Max_power, '--', 'Max Power', 'LineWidth',1,  'Color', 'b');
			end
			if min(min(u)) <= (obj.core_min_power+(0.34*obj.core_Max_power))
				yline(obj.core_min_power, '--','min Power', 'LineWidth',1,  'Color', 'b');
			end
			
			ax3=subplot(3,1,3);
			plot(t1, obj.zrplot);
			grid on,xlabel('Time [s]'),ylabel('Input Noise [W]'), hold on;
			
			linkaxes([ax1,ax2,ax3],'x');
			
			fig2 = figure('Name', 'Domains Power Consumption');
			title('Domains Power Consumption')
			ax1 = subplot(2,1,1);
			gr = sum(u,2);
			plot(t1,gr);
			grid on,title('Total Power'),ylabel('[W]'),xlabel('Time [s]'),hold on;
			%if max(gr) >= obj.usum(1)*0.2
			plot([0:length(obj.tot_pw_budget)]*(obj.tsim/length(obj.tot_pw_budget)),[obj.tot_pw_budget(1); obj.tot_pw_budget], '--', 'LineWidth',1,  'Color', '#CC0000');
			ax2 = subplot(2,1,2);
			p = plot(t1,u*obj.VDom);
			grid on,title('Domain Power'),ylabel('[W]'),xlabel('Time [s]'),hold on;
			for pl=1:obj.vd
				plot([0:size(obj.quad_pw_budget,1)]*(obj.tsim/size(obj.quad_pw_budget,1)), [obj.quad_pw_budget(1,:); obj.quad_pw_budget], '--', 'LineWidth',1,  'Color', p(pl).Color);
			end
			%TODO, mettere gli if a yline
			
			linkaxes([ax1,ax2],'x');
		end
		function [fig] = powerconstrplot(obj,u1,u2)
						
			fig = figure('Name', 'Power Constraints Compliance');
			movegui(fig, 'southeast');
			title('Power Constraints analysis')
			
			t1 = [0:1:length(u1)-1]*obj.Ts;
			
			ax1 = subplot(3,4,[1:4]);
			gr = sum(u1,2);
			plot(t1,gr, '.');
			grid on,title('Total Power'),ylabel('[W]'),xlabel('Time [s]'),hold on;
			if ((nargin >=4) && (isempty(u2)==0))
				gr = sum(u2,2);
				plot(t1,gr);
			end
			%if max(gr) >= obj.usum(1)*0.2
			plot([0:max(length(obj.tot_pw_budget)-1,1)]*(obj.tsim/max(length(obj.tot_pw_budget)-1,1)),[obj.tot_pw_budget], '--', 'LineWidth',1,  'Color', '#CC0000');
			ax2 = subplot(3,4,[5:8]);
			p = plot(t1,u1*obj.VDom, '.');
			grid on,title('Domain Power'),ylabel('[W]'),xlabel('Time [s]'),hold on;
			if ((nargin >=4) && (isempty(u2)==0))
				p = plot(t1,u2*obj.VDom);
			end
			for pl=1:obj.vd
				plot([0:max(size(obj.quad_pw_budget,1)-1,1)]*(obj.tsim/max(size(obj.quad_pw_budget,1)-1,1)), [obj.quad_pw_budget], '--', 'LineWidth',1,  'Color', p(pl).Color);
			end
			
			u1 = u1(2:end,:);	
			dim = size(u1,1);
			%dim2 = dim
			
			pbt = repelem(obj.tot_pw_budget(min(2, length(obj.tot_pw_budget)):end),max(round(size(u1,1)/max(length(obj.tot_pw_budget)-1,1)),1),1);
			pbq = repelem(obj.quad_pw_budget(min(2, size(obj.quad_pw_budget,1)):end,:),max(round(size(u1,1)/max(size(obj.quad_pw_budget,1)-1,1)),1),1 );
			gr1 = [sum(sum(u1,2) > pbt)*obj.Ts; ...
					sum( (u1*obj.VDom) > pbq )'*obj.Ts];
			
			ttqt = [sum( (sum(u1,2) > pbt)>0 ); sum( ((u1*obj.VDom) > pbq )>0 )'];
			ttqt(ttqt<1) = 1;
			gp1 = [ sum((sum(u1,2) - pbt).*(sum(u1,2) > pbt)) / ttqt(1); ...
				sum( ((u1*obj.VDom) - pbq).*((u1*obj.VDom) > pbq))' ./  ttqt(2:end)];
			
			if (nargin < 4) || isempty(u2)
				BAR1 = [gr1(1) 0];
				BAR2 = [0 gp1(1)];
			else
				%TODO:
				%u2 = u2(2:end,:);
				%gr2 = [sum(sum(u2,2) > obj.tot_pw_budget)*obj.Ts; sum( (u2*obj.VDom)' > obj.quad_pw_budget, 2 )*obj.Ts];
				%BAR1 = [gr1(1) gr2(1) 0 0];
				%BAR2 = [0 0 gp1(1) gp2(1)];
			end
			
			subplot(3,4,9);
			%yyaxis left
			%axl = gca; % current axes
			b = bar(1,BAR1);
			%axl.Color = 'none';
			grid on,xlabel('Total'),ylabel('Time [s]');

			for pb = 1:2:floor(length(b)/2)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData/(dim*obj.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end			
			
			yyaxis right
			b2 = bar(1, BAR2);
			grid on,ylabel('\DeltaPw [W]');
			%axr = gca;
			%bar(BAR/(dim*obj.Ts)*100, 'FaceColor', "#0072BD");
			%axr.Color = 'none';
			%axr.YTick = axl.YTick/(dim*obj.Ts)*100;
			%ylabel('[%]');

			
			subplot(3,4,[10:12]);
			if (nargin < 4) || isempty(u2)
				BAR1 = [gr1(2:end) zeros(length(gr1)-1,1)];
				BAR2 = [zeros(length(gp1)-1,1) gp1(2:end)];
			else
				%TODO:
				%BAR = [gr1(2:end) gr2(2:end)];
			end
			%yyaxis left
			%axl = gca; % current axes
			b = bar(1:obj.vd, BAR1);
			%axl.Color = 'none';
			grid on,xlabel('Domain'),ylabel('Time [s]');
			
			for pb = 1:2:floor(length(b)/2)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData/(dim*obj.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			
			yyaxis right
			b2 = bar(1:obj.vd, BAR2);
			grid on,ylabel('\DeltaPw [W]');
			
			%axr = gca;
			%bar(BAR/(dim*obj.Ts)*100, 'FaceColor', "#0072BD");
			%axr.Color = 'none';
			%axr.YTick = axl.YTick/(dim*obj.Ts)*100;
			%ylabel('[%]');

		end
		function [fig] = tempconstrplot(obj,x1,x2)
			
			x1 = x1(2:end,:);
			dim = size(x1,1);
			%dim2 = dim
			gr1=obj.C(1:obj.Nc,:)*sum(x1 > obj.core_crit_temp)';
			tt1 = sum(x1 > obj.core_crit_temp);
			tt1(tt1==0) = 1;
			gt1=obj.C(1:obj.Nc,:)*(sum( (x1 - obj.core_crit_temp) .* (x1 > obj.core_crit_temp)) ./ tt1 )';

			if (nargin < 3) || isempty(x2)
				BAR1 = [sum(obj.Ts*gr1)/obj.Nc 0];
				BAR2 = [0 sum(gt1)/sum(tt1>1)];
			else
				x2 = x2(2:end,:);
				gr2=obj.C(1:obj.Nc,:)*sum(x2 > obj.core_crit_temp)';
				tt2 = sum(x2 > obj.core_crit_temp);
				tt2(tt2==0) = 1;
				gt2=obj.C(1:obj.Nc,:)*(sum( (x2 - obj.core_crit_temp) .* (x2 > obj.core_crit_temp)) ./ tt2 )';
				BAR1 = [sum(obj.Ts*gr1)/obj.Nc sum(obj.Ts*gr2)/obj.Nc 0 0 ];
				BAR2 = [0 0 sum(gt1)/sum(tt1>1) sum(gt2)/sum(tt2>1)];
			end			
			
			fig = figure('Name', 'Temperature Constraints Compliance');
			movegui(fig, 'southwest');
			title('Temperature Constraints analysis');
			subplot(5,4,[1,5,9])			
			%yyaxis left
			b = bar(1,BAR1);
			ylabel('Time [s]'),xlabel('Average'),grid on;
			for pb = 1:2:floor(length(b)/2)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData/(dim*obj.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			yyaxis right
			b2 = bar(1,BAR2);
			yl = ylabel('\DeltaT_{crit} [°C]');grid on;
			%yl.Position(2) = 0;
			%BAR = [(sum(gr1)/obj.Nc/(dim)*100)' (sum(gr2)/obj.Nc/(dim)*100)'];
			%bar(BAR,'FaceColor', "#0072BD");
			%ylabel('[%]');

			subplot(5,4,[2:4, 6:8, 10:12])
			if (nargin < 3) || isempty(x2)
				BAR1 = [obj.Ts*gr1 zeros(obj.Nc,1)];
				BAR2 = [zeros(obj.Nc,1) gt1];
			else
				BAR1 = [obj.Ts*gr1 obj.Ts*gr2 zeros(obj.Nc,1) zeros(obj.Nc,1)];
				BAR2 = [zeros(obj.Nc,1) zeros(obj.Nc,1) gt1 gt2];
			end
			%yyaxis left
			b = bar(1:obj.Nc, BAR1);
			ylabel('Time [s]'),xlabel('Core'),grid on;
			%{
			for pb = 1:2:floor(length(b)/2)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData/(dim*obj.Ts)*100,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			%}
			yyaxis right
			b2 = bar(1:obj.Nc, BAR2);
			ylabel('\DeltaT_{crit} [°C]'),grid on;
			%BAR = [(gr1/(dim)*100)' (gr2/(dim)*100)'];
			%bar(BAR,'FaceColor', "#0072BD");
			%ylabel('[%]');
			%total_temp_violation = (obj.Ts*sum(gr))
			
			subplot(5,4,[13:20])
			if (nargin < 3) || isempty(x2)
				BAR1 = [obj.C(1:obj.Nc,:)*max(x1,[],1)' - 273.15, zeros(obj.Nc,1)];
				BAR2 = [zeros(obj.Nc,1), obj.C(1:obj.Nc,:)*mean(x1)' - 273.15];
				ymlj = min(x1, [], "all");
			else
				BAR = [obj.C(1:obj.Nc,:)*max(x1,[],1)' obj.C(1:obj.Nc,:)*max(x2,[],1)'] - 273.15;
				ymlj = min(min(x1, [], "all"),  min(x2, [], "all"));
			end
			
			b = bar(1:obj.Nc, BAR1);
			yl = ylim;
			ylim([min(min(obj.x_init), ymlj)-273.15, yl(2)]);
			ylabel('T_{Max} [°C]'),xlabel('Core'),grid on;
			hold on;
			plot(xlim, [obj.core_crit_temp obj.core_crit_temp]-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			
			ylp = ylim;
			
			yyaxis right
			b2 = bar(1:obj.Nc, BAR2);
			ylabel('T_{Average} [°C]'),grid on;
			ylim(ylp);

		end
		function [fig] = obsplot(obj, xl, x)
			smxl = round((size(x,1)-1)/(size(xl,1)-1));
			smx = round((size(xl,1)-1)/(size(x,1)-1));
			
			gr = repelem(x(2:end,:),max(smx,1),1) - repelem(xl(2:end,:),max(smxl,1),1);
			
			fig = figure('Name', 'Observer Performance');
			%title('Domains Power Consumption')
			ax1 = subplot(3,1,1);
			plot(obj.tsim*[1/size(gr,1):(1/size(gr,1)):1]', gr(:,end));
			grid on,xlabel('Time [s]'),ylabel('Temp Rack [°C]'), hold on;
			ax2 = subplot(3,1,2);			
			plot(obj.tsim*[1/size(gr,1):(1/size(gr,1)):1]',gr(:,2:2:end-1));
			grid on,title('Copper State'),xlabel('Time [s]'),ylabel('Difference [°C]'),hold on;
			ax3 = subplot(3,1,3);
			plot(obj.tsim*[1/size(gr,1):(1/size(gr,1)):1]',gr(:,1:2:end-1));
			grid on,title('Silicon State'),xlabel('Time [s]'),ylabel('Difference [°C]'),hold on;
			
			linkaxes([ax1,ax2,ax3],'x');
		end
		function [fig] = fvplot(obj, t, f, v)
			
			fig = figure('Name', 'Frequency-Voltage-Domains Graphs');
			
			% TODO: NOT FULLY TESTED!
			rp = 2;
			cp = 2;
			taken = 1;
			stop = 0;
			for i=1:obj.vd
				rp = 1+i;
				if (mod(rp,4)==0) || (mod(rp,3)==0) || (rp == 2)
					taken = min([ceil(rp/3) ceil(rp/2) ceil(rp/4)]);
					for j=2:2:rp			
						if (rp-taken)*j >= obj.vd
							stop = 1;
							cp = j;
							break;
						end			
					end		
				end
				if (stop)
					break;
				end
			end
			
			posif = [];
			posiv = [];
			for i=1:taken
				ind = (i-1)*cp+1;
				posif = [posif, ind:(ind-1)+cp/2];
				posiv = [posiv ind+cp/2:(ind-1)+cp];
			end
			
			ax1 = subplot(rp,cp,posif);
			plot(t, f);
			grid on,xlabel('Time [s]'), ylabel('Applied Freq [GHz]');
			xlim([0 obj.tsim]);
			
			ax2 = subplot(rp,cp,posiv);
			plot(t, v);
			grid on,xlabel('Time [s]'), ylabel('Applied Voltage [V]');
			xlim([0 obj.tsim]);
			
			axi = [];
			dim = sum(obj.VDom);
			for d=1:obj.vd
				axi = [axi subplot(rp,cp,taken*cp+d)];
				hold on;
				idx = obj.VDom(:,d).*[1:1:obj.Nc]';
				idx = idx(idx>0);
				plot(t, f(:,idx));
				idx = sum((v(:,d) * ones(1,obj.FV_levels)) > obj.FV_table(:,1)',2) + 1;
				plot(t, obj.FV_table(idx,3), '--', 'LineWidth',1, 'Color', "#0000FF");
				grid on,xlabel('Time [s]'), ylabel(strcat("Applied Freq [GHz] - Domain: ", int2str(d)));
				xlim([0 obj.tsim]);
			end

			linkaxes([ax1,ax2, axi],'x');
		end
		function [fig] = perfplot(obj, f, w)		

			font_size = 6;
			
			fig = figure('Name', 'Performance Comparison');
			movegui(fig, 'northeast');
			
			ax1 = subplot(3,4,[1:4]);
			smref = round((size(f,1)-1)/(size(obj.frplot,1)-1));
			smf = round((size(obj.frplot,1)-1)/(size(f,1)-1));
			
			gr = repelem(f(2:1:end,:),max(smf,1),1) - repelem(obj.frplot(2:end,:),max(smref,1),1);
			plot(obj.tsim*[0:(1/size(gr,1)):1]', [zeros(1,obj.Ni_c); gr]);
			grid on,xlabel('Time [s]'), ylabel('Reference Difference [GHz]');
			
			subplot(3,4,[5:6]);
			av = sum(obj.frplot(2:end,:)) / size(obj.frplot(2:end,:),1);
			b = bar(sum(gr)/size(gr,1)./av*100);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom', 'FontSize', font_size);
			end
			grid on,xlabel('Cores'), ylabel('Mean Difference [%]');
			
			subplot(3,4,7);
			%av = sum(obj.frplot(2:end,:)) / size(obj.frplot(2:end,:),1);
			b = bar(sum(gr)*obj.VDom/size(gr,1)./(av*obj.VDom)*100);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			grid on,xlabel('Domains'), ylabel('Mean Difference [%]');
			
			subplot(3,4,8);
			%av = sum(obj.frplot(2:end,:)) / size(obj.frplot(2:end,:),1);
			b = bar(sum(sum(gr))/size(gr,1)/sum(av)*100);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			grid on,xlabel('Total'), ylabel('Total Mean Difference [%]');
			
			% % %
			gr = w / (size(obj.wrplot,3)-1) * 100;
			subplot(3,4,[9:10]);
			b = bar(gr);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom', 'FontSize',font_size);
			end
			grid on,xlabel('Cores'), ylabel('Application Completion [%]');
			
			subplot(3,4,11);
			grv = diag(gr)*obj.VDom;
			tt = grv>0;
			min_grv = [];
			for v=1:obj.vd
				min_grv = [min_grv min(grv(tt(:,v),v))];
			end
			b = bar([obj.VDom'*gr ./ sum(obj.VDom)' min_grv']);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom', 'FontSize',font_size);
			end
			if obj.vd == 1
				b.FaceColor = 'flat';
				b.CData(2,:) = [0.8500 0.3250 0.0980];
				xticklabels({'Mean', 'Min'});
			end
			grid on,xlabel('Domains'), ylabel('Application Completion (Mean - Min) [%]');
			
			subplot(3,4,12);
			b = bar([sum(gr)/obj.Nc min(min_grv)]);
			for pb = 1:length(b)
				xtips = b(pb).XEndPoints;
				ytips = b(pb).YEndPoints;
				labels = string(round(b(pb).YData,2))+'%';
				text(xtips,ytips,labels,'HorizontalAlignment','center',...
					'VerticalAlignment','bottom');
			end
			b.FaceColor = 'flat';
			b.CData(2,:) = [0.8500 0.3250 0.0980];
			grid on,xlabel(''), ylabel('Application Completion [%]'), xticklabels({'Mean', 'Min'});
			
			%linkaxes([ax1,ax2,ax3],'x');
			
		end
		%TODO remove change these 2
		function [fig] = perfuplot(obj, u, t)
			
			sim_mul = round(obj.Ts_input/obj.Ts);					
			
			% Perf
			fig = figure('Name', 'Performance Comparison')
			subplot(3,4,[1:4]);
			umpc = u-obj.zrplot;
			gr = obj.urplot(2:end,:) - umpc(sim_mul+1:sim_mul:end,:);
			plot(t, [zeros(1,obj.Ni_c); gr]);
			grid on,xlabel('Time [s]'), ylabel('Reference Difference [W]');
			subplot(3,4,[5:8]);
			plot(t, [0; sum(gr,2)]);
			grid on,xlabel('Time [s]'), ylabel('Total Difference [W]');
			subplot(3,4,[9,10]);
			av = sum(obj.urplot(2:end,:)) / size(gr,1);
			bar(sum(gr)/size(gr,1)./av*100);
			grid on,xlabel('Cores'), ylabel('Mean Difference [%]');
			subplot(3,4,11);
			bar((sum(gr)*obj.VDom)/size(gr,1)./(av*obj.VDom)*100);
			grid on,xlabel('Domains'), ylabel('Mean Difference [%]');
			subplot(3,4,12);
			av = sum(obj.urplot(2:end,:)) / size(gr,1);
			bar(sum(sum(gr))/size(gr,1)/sum(av)*100);
			grid on,xlabel(''), ylabel('Total Mean Difference [%]');			
		end
		function [fig] = perfuplotcomp(obj, u1,u2, t)
			
			sim_mul = round(obj.Ts_input/obj.Ts);					
			
			% Perf
			fig = figure('Name', 'Performance Comparison')
			subplot(3,4,[1:2]);
			umpc = u1-obj.zrplot;
			gr1 = obj.urplot(2:end,:) - umpc(sim_mul+1:sim_mul:end,:);
			plot(t, [zeros(1,obj.Ni_c); gr1]);
			grid on,xlabel('Time [s]'), ylabel('Reference Difference [W]');
			subplot(3,4,[3:4]);
			plot(t, [0; sum(gr1,2)]);
			grid on,xlabel('Time [s]'), ylabel('Total Difference [W]');
			
			subplot(3,4,[5:6]);
			umpc = u2-obj.zrplot;
			gr2 = obj.urplot(2:end,:) - umpc(sim_mul+1:sim_mul:end,:);
			plot(t, [zeros(1,obj.Ni_c); gr2]);
			grid on,xlabel('Time [s]'), ylabel('Reference Difference [W]');
			subplot(3,4,[7:8]);
			plot(t, [0; sum(gr2,2)]);
			grid on,xlabel('Time [s]'), ylabel('Total Difference [W]');
			
			subplot(3,4,[9,10]);
			av1 = sum(obj.urplot(2:end,:)) / size(gr1,1);
			av2 = sum(obj.urplot(2:end,:)) / size(gr2,1);
			BAR = [(sum(gr1)/size(gr1,1)./av1*100)' ...
				(sum(gr2)/size(gr2,1)./av2*100)'];
			bar(1:length(sum(gr1)),BAR);
			grid on,xlabel('Cores'), ylabel('Mean Difference [%]');
			subplot(3,4,11);
			BAR = [((sum(gr1)*obj.VDom)/size(gr1,1)./(av1*obj.VDom)*100)' ...
				((sum(gr2)*obj.VDom)/size(gr2,1)./(av2*obj.VDom)*100)' ];
			bar(1:length((sum(gr2)*obj.VDom)), BAR);
			grid on,xlabel('Domains'), ylabel('Mean Difference [%]');
			subplot(3,4,12);
			av1 = sum(obj.urplot(2:end,:)) / size(gr1,1);
			av2 = sum(obj.urplot(2:end,:)) / size(gr2,1);
			BAR = [sum(sum(gr1))/size(gr1,1)/sum(av1)*100 ...
				sum(sum(gr2))/size(gr2,1)/sum(av2)*100];
			bar(1, BAR);
			grid on,xlabel(''), ylabel('Total Mean Difference [%]');			
		end
		function [oc] = stats_analysis(obj, x,u,f,v,w)
			% Since workload is changing (even if it is averaged), takings values for
			% each core will give us high variability. So I'm taking average chip
			% values:

			% exceeded P Budget (W and %)
			xcp = x(2:end,:)*obj.C(1:obj.Nc,:)';
			ucp = u(2:end,:);
			fcp = f(2:end,:);
			vcp = v(2:end,:);
			
			%%%%%%%%%%%%%%%%%%%%%%%%%
			% Power

			pwct = max(round(size(ucp,1)/max(length(obj.tot_pw_budget)-1,1)),1);
			pwtb = repelem(obj.tot_pw_budget(min(2, length(obj.tot_pw_budget)):end), pwct ,1);

			% here it is difficult, because if I make it ceil, it will always be ok for
			% pw budget going down, but never for power budget going up. If I take
			% floor it will be the opposite. Hope that round will be ok!
			delay = round((obj.delay_F_max+obj.delay_V_max) / obj.Ts);
			sidx = delay + pwct;
			%
			pwgr = sum(ucp(1+sidx:end,:),2) - pwtb(1:end-sidx);
			pwgrP = pwgr ./ pwtb(1:end-sidx);
			
			pwgr = pwgr(pwgr > 0);
			pwgrP = pwgrP(pwgrP > 0);
			
			if isempty(pwgr)
				pwgr = 0;
			end
			if isempty(pwgrP)
				pwgrP = 0;
			end

			power.exMaxW = max(pwgr);
			power.exMaxP = max(pwgrP);
			power.ex95W = prctile(pwgr, 95);
			power.ex95P = prctile(pwgrP, 95);
			power.exAvW = mean( pwgr );
			power.exAvP = mean( pwgrP );
			power.exSDW = std(pwgr);
			power.exSkew = skewness(pwgr);
			power.exKurt = kurtosis(pwgr);
			%
			power.exTime = sum( sum(ucp(1+sidx:end,:),2) > pwtb(1:end-sidx) ) *obj.Ts; 

			%%%%%%%%%%%%%%%%%%%%%%%%%
			% Temperature

			temp.Max = max(xcp(:)) -273.15;
			temp.p95 = prctile(xcp(:),95) - 273.15;
			temp.Av = mean(xcp(:)) -273.15;
			temp.AvSD = std(mean(xcp) - 273.15);
			temp.AvSkew = skewness(mean(xcp) - 273.15);
			temp.AvKurt = kurtosis(mean(xcp) - 273.15);
			temp.SkewMean = mean(skewness(xcp));
			temp.KurtMean = mean(kurtosis(xcp));
			%
			exTC = (xcp-obj.core_crit_temp) .* (xcp > obj.core_crit_temp);
			exTC = exTC(exTC > 0);
			if isempty(exTC)
				exTC = 0;
			end
			temp.exMaxCr = max(max( exTC ));
			temp.ex95pCr = mean(prctile( exTC ,95));
			temp.ex80pCr = mean(prctile( exTC ,80));
			temp.exAvCr = mean(mean( exTC ));
			temp.exSDCr = std(mean( exTC ));
			temp.exSkewCr = skewness(mean( exTC ));
			temp.exKurtCr = kurtosis(mean( exTC ));
			%
			exTn = (xcp-obj.core_crit_temp) .* (xcp > obj.core_crit_temp-obj.T_margin);
			exTn = exTn(exTn > 0);
			if isempty(exTn)
				exTn = 0;
			end
			temp.exAvMn = mean(mean( exTn ));
			temp.exSDMn = std(mean( exTn ));
			temp.exSkewMn = skewness(mean( exTn ));
			temp.exKurtMn = kurtosis(mean( exTn ));
			%
			temp.exTotTimeCr = sum(sum( xcp > obj.core_crit_temp ) *obj.Ts);
			temp.exSDTimeCr = std(sum( xcp > obj.core_crit_temp ) *obj.Ts);
			temp.exSkewTimeCr = skewness(sum( xcp > obj.core_crit_temp ) *obj.Ts);
			temp.exKurtTimeCr = kurtosis(sum( xcp > obj.core_crit_temp ) *obj.Ts);
			%
			temp.exTotTimeMn = sum(sum( xcp > obj.core_crit_temp -obj.T_margin) *obj.Ts);
			temp.exSDTimeMn = std(sum( xcp > obj.core_crit_temp -obj.T_margin) *obj.Ts);
			temp.exSkewTimeMn = skewness(sum( xcp > obj.core_crit_temp -obj.T_margin) *obj.Ts);
			temp.exKurtTimeMn = kurtosis(sum( xcp > obj.core_crit_temp -obj.T_margin) *obj.Ts);
			% no need these two: I can divide per obj.Nc
			temp.exAvTimeCr = mean(sum( xcp > obj.core_crit_temp ) *obj.Ts);
			temp.exAvTimeMn = mean(sum( xcp > obj.core_crit_temp - obj.T_margin) *obj.Ts);

			%%%%%%%%%%%%%%%%%%%%%%%%%
			% freq

			freq.Max = max(fcp(:));
			freq.min = min(fcp(:));
			freq.Av = mean(fcp(:));
			freq.AvSD = std(mean(fcp));
			freq.AvSkew = skewness(mean(fcp));
			freq.AvKurt = kurtosis(mean(fcp));
			freq.SkewMean = mean(skewness(fcp));
			freq.KurtMean = mean(kurtosis(fcp));

			% oscillations:
			%want to study for 1ms, 10ms, 100ms, 1s
			mdim = [1e-3 1e-2 1e-1 1] / obj.Ts_ctrl;
			tt = mdim > obj.tsim / obj.Ts_ctrl;
			mdim(tt) = obj.tsim / obj.Ts_ctrl;
			
			for md = 1:length(mdim)
				dd = mdim(md);
				ve = zeros(size(fcp) - [dd 0]);
				for i=1:length(fcp)-dd+1
					idx = i+dd-1;
					ve(i,:) = std(fcp(i:idx,:));
				end

				freq.OscMax(md) = max(ve(:));
				freq.OscMaxMean(md) = mean(max(ve));
				freq.OscMaxSD(md) = std(max(ve));
				freq.OscMaxSkew(md) = skewness(max(ve));
				freq.OscMaxKurt(md) = kurtosis(max(ve));
				freq.OscAv(md) = mean(ve(:));
				freq.OscAvSD(md) = std(mean(ve));
				freq.OscSDMean(md) = mean(std(ve));
				freq.OscAvSkew(md) = skewness(mean(ve));
				freq.OscAvKurt(md) = kurtosis(mean(ve));
				freq.OscSkewMean(md) = mean(skewness(ve));
				freq.OscKurtMean(md) = mean(kurtosis(ve));	
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%
			% Voltage

			vdd.Max = max(vcp(:));
			vdd.min = min(vcp(:));
			vdd.Av = mean(vcp(:));
			vdd.AvSD = std(mean(vcp));
			vdd.AvSkew = skewness(mean(vcp));
			vdd.AvKurt = kurtosis(mean(vcp));
			vdd.SkewMean = mean(skewness(vcp));
			vdd.KurtMean = mean(kurtosis(vcp));

			% oscillations:
			%want to study for 1ms, 10ms, 100ms, 1s
			mdim = [1e-3 1e-2 1e-1 1] / obj.Ts_ctrl;
			tt = mdim > obj.tsim / obj.Ts_ctrl;
			mdim(tt) = obj.tsim / obj.Ts_ctrl;
			
			for md = 1:length(mdim)
				dd = mdim(md);
				ve = zeros(size(vcp) - [dd 0]);
				for i=1:length(vcp)-dd+1
					idx = i+dd-1;
					ve(i,:) = std(vcp(i:idx,:));
				end

				vdd.OscMax(md) = max(ve(:));
				vdd.OscMaxMean(md) = mean(max(ve));
				vdd.OscMaxSD(md) = std(max(ve));
				vdd.OscMaxSkew(md) = skewness(max(ve));
				vdd.OscMaxKurt(md) = kurtosis(max(ve));
				vdd.OscAv(md) = mean(ve(:));
				vdd.OscAvSD(md) = std(mean(ve));
				vdd.OscSDMean(md) = mean(std(ve));
				vdd.OscAvSkew(md) = skewness(mean(ve));
				vdd.OscAvKurt(md) = kurtosis(mean(ve));
				vdd.OscSkewMean(md) = mean(skewness(ve));
				vdd.OscKurtMean(md) = mean(kurtosis(ve));	
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%
			% Perf
			smref = round((size(fcp,1))/(size(obj.frplot(2:end,:),1)));
			smf = round((size(obj.frplot(2:end,:),1))/(size(fcp,1)));

			fd = repelem(obj.frplot(2:end,:),max(smref,1),1) - repelem(fcp,max(smf,1),1);

			n2 = reshape(fd,[],1);
			perf.fd2norm = norm(n2) / sqrt(obj.Nc) / sqrt(ceil(obj.tsim / obj.Ts_ctrl));
			%perf.fdAv2norm = mean(norm
			perf.fdMax = max(fd(:));
			perf.fdmin = min(fd(:));
			perf.fdAv = mean(fd(:));
			perf.fdSum = sum(fd(:));
			perf.fdAvSD = std(mean(fd));
			perf.fdAvSkew = skewness(mean(fd));
			perf.fdAvKurt = kurtosis(mean(fd));

			perf.fdSkewMean = mean(skewness(fd));
			perf.fdKurtMean = mean(kurtosis(fd));
			perf.fdSDMean = mean(std(fd));


			perf.wlMax = max(w);
			perf.wlAv = mean(w);
			perf.wlmin = min(w);
			perf.wlSD = std(w);
			perf.wlSkew = skewness(w);
			perf.wlKurt = kurtosis(w);
			
			%%%%%%%%%%%%%%%%%%%%%%%%%
			% out
			
			oc.power = power;
			oc.temp = temp;
			oc.freq = freq;
			oc.vdd = vdd;
			oc.perf = perf;			
		end
		function [fig] = paperTPplot(obj,t,x,u)
			add_state = obj.Ns-obj.Nc*2; %TODO
			fig = figure('Name', 'Metrics'); %, 'Position', get(0, 'Screensize'));
			movegui(fig, 'northwest');
			ax1=subplot(4,1,1);
			plot(t, x(:,1+end-add_state:end)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Temp Rack [°C]'), hold on;
			legend("Heat-sink", "PCB", "Motherboard", "Air");
			%yline(obj.core_crit_temp-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			ax2=subplot(4,1,2);
			plot(t, x(:,2:2:end-add_state)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Copper [°C]'),hold on;
			if max(max(x(:,2:2:end-add_state))) >= (obj.core_crit_temp-(0.1*(obj.core_crit_temp-273)))
				yline(obj.core_crit_temp-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			end
			ax3=subplot(4,1,3);
			plot(t, x(:,1:2:end-add_state)-273.15);
			grid on,xlabel('Time [s]'),ylabel('Silicon [°C]'), hold on;
			if max(max(x(:,1:2:end-1))) >= (obj.core_crit_temp-(0.1*(obj.core_crit_temp-273)))
				yline(obj.core_crit_temp-273.15, '--', 'LineWidth',1,   'Color', '#CC0000');
			end
			
			ax4=subplot(4,1,4);
			gr = sum(u,2);
			plot(t,gr, '.');
			grid on,title('Total Power'),ylabel('[W]'),xlabel('Time [s]'),hold on;
			%if max(gr) >= obj.usum(1)*0.2
			plot([0:max(length(obj.tot_pw_budget)-1,1)]*(obj.tsim/max(length(obj.tot_pw_budget)-1,1)),[obj.tot_pw_budget], '--', 'LineWidth',1,  'Color', '#CC0000');
			
			linkaxes([ax1,ax2,ax3,ax4],'x');
		end
		function [] = savetofile(obj, fig, path_name, bmpres)
			
			% Mange screens
			graph_dpi = 416;
			file_ext = ".png";
			%bmpres = 1; %1.5;
						
			%save as fig
			saveas(fig, strcat(path_name,".fig"));
			%change res for bmp:
			fig.Visible = 'off';
			fig.Position = [1, 1, 1920*bmpres,1080*bmpres];
			exportgraphics(fig, strcat(path_name,file_ext), 'Resolution', graph_dpi);
			close(fig);			
		end
		function [] = saveall(obj, xop, uop, fop, vop, wop, t, path, name)
			
			if nargin < 8 || (strlength(path)<1)
				if isunix
					path = "/tmp/MATLAB-Figures";
				else
					path = "C:\temp\MATLAB-Figures";
				end
			end
			if nargin < 9 || (strlength(name)<1)
				name = string(datetime('now','Format','dd-MMM-yyyy HH_mm_ss_SSS'));
			end

			if ~exist(path, 'dir')
				mkdir(path);
			end
			if isunix
				separator_os = "/";
			else
				separator_os = "\";
			end
			
			fig = [];
			namefig = "";
			res = 1;

			for i=1:5
				switch i
					case 1
						fig = obj.xutplot(xop, uop);
						namefig = "TP";
						res = 1.5;
					case 2
						fig = obj.powerconstrplot(uop);
						namefig = "Pw";
						res = 1;
					case 3
						fig = obj.tempconstrplot(xop);
						namefig = "T";
						res = 1;
					case 4
						fig = obj.perfplot(fop, wop * (size(obj.wrplot,3)-1) / 100 );
						namefig = "Perf";
						res = 1;
					case 5
						fig = obj.fvplot(t, fop,vop);
						namefig = "FV";
						res = 1.5;
				end
				if ~isempty(fig)
					path_name = strcat(path,separator_os,namefig, "-", name);
					obj.savetofile(fig, path_name, res);
				end
			end
		end
		
	end %methods
	
	% ==============
	% Sim
	% ==============
	methods
		obj = generate_reference(obj, show, alpha_u, alpha_f, alpha_w)			
		[] = simulate_aut(obj, lin, x_init, ts)			
		[wlop] = max_exec_test(obj)
		
		% Simulations
		[xplot, uplot] = launch_lin_mpc_sim(obj, robust, obs, comp)
		[cpxplot, cpuplot, cpfplot, cpvplot, wlop] = launch_cp_sim(obj, robust, mode, show)
		[cpxplot,cpuplot, cpfplot, cpvplot, xlplot, wlop] = launch_cpmpc_sim(obj, obs)
		[cpxplot, cpuplot, cpfplot, cpvplot, wlop] = launch_ncp_sim(obj, robust, show)
		[cpxplot, cpuplot, cpfplot, cpvplot, wlop] = launch_occ_sim(obj, Tall, Ttick, robust, show)
		[cpxplot, cpuplot, cpfplot, cpvplot, wlop] = launch_arm_sim(obj, robust, mode, show)
		[cpxplot,cpuplot, cpfplot, cpvplot, xlplot, wlop] = launch_nmpc_sim(obj, obs)
	end %methods

	%% Setup
	methods
		function obj = create_thermal_model(obj, th_model_ver)
			
			if nargin < 2
				th_model_ver = obj.thermal_model_ver;
			end
			
			% ==============
			% C & D
			% ==============
			obj.C = zeros(obj.Nout,obj.Ns);
			k=-2;
			for i=1:obj.Nc
				k=k+2;
				obj.C(i,k+1)=1;
			end
			obj.C(end, end) = 1;
			obj.D = zeros(obj.Nout, obj.Ni);
			
			[obj.Ac_nom, obj.Bc_nom] = obj.lin_model_create(0, th_model_ver, 0);
			[obj.Ac_true, obj.Bc_true] = obj.lin_model_create(1, th_model_ver, 0);
			
			%obj.x_init = ones(obj.Ns,1)*obj.temp_amb;
			% TODO: remove this after leakage is fixed:
			leak_store = obj.exp_leakage;
			obj.exp_leakage = 0;
			%obj.min_pw_red = obj.power_compute(ones(obj.Nc,1)*obj.F_min,ones(obj.Nc,1)*obj.V_min,(obj.core_crit_temp)*ones(obj.Nc,1),[zeros(obj.Nc,obj.ipl-1) ones(obj.Nc,1)],ones(obj.Nc,1));
			obj.min_pw_red = obj.power_compute(ones(obj.Nc,1)*obj.F_min,ones(obj.Nc,1)*obj.V_Max,(obj.core_crit_temp)*ones(obj.Nc,1),[zeros(obj.Nc,obj.ipl-1) ones(obj.Nc,1)],ones(obj.Nc,1));
			obj.exp_leakage = leak_store;
			
			if rank([obj.C(1:obj.Nc,:);obj.C(1:obj.Nc,:)*obj.Ac_nom;obj.C(1:obj.Nc,:)*obj.Ac_nom*obj.Ac_nom]) == obj.Ns
			%This is not working
			%(rank(obsv(obj.Ac_nom, obj.C)) == obj.Ns)
				disp("The system is Observable!")
				obj.observable = 1;
			else
				disp("MAY be NOT Observable!")
				obj.observable = 0;
			end
		end
		function obj = create_core_pw_noise(obj)
			Covr = chol(obj.pw_gvar);
			obj.core_pw_noise_char = ones(obj.Nc, 1);
			for c=1:obj.Nc
				stop = 0;
				while(~stop)	
					z = repmat(obj.pw_gmean,1,1) + randn(1,1)*Covr * obj.pw_3sigma_on3/3;	
					if (z>obj.pw_glim(1)) && (z<obj.pw_glim(2))
						stop = 1;
					end
				end %while
				obj.core_pw_noise_char(c) = obj.core_pw_noise_char(c) + z;
			end %for
		end %func
		function obj = create_thermal_model_noise(obj)
			ddiv = 30;
			dmean = 1.01;
			obj.ThMoNoise = (randn(obj.Nc, 8) / ddiv + dmean);
		end
	end
	
	%% Constructor & Utilities
	methods
		function obj = hpc_system(Nc,Nh,Nv,vd)
			
			% Pre Initialization %%
			% Any code not using output argument (obj)
			if (nargin < 3)
				obj.VDom = zeros(obj.Nc, obj.vd);

				obj.VDom([1:3 7:9 13:15], 1) = 1;
				obj.VDom([4:6 10:12 16:18], 2) = 1;
				obj.VDom([19:21 25:27 31:33], 3) = 1;
				obj.VDom([22:24 28:30 34:36], 4) = 1;
			end
			switch nargin
				case 3
					obj.Nc = Nc;
					obj.Nh = Nh;
					obj.Nv = Nv;
                    obj.VDom = eye(obj.Nc);
                    obj.vd = obj.Nc;
				case 4
					obj.Nc = Nc;
					obj.Nh = Nh;
					obj.Nv = Nv;
					obj.vd = vd;
				%else
					%error
			end

			% Object Initialization %%
			% Call superclass constructor before accessing object
			% You cannot conditionalize this statement
			%obj = obj@BaseClass1(args{:});

			% Post Initialization %%
			% Any code, including access to object
			%obj = obj.create_thermal_model();
			%obj = obj.init_mpc();
			%obj = obj.generate_reference();
			%obj.customColormap = 
			load('-mat', "CustomColorMap3.mat");
			obj.customColormap = CustomColormap3;
			%TODO:
			obj = obj.create_core_pw_noise();
			obj = obj.create_thermal_model_noise();
			
			obj.tot_pw_budget = ((obj.core_Max_power-0.5)*obj.Ni_c-obj.Ni_c)*ones(2,1);
			obj.quad_pw_budget = ((obj.core_Max_power-0.5)*obj.Ni_c-obj.Ni_c)/obj.vd*ones(2,obj.vd);
			
		end
		function totSize = get_size(obj) 
		   props = properties(obj); 
		   totSize = 0; 

		   for ii=1:length(props) 
			  currentProperty = getfield(obj, char(props(ii))); 
			  s = whos('currentProperty'); 
			  totSize = totSize + s.bytes; 
		   end

		   fprintf(1, '%d bytes\n', totSize); 
		end
	end
	
	%% Dependent Variables
	methods
		% State Dimension
		function value = get.Ns(obj)
			value = obj.Nc*2+4;
		end
		% Input Dimension
		function value = get.Ni(obj)
			value = obj.Nc+1;
		end
		function value = get.Ni_nl(obj)
			value = obj.Nc+obj.vd;
		end
		% Controllable Inputs
		function value = get.Ni_c(obj)
			value = obj.Nc;
		end
		% Uncontrollable Inputs
		function value = get.Ni_nc(obj)
			value = obj.Ni-obj.Ni_c;
		end
		% Outputs
		function value = get.Nout(obj)
			value = obj.Nc+1;
		end	
		
		function value = get.static_pw_coeff(obj)
			value = [obj.leak_vdd, obj.leak_temp, obj.leak_process];
		end
		function value = get.FV_levels(obj)
			value = size(obj.FV_table,1);
		end
		function value = get.V_Max(obj)
			value = obj.FV_table(end,1);
		end
		function value = get.V_min(obj)
			value = obj.FV_table(1,1);
		end
		function value = get.F_min(obj)
			value = obj.FV_table(1,2);
		end
		function value = get.F_Max(obj)
			value = obj.FV_table(end,3);
		end
		function value = get.core_Max_power(obj)
			value = 10.5348;
			%TODO: here put the call to the function model
		end
		function value = get.core_min_power(obj)
			value = 0.8939;
			%TODO: here put the call to the function model
		end
		function value = get.ipl(obj)
			value = size(obj.ceff_pw_coeff,2);
		end
		function value = get.polyFV(obj)
			value = polyfit(obj.FV_table(:,3),obj.FV_table(:,1), 3);
		end
		function value = get.quantum_instr(obj)
			value = obj.F_Max * 1e9 * (obj.quantum_us * 1e-6);
			%value = 3 * 1e9 * (obj.quantum_us * 1e-6);
		end
		
		% Matrices
		function value = get.Ad_nom(obj)
			disc = c2d(ss(obj.Ac_nom, obj.Bc_nom, obj.C, obj.D), obj.Ts);
			value = disc.A;
		end
		function value = get.Bd_nom(obj)
			disc = c2d(ss(obj.Ac_nom, obj.Bc_nom, obj.C, obj.D), obj.Ts);
			value = disc.B;
		end
		function value = get.Ad_true(obj)
			disc = c2d(ss(obj.Ac_true, obj.Bc_true, obj.C, obj.D), obj.Ts);
			value = disc.A;
		end
		function value = get.Bd_true(obj)
			disc = c2d(ss(obj.Ac_true, obj.Bc_true, obj.C, obj.D), obj.Ts);
			value = disc.B;
		end
		function value = get.Ad_mpc(obj)
			disc = c2d(ss(obj.Ac_nom, obj.Bc_nom, obj.C, obj.D), obj.Ts_mpc);
			value = disc.A;
		end
		function value = get.Bd_mpc(obj)
			disc = c2d(ss(obj.Ac_nom, obj.Bc_nom, obj.C, obj.D), obj.Ts_mpc);
			value = disc.B;
		end
		function value = get.Ad_ctrl(obj)
			disc = c2d(ss(obj.Ac_nom, obj.Bc_nom, obj.C, obj.D), obj.Ts_ctrl);
			value = disc.A;
		end
		function value = get.Bd_ctrl(obj)
			disc = c2d(ss(obj.Ac_nom, obj.Bc_nom, obj.C, obj.D), obj.Ts_ctrl);
			value = disc.B;
		end
		function value = get.Ad_obs(obj)
			disc = c2d(ss(obj.Ac_nom, obj.Bc_nom, obj.C, obj.D), obj.Ts_obs);
			value = disc.A;
		end
		function value = get.Bd_obs(obj)
			disc = c2d(ss(obj.Ac_nom, obj.Bc_nom, obj.C, obj.D), obj.Ts_obs);
			value = disc.B;
		end
		
		
		%% SET
		function obj = set.model_variation(obj, val)
			%check modifications
			if val ~= obj.model_variation 
				if val
					obj = obj.create_core_pw_noise();
					obj = obj.create_thermal_model_noise();
				else
					%reset: 
					obj.core_pw_noise_char = ones(obj.Nc,1);
					obj.ThMoNoise = ones(obj.Nc, 8);					
				end
			end
			obj.model_variation = val;
		end
		
		
	end

end
