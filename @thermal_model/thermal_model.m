classdef thermal_model
	%THERMAL_MODEL Summary of this class goes here
	%   Detailed explanation goes here
	
	properties		
		%% System Structure
		Nc (1,1) {mustBePositive, mustBeInteger} ...
			= 36;		% Number of Cores
		Nh (1,1) {mustBePositive, mustBeInteger} ...
			= 6;		% Number of rows
		Nv (1,1) {mustBePositive, mustBeInteger} ...
			= 6;		% Number of cols
		vd (1,1) {mustBePositive, mustBeInteger} ...
			= 4;		% Voltage/Power Domains
		VDom {mustBeNonnegative, mustBeNumericOrLogical, mustBeNonempty, mustBeLessThanOrEqual(VDom,1)} ...
			= 1;		% Structure of the Voltage Domains Nc x vd
		
		%% Matlab/Simulation
		model_ver (1,1) {mustBeNonnegative, mustBeInteger} ...
			= 0;				% Version of the Thermal Model
		model_deviations (1,1) {mustBeNonnegative, mustBeNumericOrLogical} ...
			= 1;				% If the true system should be != from nominal
		%measure_noise = 1;
		sensor_noise (1,1) {mustBeNonnegative, mustBeNumericOrLogical} ...
			= 1;
		%T_noise_max = 1.0;
		sensor_noise_amplitude (3,1) {mustBeNonnegative, mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1.0; 1.0; 1.0];	% Process, Voltage, Temperature amplitudes.
		Ts (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 250e-6;			% Discretization Time
		
		Ac_nom {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1];
		Bc_nom {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1];
		Ac_true {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1];
		Bc_true {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1];
		
		C {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1];
		D {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1];

		temp_amb (1,1) {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 25.0 + 273.15;	% External Ambient Temperature

		param_dev_per {mustBeNonnegative, mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 1;

		% fp = FloorPlan
		RC_fp_dim (:,:,4) {mustBeNonnegative, mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= ones(1,1,4);
		%		(:,:,1) = North
		%		(:,:,2) = East
		%		(:,:,3) = South
		%		(:,:,4) = West
		CPw_fp_dim (:,:,4) {mustBeNonnegative, mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= ones(1,1,4);
		%		(:,:,1) = North
		%		(:,:,2) = East
		%		(:,:,3) = South
		%		(:,:,4) = West
		R_fp_material {mustBeNonnegative, mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 1;
		%		(:,:,1) = Si
		%		(:,:,2) = Cu
		C_fp_material {mustBeNonnegative, mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= 1;
		%		(:,:,1) = Si
		%		(:,:,2) = Cu

		% Thickness
		t_si (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 5e-4; %3e-4;
		t_cu (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 1e-3; %1.1e-3;
		t_comp {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= [1.5e-2, 1e-3, 2e-3, 10e-2];
		%al, pcb, mb, air (%reverse order of pos)

		% Width
		wid_comp {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= [0.0416, 0.0366, 10e-2, 20e-2];
		%al, pcb, mb, air (%reverse order of pos)

		% Length
		len_comp {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= [0.0416, 0.0366, 10e-2, 20e-2];
		%al, pcb, mb, air (%reverse order of pos)


		% Additional vertical resistance
		R_TIM1 (1,1) {mustBeNonnegative, mustBeNumeric, mustBeFinite} ...
			= 1;%0.25; %6 4.5
		R_TIM2 (1,1) {mustBeNonnegative, mustBeNumeric, mustBeFinite} ...
			= 2.5; %6; %4; 0.75;
		case_fan_dis (1,1) {mustBeNonnegative, mustBeNumeric, mustBeFinite} ...
			= 0.25; %0.85 %0.1 - 1
		case_fan_nom_speed (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 1000; %fan_nom_speed
		al_fan_dis (1,1) {mustBeNonnegative, mustBeNumeric, mustBeFinite} ...
			= 10; %heatsink_fan_dis
		al_fins_coeff (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 3; %fins_coeff

		air_factor (1,1) {mustBePositive, mustBeNumeric, mustBeFinite} ...
			= 100;
		%air_time_factor = 10; %20

		pw2therm_coeff = 1.05; %0.90 %0.41; 0.7;

		add_sensor_outputs (1,1) {mustBePositive, mustBeInteger} ...
			= 0;

		sensor_active {mustBeNonnegative, mustBeNumericOrLogical, mustBeVector} ...
			= [0 0 0 0];

end

	properties(Dependent)
		%% Matlab/Simulation
		Ns;							% Number of States
		Ni;							% Number of inputs
		Nout;						% Number of "Outputs"

		% Matlab Matrices
		Ad_nom;
		Bd_nom;
		Ad_true;
		Bd_true;
	end

	properties(SetAccess=protected, GetAccess=public)
		observable = 0;
	end
	
	properties(SetAccess=immutable, GetAccess=public)

		%%% Value that, if changed, I need to change the code

		add_states		= 4;
		add_inputs		= 1;
		param_dev_dim2	= 8;

		full_model_layers = 2;

		extt_rows = 1;
		extb_rows = 1;
		extl_cols = 1;
		extr_cols = 1;

		air_pos		= 0;
		mb_pos		= 1;
		pcb_pos		= 2;
		al_pos		= 3;	

		north_pos	= 1;
		east_pos	= 2;
		south_pos	= 3;
		west_pos	= 4;

		si_pos		= 1;
		cu_pos		= 2;

		si_pcb_fact = 0.1;
		pcb_mb_fact = 15;
		mb_air_fact = 10;

		%%% PHYSICAL values

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

		%{
		A_dim
		B_dim
		C_dim
		D_dim
		%}
	end
	
	%% Gerenal Methods
	methods		
		function obj = thermal_model(Nc,Nh,Nv)
			
			% Pre Initialization %%
			% Any code not using output argument (obj)
			if (nargin < 3)
				warning("[TM Error]: missing one or more values. Default Parameters will be used")
			else
				obj.Nc = Nc;
				obj.Nh = Nh;
				obj.Nv = Nv;
			end

			% Object Initialization %%
			% Call superclass constructor before accessing object
			% You cannot conditionalize this statement
			%obj = obj@BaseClass1(args{:});
			obj.sensor_active = zeros(add_states,1);
			obj = default_floorplan_config();

			% Post Initialization %%
			% Any code, including access to object
			obj = create_thermal_model_noise();
			obj = create_thermal_model();


			get_size(obj);
		end
	end
	methods(Static)
		function totSize = get_size(class) 
			props = properties(class); 
			totSize = 0;
			
			for ii=1:length(props) 
  			currentProperty = getfield(class, char(props(ii))); 
  			s = whos('currentProperty'); 
  			totSize = totSize + s.bytes; 
			end
			
			disp(strcat("Size: ", num2str(totSize), " bytes")); 
		end
	end

	%% Models
	methods(Static)
		function R = spreading_r_computation(source_area, plate_area, plate_thickness, k_source, R_plate)	
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
	end

	methods
		function [A, B] = lin_model_create(obj, T, pdev, tm_ver)

			% TODO: T
			% TODO: checks on inputs, etc. like Nv*Nh = Nc?
			
			lNc = obj.Nc;
			lNh = obj.Nh;
			lNv = obj.Nv;
			
			debug = 0;
			
			lAirPos = obj.air_pos;
			lMbPos = obj.mb_pos;
			lPcbPos = obj.pcb_pos;
			lAlPos = obj.al_pos;		

			North = obj.north_pos;
			East = obj.east_pos;
			South = obj.south_pos;
			West = obj.west_pos;

			Si = obj.si_pos;
			Cu = obj.cu_pos;
			
			chiplet_width = max( sum( squeeze(obj.RC_fp_dim(1+obj.extt_rows:end-obj.extb_rows,1+obj.extl_cols:end-obj.extr_cols, East)), 2 ) + ...
				sum( squeeze(obj.RC_fp_dim(1+obj.extt_rows:end-obj.extb_rows,1+obj.extl_cols:end-obj.extr_cols, West)), 2 ) );
			chiplet_length = max( sum( squeeze(obj.RC_fp_dim(1+obj.extt_rows:end-obj.extb_rows,1+obj.extl_cols:end-obj.extr_cols, North)), 2 ) + ...
				sum( squeeze(obj.RC_fp_dim(1+obj.extt_rows:end-obj.extb_rows,1+obj.extl_cols:end-obj.extr_cols, South)), 2 ) );
	
			%%% Components of non-fully modeled Layers
			A_comp = zeros(obj.add_states, 1);
			for i=1:obj.add_states
				A_comp(i) = obj.wid_comp(i) * obj.len_comp(i);
			end
			
			% Capacitances
			C_al = obj.c_al * A_comp(end-lAlPos) * obj.t_comp(end-lAlPos);
			C_pcb = obj.c_pcb * A_comp(end-lPcbPos) * obj.t_comp(end-lPcbPos);
			C_mb = obj.c_mb * A_comp(end-lMbPos) * obj.t_comp(end-lMbPos);
			C_air = obj.c_air * A_comp(end-lAirPos) * obj.t_comp(end-lAirPos);
			
			% Resistances
			%TODO: Move this to functions to make it less WOT
			%partial
			Ri_air_tot_v = obj.t_comp(end-lAirPos) / A_comp(end-lAirPos) / obj.k_air / obj.air_factor;
			air_al_t = obj.t_comp(end-lAirPos) - (obj.t_comp(end-lAlPos)+obj.t_cu+obj.t_si+obj.t_comp(end-lAirPos));
			Ri_air_al_t = air_al_t / A_comp(end-lAirPos) / obj.k_air / obj.air_factor;
			Ri_air_al_len = (obj.len_comp(end-lAirPos)-obj.len_comp(end-lAlPos)) / 2 ...
							/ (obj.wid_comp(end-lAirPos) * obj.t_comp(end-lAirPos)) ...
							/ obj.k_air / obj.air_factor;
			Ri_air_al_wid = (obj.wid_comp(end-lAirPos)-obj.wid_comp(end-lAlPos)) / 2 ...
							/ (obj.len_comp(end-lAirPos) * obj.t_comp(end-lAirPos)) ...
							/ obj.k_air / obj.air_factor;
			%
			Ri_al_v = obj.t_comp(end-lAlPos) / A_comp(end-lAlPos) / obj.k_al;
			Ri_al_len = obj.len_comp(end-lAlPos) / (obj.wid_comp(end-lAlPos) * obj.t_comp(end-lAlPos)) / obj.k_al;
			Ri_al_wid = obj.wid_comp(end-lAlPos) / (obj.len_comp(end-lAlPos) * obj.t_comp(end-lAlPos)) / obj.k_al;
			%
			Ri_pcb_v = obj.t_comp(end-lPcbPos) / A_comp(end-lPcbPos) / obj.k_pcb;
			Ri_mb_v = obj.t_comp(end-lMbPos) / A_comp(end-lMbPos) / obj.k_mb;
			%
			RaL = 1.948E+09 * air_al_t^3;
			h_alair_cond = RaL^(1/4) * 0.54 * obj.k_air / air_al_t;
			R_alair_cond = 1 / (h_alair_cond * A_comp(end-lAlPos) * obj.al_fins_coeff);
			
			RaL = 1.948E+09 * (air_al_t + obj.t_comp(end-lAlPos)/2)^3;
			h_alair_cond = obj.k_air / (air_al_t + obj.t_comp(end-lAlPos)/2) * ...
				(0.825 + (0.387 * RaL^(1/6))/(1+(0.492/7.039E-01)^(9/16))^(8/27))^2;
			R_alair_len_cond = 1 / (h_alair_cond * (obj.wid_comp(end-lAlPos) * obj.t_comp(end-lAlPos)));
			R_alair_wid_cond = 1 / (h_alair_cond * (obj.len_comp(end-lAlPos) * obj.t_comp(end-lAlPos)));
			
			RaL = 1.948E+09 * (obj.t_comp(end-lAirPos))^3;
			h_mbair_cond = RaL^(1/4) * 0.54 * obj.k_air / (obj.t_comp(end-lAirPos));
			R_mbair_cond = 1 / (h_mbair_cond * (A_comp(end-lMbPos) - A_comp(end-lAlPos)));
			
			%source_area, plate_area, plate_thickness, k_source, R_plate
			R_spr = obj.spreading_r_computation(A_comp(end-lAlPos), A_comp(end-lAirPos), ...
				obj.t_comp(end-lAirPos), obj.k_al, Ri_air_al_t);
			R_alair_v = Ri_al_v + R_alair_cond + R_spr;

			R_spr = obj.spreading_r_computation((obj.wid_comp(end-lAlPos) * obj.t_comp(end-lAlPos)), ...
				(obj.wid_comp(end-lAirPos) * obj.t_comp(end-lAirPos)), obj.len_comp(end-lAirPos), obj.k_al, Ri_air_al_len);
			R_alair_len = Ri_al_len + R_alair_len_cond + R_spr;
			
			R_spr = obj.spreading_r_computation(obj.len_comp(end-lAlPos) * obj.t_comp(end-lAlPos), ...
				(obj.len_comp(end-lAirPos) * obj.t_comp(end-lAirPos)), obj.wid_comp(end-lAirPos), obj.k_al, Ri_air_al_wid);
			R_alair_wid = Ri_al_wid + R_alair_wid_cond + R_spr;
			%
			R_spr = obj.spreading_r_computation(A_comp(end-lPcbPos), A_comp(end-lMbPos), ...
				obj.t_comp(end-lMbPos), obj.k_pcb, Ri_mb_v);
			R_pcbmb_v = Ri_pcb_v + Ri_mb_v + R_spr;
			R_spr = obj.spreading_r_computation(A_comp(end-lMbPos), A_comp(end-lAirPos), ...
				obj.t_comp(end-lAirPos), obj.k_mb, Ri_air_tot_v);
			R_mbair_v = Ri_mb_v + R_mbair_cond + R_spr;
			
			%TODO NOWWW
			%R_air_v = R_air_v * obj.ThMoNoise(fix(i/2)+1,7);
			%C_air = Ci_air * obj.ThMoNoise(fix(i/2)+1,7);

			A = zeros(obj.Ns);
			
			for i=1:2*lNc
				if (mod(i,2)>0)	%For each core:

					% Core Position
					ci = fix((i-1)/2)+1;
					irow = fix((ci-1)/obj.Nv)+1 + obj.extt_rows;
					icol = mod((ci-1),obj.Nv)+1 + obj.extl_cols;

					% Computing Core R 
					li = obj.RC_fp_dim(irow,icol, North) + obj.RC_fp_dim(irow,icol, South);
					wi = obj.RC_fp_dim(irow,icol, East) + obj.RC_fp_dim(irow,icol, West);
					%Si:
					Ri_si_hn = obj.RC_fp_dim(irow,icol, North)/(wi*obj.t_si)/obj.R_fp_material(irow,icol, Si);
					Ri_si_hs = obj.RC_fp_dim(irow,icol, South)/(wi*obj.t_si)/obj.R_fp_material(irow,icol, Si);
					Ri_si_he = obj.RC_fp_dim(irow,icol, East)/(li*obj.t_si)/obj.R_fp_material(irow,icol, Si);
					Ri_si_hw = obj.RC_fp_dim(irow,icol, West)/(li*obj.t_si)/obj.R_fp_material(irow,icol, Si);
					%
					Ri_si_v = (obj.t_si/2)/(wi*li)/obj.R_fp_material(irow,icol, Si);

					%Cu:
					Ri_cu_hn = obj.RC_fp_dim(irow,icol, North)/(wi*obj.t_cu)/obj.R_fp_material(irow,icol, Cu);
					Ri_cu_hs = obj.RC_fp_dim(irow,icol, South)/(wi*obj.t_cu)/obj.R_fp_material(irow,icol, Cu);
					Ri_cu_he = obj.RC_fp_dim(irow,icol, East)/(li*obj.t_cu)/obj.R_fp_material(irow,icol, Cu);
					Ri_cu_hw = obj.RC_fp_dim(irow,icol, West)/(li*obj.t_cu)/obj.R_fp_material(irow,icol, Cu);
					%
					Ri_cu_v = (obj.t_cu/2)/(wi*li)/obj.R_fp_material(irow,icol, Cu);

					%Others: Al + PCB
					R_spr_cual = obj.spreading_r_computation((wi*li), A_comp(end-lAlPos), ...
						obj.t_comp(end-lAlPos), obj.k_cu, Ri_al_v);
					R_spr_sipcb = obj.spreading_r_computation((wi*li), A_comp(end-lPcbPos), ...
						obj.t_comp(end-lPcbPos), obj.k_si, Ri_pcb_v);

					% Computing Core C
					Ci_si = obj.C_fp_material(irow, icol, Si) * li*wi*obj.t_si;
					Ci_cu = obj.C_fp_material(irow, icol, Cu) * li*wi*obj.t_cu;
					
					% Computing Neighbourhood R
					wi = obj.RC_fp_dim(irow-1,icol, East) + obj.RC_fp_dim(irow-1,icol, West);
					Rin_si_hn = obj.RC_fp_dim(irow-1,icol, South) / (wi*obj.t_si) / obj.R_fp_material(irow-1,icol, Si);
					Rin_cu_hn = obj.RC_fp_dim(irow-1,icol, South) / (wi*obj.t_cu) / obj.R_fp_material(irow-1,icol, Cu);
					wi = obj.RC_fp_dim(irow+1,icol, East) + obj.RC_fp_dim(irow+1,icol, West);
					Rin_si_hs = obj.RC_fp_dim(irow+1,icol, North) / (wi*obj.t_si) / obj.R_fp_material(irow+1,icol, Si);
					Rin_cu_hs = obj.RC_fp_dim(irow+1,icol, North) / (wi*obj.t_cu) / obj.R_fp_material(irow+1,icol, Cu);
					li = obj.RC_fp_dim(irow,icol+1, North) + obj.RC_fp_dim(irow,icol+1, South);
					Rin_si_he = obj.RC_fp_dim(irow,icol+1, West) / (li*obj.t_si) / obj.R_fp_material(irow,icol+1, Si);
					Rin_cu_he = obj.RC_fp_dim(irow,icol+1, West) / (li*obj.t_cu) / obj.R_fp_material(irow,icol+1, Cu);
					li = obj.RC_fp_dim(irow,icol-1, North) + obj.RC_fp_dim(irow,icol-1, South);
					Rin_si_hw = obj.RC_fp_dim(irow,icol-1, East) / (li*obj.t_si) / obj.R_fp_material(irow,icol-1, Si);
					Rin_cu_hw = obj.RC_fp_dim(irow,icol-1, East) / (li*obj.t_cu) / obj.R_fp_material(irow,icol-1, Cu);

					% Series
					R_si_hn = Ri_si_hn + Rin_si_hn;
					R_si_hs = Ri_si_hs + Rin_si_hs;
					R_si_he = Ri_si_he + Rin_si_he;
					R_si_hw = Ri_si_hw + Rin_si_hw;

					R_sicu_v = Ri_si_v + Ri_cu_v + obj.R_TIM1;

					R_cu_hn = Ri_cu_hn + Rin_cu_hn;
					R_cu_hs = Ri_cu_hs + Rin_cu_hs;
					R_cu_he = Ri_cu_he + Rin_cu_he;
					R_cu_hw = Ri_cu_hw + Rin_cu_hw;					
					
					R_cual_v = Ri_cu_v + obj.R_TIM2 + Ri_al_v + R_spr_cual;
					R_sipcb_v = Ri_si_v + Ri_pcb_v + R_spr_sipcb;
					
					C_si = Ci_si;
					C_cu = Ci_cu;

					% Noise
					if pdev==1
					R_si_hn = R_si_hn * obj.obj.param_dev_per(ci,1);
					R_si_hs = R_si_hs * obj.obj.param_dev_per(ci,1);
					R_si_he = R_si_he * obj.obj.param_dev_per(ci,1);
					R_si_hw = R_si_hw * obj.obj.param_dev_per(ci,1);

					R_cu_hn = R_cu_hn * obj.obj.param_dev_per(ci,2);
					R_cu_hs = R_cu_hs * obj.obj.param_dev_per(ci,2);
					R_cu_he = R_cu_he * obj.obj.param_dev_per(ci,2);
					R_cu_hw = R_cu_hw * obj.obj.param_dev_per(ci,2);

					R_sicu_v = R_sicu_v * obj.obj.param_dev_per(ci,3);
					R_cual_v = R_cual_v * obj.obj.param_dev_per(ci,4);

					C_si = Ci_si * obj.obj.param_dev_per(ci,5);
					C_cu = Ci_cu * obj.obj.param_dev_per(ci,6);
					end
					
					switch tm_ver
						case 10
							% R_si_v=1.6*2; %K/W silicon vertical (towards heat spreader) thermal resistance
							% R_cu_v=290; % K/W copper vertical (towards ambient) thermal resistance
							% R_si_h=31.9;%22.9; % K/W silicon horizontal (towards neighbours) thermal resistance
							% R_cu_h=1.2; % % K/W copper horizontal (towards neighbours) thermal resistance
							R_sicu_v=5; %K/W silicon vertical (towards heat spreader) thermal resistance
							R_cual_v=290; % K/W copper vertical (towards ambient) thermal resistance
							R_si_h=42;%22.9; % K/W silicon horizontal (towards neighbours) thermal resistance
							R_si_hn = R_si_h;
							R_si_hs = R_si_h;
							R_si_he = R_si_h;
							R_si_hw = R_si_h;
							R_cu_h=2*1.2; % % K/W copper horizontal (towards neighbours) thermal resistance
							R_cu_hn = R_cu_h;
							R_cu_hs = R_cu_h;
							R_cu_he = R_cu_h;
							R_cu_hw = R_cu_h;
							
							C_si=1e-3; % Silicon part thermal capacitance
							C_cu=1.2e-2; % copper part thermal capacitance

							% Noise
							if pdev==1
							R_si_hn = R_si_hn * obj.obj.param_dev_per(ci,1);
							R_si_hs = R_si_hs * obj.obj.param_dev_per(ci,1);
							R_si_he = R_si_he * obj.obj.param_dev_per(ci,1);
							R_si_hw = R_si_hw * obj.obj.param_dev_per(ci,1);
		
							R_cu_hn = R_cu_hn * obj.obj.param_dev_per(ci,2);
							R_cu_hs = R_cu_hs * obj.obj.param_dev_per(ci,2);
							R_cu_he = R_cu_he * obj.obj.param_dev_per(ci,2);
							R_cu_hw = R_cu_hw * obj.obj.param_dev_per(ci,2);
		
							R_sicu_v = R_sicu_v * obj.obj.param_dev_per(ci,3);
							R_cual_v = R_cual_v * obj.obj.param_dev_per(ci,4);
		
							C_si = Ci_si * obj.obj.param_dev_per(ci,5);
							C_cu = Ci_cu * obj.obj.param_dev_per(ci,6);
							end

						case 1				
							obj.pw2therm_coeff = obj.pw2therm_coeff * 0.85;
							
							C_si = C_si*1.1;
							C_cu = C_cu*1.1;

							g1=obj.gaussian_filter(max(obj.Nv, obj.Nh),2.5)*16;
							%g1 = g1 + (1-0.025 - g1(1));
							g1 = 1./g1 * 0.75;						
							
							R_si_hn = R_si_hn * g1(irow,icol); % * 2.67;
							R_si_hs = R_si_hs * g1(irow,icol); % * 2.67;
							R_si_he = R_si_he * g1(irow,icol); % * 2.67;
							R_si_hw = R_si_hw * g1(irow,icol); % * 2.67;							

							R_sicu_v = R_sicu_v / 2.2; %1.5;

							R_cu_hn = R_cu_hn * g1(irow,icol) * 2; % * 4;
							R_cu_hs = R_cu_hs * g1(irow,icol) * 2; % * 4;
							R_cu_he = R_cu_he * g1(irow,icol) * 2; % * 4;
							R_cu_hw = R_cu_hw * g1(irow,icol) * 2; % * 4;

							%R_cual_v = R_cual_v; %* 1.5;% * 1.3;							
							%R_sipcb_v = R_sipcb_v; % * 1.5;
							%R_pcbmb_v = R_pcbmb_v;						
							
						case 2 
							C_si = C_si*1.1;
							C_cu = C_cu*1.1;
							
							%pw2therm_coeff = 0.9;

							g1=obj.gaussian_filter(max(obj.Nv, obj.Nh),2.5)*16;
							%g1 = g1 + (1-0.025 - g1(1));
							g1 = 1./g1 * 0.7;						
							
							R_si_hn = R_si_hn * g1(irow,icol); % * 2.67;
							R_si_hs = R_si_hs * g1(irow,icol); % * 2.67;
							R_si_he = R_si_he * g1(irow,icol); % * 2.67;
							R_si_hw = R_si_hw * g1(irow,icol); % * 2.67;							

							R_sicu_v = R_sicu_v / 1.2;

							R_cu_hn = R_cu_hn * g1(irow,icol) * 2; % * 4;
							R_cu_hs = R_cu_hs * g1(irow,icol) * 2; % * 4;
							R_cu_he = R_cu_he * g1(irow,icol) * 2; % * 4;
							R_cu_hw = R_cu_hw * g1(irow,icol) * 2; % * 4;

							R_cual_v = R_cual_v / 1.5;% * 1.3;							
							R_sipcb_v = R_sipcb_v / 1.5;
							R_pcbmb_v = R_pcbmb_v / 1.07;
							
							air_t_couple = 0.2; %exp(obj.Ts / (C_al*Ri_al_v) ) - 1; %0.4;
					end
					
					if debug
						clc;
						drow = fix((ci-1)/lNv)+1;
						dcol = mod((ci-1),lNh)+1;
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
					A(i,i) = -1/(R_sicu_v*C_si) -1/C_si * (1/R_si_he+1/R_si_hn+1/R_si_hw+1/R_si_hs) - 1/(C_si*R_sipcb_v) * obj.si_pcb_fact;
					
					%PCB:
					A(i, end-lPcbPos) = 1/(C_si*R_sipcb_v) * obj.si_pcb_fact;
					A(end-lPcbPos, i) = 1/(C_pcb*R_sipcb_v) * obj.si_pcb_fact;
				% Heat-Spread (Cooper dyn)
				else %if (mod(i,2)==0)
					A(i,i) = -1/(R_sicu_v*C_cu) -1/(R_cual_v*C_cu) -1/C_cu * (1/R_cu_he+1/R_cu_hn+1/R_cu_hw+1/R_cu_hs);
					
					%heat sink:
					A(i, end-lAlPos) = 1/(C_cu*R_cual_v);
					A(end-lAlPos, i) = 1/(C_al*R_cual_v);
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
						A(i,end-lAirPos)=1/C_cu * (1/R_cu_hn + 1/R_cu_hw);
						A(end-lAirPos,i)=1/C_air * (1/R_cu_hn + 1/R_cu_hw);
					elseif (i==2*lNv)
						A(i,end-lAirPos)=1/C_cu * (1/R_cu_hn + 1/R_cu_he);
						A(end-lAirPos,i)=1/C_air * (1/R_cu_hn + 1/R_cu_he);
					elseif (i==(lNh-1)*2*lNv+2) 
						A(i,end-lAirPos)=1/C_cu * (1/R_cu_hs + 1/R_cu_hw);
						A(end-lAirPos,i)=1/C_air * (1/R_cu_hs + 1/R_cu_hw);
					elseif (i==2*lNc)
						A(i,end-lAirPos)=1/C_cu * (1/R_cu_hs + 1/R_cu_he);
						A(end-lAirPos,i)=1/C_air * (1/R_cu_hs + 1/R_cu_he);
						
					% diagonal terms for vertical (west and east) boundary nodes
					elseif ((mod(i,2*lNv)==2)&&(i~=(lNh-1)*2*lNv+2)&&(i~=2))
						A(i,end-lAirPos)=1/(C_cu*R_cu_hw);
						A(end-lAirPos,i)=1/(C_air*R_cu_hw);
					elseif ((mod(i,2*lNv)==0) && (i~=2*lNc) && (i~=2*lNv))
						A(i,end-lAirPos)=1/(C_cu*R_cu_he);
						A(end-lAirPos,i)=1/(C_air*R_cu_he);
					% diagonal terms for horizontal (upper and lower) boundary nodes
					elseif ((i~=2)&& (i<2*lNv))
						A(i,end-lAirPos)=1/(C_cu*R_cu_hn);
						A(end-lAirPos,i)=1/(C_air*R_cu_hn);
					elseif ((i>(lNh-1)*2*lNv+2) && (i<2*lNc))
						A(i,end-lAirPos)=1/(C_cu*R_cu_hs);
						A(end-lAirPos,i)=1/(C_air*R_cu_hs);
						
					% diagonal terms corresponding to other nodes (internal)
					%else
						%A(i,end)=0;
						%A(end,i)=0;
					end
					
				% Silicium
				else
					% diagonal terms for vertexes
					if (i==1)
						A(i,end-lAirPos)=1/C_si * (1/R_si_hn + 1/R_si_hw);
						A(end-lAirPos,i)=1/C_air * (1/R_si_hn + 1/R_si_hw);					
					elseif (i==2*lNv-1)
						A(i,end-lAirPos)=1/C_si * (1/R_si_hn + 1/R_si_he);
						A(end-lAirPos,i)=1/C_air * (1/R_si_hn + 1/R_si_he);	
					elseif (i==(lNh-1)*2*lNv+1)
						A(i,end-lAirPos)=1/C_si * (1/R_si_hs + 1/R_si_hw);
						A(end-lAirPos,i)=1/C_air * (1/R_si_hs + 1/R_si_hw);						
					elseif (i==2*lNc-1)
						A(i,end-lAirPos)=1/C_si * (1/R_si_hs + 1/R_si_he);
						A(end-lAirPos,i)=1/C_air * (1/R_si_hs + 1/R_si_he);	

					% diagonal terms for vertical (west and east) boundary nodes
					elseif (mod(i,2*lNv)==1)&&(i~=(lNh-1)*2*lNv+1)&&(i~=1)
						A(i,end-lAirPos)=1/(C_si*R_si_hw);
						A(end-lAirPos,i)=1/(C_air*R_si_hw);
					elseif (mod(i+1,2*lNv)==0) && (i~=2*lNc-1) && (i~=2*lNv-1)
						A(i,end-lAirPos)=1/(C_si*R_si_he);
						A(end-lAirPos,i)=1/(C_air*R_si_he);
					% diagonal terms for horizontal (upper and lower) boundary nodes
					elseif (i~=1)&& (i<2*lNv-2)
						A(i,end-lAirPos)=1/(C_si*R_si_hn);
						A(end-lAirPos,i)=1/(C_air*R_si_hn);
					elseif (i>(lNh-1)*2*lNv+1) && (i<2*lNc-1)
						A(i,end-lAirPos)=1/(C_si*R_si_hs);
						A(end-lAirPos,i)=1/(C_air*R_si_hs);

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
					if (tm_ver == 2)
						%for j=1:lNc
						j=i/2; %remember
						row = fix((j-1)/lNv)+1;
						%col = mod((j-1),lNh)+1;

						%1: Add/change the final input
							A(2*j,end-lAirPos) = A(2*j,end-lAirPos)-1/(R_cual_v*C_cu) + ...
								(1/(R_cual_v*C_cu))*(1-air_t_couple)^(row-1);
						%2: add previous rows part
						for k=1:row-1
							id2 = j - (obj.Nv*(row-k));
							A(2*j, 2*id2) = A(2*j, 2*id2) + ...
								air_t_couple * (1-air_t_couple)^(row-1-k) * ...
								1/(R_cual_v*C_cu);
						end
						%end
						
						%3: fix self
						%This is not needed beacuse all coefficients stay
						%the same!
						A(2*j, 2*j) = -(sum(A(2*j,:)) - A(2*j, 2*j));
						
					end	%(tm_ver == 2)
				end %mod(i,2)==0

			end % for i

			%HeatSink
			A(end-lAlPos, end-lAirPos) = 1/C_al * (1/R_alair_v + 2/R_alair_len + 2/R_alair_wid) + heatsink_fan_dis/C_al;
			A(end-lAirPos, end-lAlPos) = 1/C_air * (1/R_alair_v + 2/R_alair_len + 2/R_alair_wid) + heatsink_fan_dis/C_al; %here C_al or C_air?
			A(end-lAlPos, end-lAlPos) = 0;
			A(end-lAlPos, end-lAlPos) = -sum(A(end-lAlPos,:));
			
			%PCB
			%negligible with air
			A(end-lPcbPos, end-lMbPos) =1/(C_pcb*R_pcbmb_v) * obj.pcb_mb_fact;
			A(end-lMbPos, end-lPcbPos) =1/(C_mb*R_pcbmb_v) * obj.ppcb_mb_fact;
			A(end-lPcbPos, end-lPcbPos) = 0;
			A(end-lPcbPos, end-lPcbPos) = -sum(A(end-lPcbPos,:));
			
			%Motherboard
			A(end-lMbPos, end-lAirPos) = 1/(C_mb * R_mbair_v) * obj.pmb_air_fact;
			A(end-lAirPos, end-lMbPos) = 1/(C_air * R_mbair_v) * obj.pmb_air_fact;
			A(end-lMbPos, end-lMbPos) = 0;
			A(end-lMbPos, end-lMbPos) = -sum(A(end-lMbPos,:));
			
			%Air
			%A(end, end) = -(lNc+(lNh-2)*2+(lNv-2)*2+4*2)/(R_cu_v*C_cu) -100*lNc/4; % -0.149075;
			A(end-lAirPos, end-lAirPos)=0;
			A(end-lAirPos, end-lAirPos) = -sum(A(end-lAirPos,:)) - obj.case_fan_dis; % -0.149075;
			
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
					irow = fix((ci-1)/obj.Nv)+1 ;
					icol = mod((ci-1),obj.Nv)+1;
					
					li = obj.CPw_fp_dim(irow,icol, North) + obj.CPw_fp_dim(irow,icol, South);
					wi = obj.CPw_fp_dim(irow,icol, East) + obj.CPw_fp_dim(irow,icol, West);
					C_core = obj.C_fp_material(irow, icol, Si) * li*wi*obj.t_si;					
					if d==1
						C_core = C_core * obj.param_dev_per(ci,5);
					end
					
					if (tm_ver == 1) %|| (tm_ver == 2)
						g2=obj.gaussian_filter(max(obj.Nv, obj.Nh),3)*16;
						g2 = g2 + (1-0.01 - g2(1));
						C_core = C_core / g2(irow-obj.extt_rows,icol-obj.extl_cols);% / 1.5;
					end							
							
					B(i,k) = 1/C_core * obj.pw2therm_coeff;
				end
			end

			B(end-lAirPos,:) = [zeros(1, lNc) obj.case_fan_dis/obj.case_fan_nom_speed];

		end %function
		function obj = create_thermal_model(obj)
			
			tm_ver = obj.model_ver;
			
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
			
			[obj.Ac_nom, obj.Bc_nom] = obj.lin_model_create(obj.temp_amb, 0, tm_ver);
			[obj.Ac_true, obj.Bc_true] = obj.lin_model_create(obj.temp_amb, 1, tm_ver);
			
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
		function obj = create_thermal_model_noise(obj)
			ddiv = 30;
			dmean = 1.01;
			obj.param_dev_per = (randn(obj.Nc, 8) / ddiv + dmean);
		end
		function obj = default_floorplan_config(obj)
			% create Distance/position matrix

			North = obj.north_pos;
			East = obj.east_pos;
			South = obj.south_pos;
			West = obj.west_pos;

			Si = obj.si_pos;
			Cu = obj.cu_pos;

			core_cols = obj.Nv;
			core_rows = obj.Nh;

			cols = core_cols + obj.extl_cols + obj.extr_cols;
			rows = core_rows + obj.extt_rows + obj.extb_rows;

			%TODO instead of initial values, add air
			obj.RC_fp_dim		= zeros(cols, rows, 4);
			obj.CPw_fp_dim		= zeros(core_cols, core_rows, 4);
			obj.R_fp_material	= zeros(cols, rows, obj.full_model_layers);
			obj.C_fp_material	= zeros(cols, rows, obj.full_model_layers);

			% Populating:
			obj.RC_fp_dim([1:obj.extt_rows, end-obj.extb_rows+1:end])		= 0.05;
			obj.R_fp_material([1:obj.extt_rows, end-obj.extb_rows+1:end])	= 1e-16;
			obj.C_fp_material([1:obj.extt_rows, end-obj.extb_rows+1:end])	= 0.025;

			for r = (1+obj.extt_rows):(rows - obj.extb_rows)
				for c = (1+obj.extl_cols):(cols - obj.extr_cols)

					obj.RC_fp_dim(r,c,North) = 2.30e-3;
					obj.RC_fp_dim(r,c,South) = 1e-3;

					obj.RC_fp_dim(r,c,West) = 0.75e-3;
					obj.RC_fp_dim(r,c,East) = 0.75e-3;					

					% Internal Mesh connection
					if (mod(c,2) == 1) && (c~=(core_cols+obj.extl_cols))
						obj.RC_fp_dim(r,c,East) = obj.RC_fp_dim(r,c,East) + 0.25e-3;
					end
					
					% Eternal Mesh connection: row
					if r == (1+obj.extt_rows)
						obj.RC_fp_dim(r,c,North) = obj.RC_fp_dim(r,c,North) + 0.3e-3;
					end
					if r == (rows - obj.extb_rows)
						obj.RC_fp_dim(r,c,South) = obj.RC_fp_dim(r,c,South) + 0.3e-3;
					end
					% Eternal Mesh connection: cols
					if c == (1+obj.extl_cols)
						obj.RC_fp_dim(r,c,West) = obj.RC_fp_dim(r,c,West) + 0.3e-3;
					end
					if c == (cols - obj.extr_cols)
						obj.RC_fp_dim(r,c,East) = obj.RC_fp_dim(r,c,East) + 0.3e-3;
					end
					
					obj.R_fp_material(r,c,Si) = obj.k_si;
					obj.R_fp_material(r,c,Cu) = obj.k_cu;
					obj.C_fp_material(r,c,Si) = obj.c_si;
					obj.C_fp_material(r,c,Cu) = obj.c_cu;
				end	
			end
			for r = 1:rows
				for c = 1:cols
					obj.CPw_fp_dim(r,c,North) = 1e-3;
					obj.CPw_fp_dim(r,c,South) = 1e-3;
					obj.CPw_fp_dim(r,c,West) = 0.75e-3;
					obj.CPw_fp_dim(r,c,East) = 0.75e-3;
				end
			end

		end
	end %Methods

	%% Dependent Variables
	methods
		% State Dimension
		function value = get.Ns(obj)
			value = obj.Nc*obj.full_model_layers+obj.add_states;
		end
		% Input Dimension
		function value = get.Ni(obj)
			value = obj.Nc+obj.add_inputs;
		end
		% Outputs
		function value = get.Nout(obj)
			value = obj.Nc+obj.add_outputs;
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

		%% SET
		function obj = set.model_deviations(obj, val)
			%check modifications
			if val
				if val ~= obj.model_deviations 
					obj = obj.create_thermal_model_noise();
				end
			else
				%reset: 
				obj.param_dev_per = ones(obj.Nc, obj.param_dev_dim2);					
			end
			obj.model_deviations = val;
		end

		% Matrices
		function obj = set.Ac_nom(obj, val)
			cmp = [obj.Ns obj.Ns];
			if all(size(val) == cmp) && ~isempty(val)
				obj.Ac_nom = val;
			else
				warning("[TM Error]: Ac_nom wrong input size.");
				disp(strcat("[TM] Ac_nom required size: [", num2str(cmp(1)), ",", num2str(cmp(2)),"]"));
			end
		end
		function obj = set.Ac_true(obj, val)
			cmp = [obj.Ns obj.Ns];
			if all(size(val) == cmp)  && ~isempty(val)
				obj.Ac_true = val;
			else
				warning("[TM Error]: Ac_true wrong input size.");
				disp(strcat("[TM] Ac_true required size: [", num2str(cmp(1)), ",", num2str(cmp(2)),"]"));
			end
		end
		function obj = set.Bc_nom(obj, val)
			cmp = [obj.Ns obj.Ni];
			if all(size(val) == cmp) && ~isempty(val)
				obj.Bc_nom = val;
			else
				warning("[TM Error]: Bc_nom wrong input size.");
				disp(strcat("[TM] Bc_nom required size: [", num2str(cmp(1)), ",", num2str(cmp(2)),"]"));
			end
		end
		function obj = set.Bc_true(obj, val)
			cmp = [obj.Ns obj.Ni];
			if all(size(val) == cmp) && ~isempty(val)
				obj.Bc_true = val;
			else
				warning("[TM Error]: Bc_true wrong input size.");
				disp(strcat("[TM] Bc_true required size: [", num2str(cmp(1)), ",", num2str(cmp(2)),"]"));
			end
		end
		function obj = set.C(obj, val)
			cmp = [obj.Nout obj.Ns];
			if all(size(val) == cmp) && ~isempty(val)
				obj.C = val;
			else
				warning("[TM Error]: C wrong input size.");
				disp(strcat("[TM] C required size: [", num2str(cmp(1)), ",", num2str(cmp(2)),"]"));
			end
		end
		function obj = set.D(obj, val)
			cmp = [obj.Nout obj.Ni];
			if all(size(val) == cmp) && ~isempty(val)
				obj.D = val;
			else
				warning("[TM Error]: D wrong input size.");
				disp(strcat("[TM] D required size: [", num2str(cmp(1)), ",", num2str(cmp(2)),"]"));
			end
		end
		
		function obj = set.param_dev_per(obj, val)
			cmp = [obj.Nc obj.param_dev_dim2];
			if all(size(val) == cmp) && ~isempty(val)
				obj.param_dev_per = val;
			else
				warning("[TM Error]: param_dev_per wrong input size.");
				disp(strcat("[TM] param_dev_per required size: [", num2str(cmp(1)), ",", num2str(cmp(2)),"]"));
			end
		end

		function obj = set.sensor_active(obj, val)
			cmp = [obj.add_states 1];
			if all(size(val) == cmp) && ~isempty(val)
				obj.sensor_active = val;
			else
				warning("[TM Error]: sensor_active wrong input size.");
				disp(strcat("[TM] sensor_active required size: [", num2str(cmp(1)), ",", num2str(cmp(2)),"]"));
			end
		end

		function obj = set.RC_fp_dim(obj, val)
			cols = obj.Nv + obj.extl_cols + obj.extr_cols;
			rows = obj.Nh + obj.extt_rows + obj.extb_rows;
			cmp = [cols, rows, 4];
			if all(size(val) == cmp) && ~isempty(val)
				obj.RC_fp_dim = val;
			else
				warning("[TM Error]: RC_fp_dim wrong input size.");
				disp(strcat("[TM] RC_fp_dim required size: [", num2str(cmp(1)), ",", num2str(cmp(2)), ",", num2str(cmp(3)), "]"));
			end
		end

		function obj = set.CPw_fp_dim(obj, val)
			cols = obj.Nv;
			rows = obj.Nh;
			cmp = [cols, rows, 4];
			if all(size(val) == cmp) && ~isempty(val)
				obj.CPw_fp_dim = val;
			else
				warning("[TM Error]: CPw_fp_dim wrong input size.");
				disp(strcat("[TM] CPw_fp_dim required size: [", num2str(cmp(1)), ",", num2str(cmp(2)), ",", num2str(cmp(3)), "]"));
			end
		end

		function obj = set.R_fp_material(obj, val)
			cols = obj.Nv + obj.extl_cols + obj.extr_cols;
			rows = obj.Nh + obj.extt_rows + obj.extb_rows;
			cmp = [cols, rows, obj.full_model_layers];
			if all(size(val) == cmp) && ~isempty(val)
				obj.R_fp_material = val;
			else
				warning("[TM Error]: R_fp_material wrong input size.");
				disp(strcat("[TM] R_fp_material required size: [", num2str(cmp(1)), ",", num2str(cmp(2)), ",", num2str(cmp(3)), "]"));
			end
		end

		function obj = set.C_fp_material(obj, val)
			cols = obj.Nv + obj.extl_cols + obj.extr_cols;
			rows = obj.Nh + obj.extt_rows + obj.extb_rows;
			cmp = [cols, rows, obj.full_model_layers];
			if all(size(val) == cmp) && ~isempty(val)
				obj.C_fp_material = val;
			else
				warning("[TM Error]: C_fp_material wrong input size.");
				disp(strcat("[TM] C_fp_material required size: [", num2str(cmp(1)), ",", num2str(cmp(2)), ",", num2str(cmp(3)), "]"));
			end
		end
		


	end


end