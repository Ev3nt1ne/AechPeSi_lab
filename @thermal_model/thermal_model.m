classdef thermal_model < handle
	%THERMAL_MODEL Summary of this class goes here
	%   Detailed explanation goes here
	
	properties		
		%% System Structure
		Nc (1,1) {mustBePositive, mustBeInteger} ...
			= 12;		% Number of Cores
		Nh (1,1) {mustBePositive, mustBeInteger} ...
			= 3;		% Number of rows
		Nv (1,1) {mustBePositive, mustBeInteger} ...
			= 4;		% Number of cols
		epos {mustBeNonnegative, mustBeNumeric, mustBeFinite} ... %empty positions
			= [];

		%% Matlab/Simulation
		model_ver (1,1) {mustBeNonnegative, mustBeInteger} ...
			= 0;				% Version of the Thermal Model
		%model_deviations (1,1) {mustBeNonnegative, mustBeNumericOrLogical} ...
		%	= 1;				% If the true system should be != from nominal
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

		sensors_active {mustBeNonnegative, mustBeNumericOrLogical, mustBeVector} ...
			= [0 0 0 0];

end

	properties(Dependent)
		%% Matlab/Simulation
		Ns;							% Number of States
		Ni;							% Number of inputs
		Nout;						% Number of "Outputs"
		add_outputs;

		% Matlab Matrices
		Ad_nom;
		Bd_nom;
		Ad_true;
		Bd_true;
	end

	properties(SetAccess=protected, GetAccess=public)
		observable = 0;
		
		Cc {mustBeNumeric, mustBeNonempty, mustBeFinite} ...
			= [1];

		var_changed = 1;
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
			obj.sensors_active = zeros(obj.add_states,1);
			obj.default_floorplan_config();

			% Post Initialization %%
			% Any code, including access to object
			obj.create_model_deviation();
			obj.model_init();
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
		function obj = anteSimCheckTM(obj)
			if obj.var_changed
				% Check if the dimensions are consistent
				if obj.Nc ~= (obj.Nh*obj.Nv - length(obj.epos))
					error("[TM] The dimensions of Nc, Nh, Nv, epos are not consistent for a model");
				else
					% Check if the floorplan has been updated
					lastwarn('');
					obj.RC_fp_dim = obj.RC_fp_dim;
					obj.CPw_fp_dim = obj.CPw_fp_dim;
					obj.R_fp_material = obj.R_fp_material;
					obj.C_fp_material = obj.C_fp_material;
					[msgstr, ~] = lastwarn;
					if ~isempty(msgstr)
						warning("[TM] The floorplan was not updated correctly. Please update RC_fp_dim, CPw_fp_dim, R_fp_material, C_fp_material or launch 'default_floorplan_config()'.");

						prompt = {['Do you want to run default_floorplan_config()? (1=yes, 0=no)' newline 'ATTENTION! This will OVERWRITE previous values.']};
						dlgtitle = '[TM] Thermal Model';
						fieldsize = [1 45];
						definput = {'1'};
						usrin = inputdlg(prompt,dlgtitle,fieldsize,definput);
						%usrin = input(['[TM] Do you want to run default_floorplan_config()? (1=yes, 0=no)' newline 'ATTENTION! This will overwrite previous values.']);
						if usrin{1}=='1'
							obj.default_floorplan_config();
						end
					end
					if size(obj.param_dev_per,1) ~= obj.Nc
						warning("[TM] The parameter deviation vector 'param_dev_per' was not updated after Nc was changed. Fix it or call 'create_model_deviation'.");
						
						prompt = {['Do you want to run create_model_deviation()? (1=yes, 0=no)' newline 'ATTENTION! This will OVERWRITE previous values.']};
						dlgtitle = '[TM] Thermal Model';
						fieldsize = [1 45];
						definput = {'1'};
						usrin = inputdlg(prompt,dlgtitle,fieldsize,definput);
						if usrin{1}=='1'
							obj.create_model_deviation();
						end
					end

					warning("[TM] Some parametric values of the Thermal model changed without running model_init()!")
					prompt = {['Do you want to run model_init()? (1=yes, 0=no)' newline 'ATTENTION! This will OVERWRITE previous values.']};
					dlgtitle = '[TM] Thermal Model';
					fieldsize = [1 45];
					definput = {'1'};
					usrin = inputdlg(prompt,dlgtitle,fieldsize,definput);
					if usrin{1} == '1'
						obj.model_init();
					else
						obj.var_changed = 0;
					end
				end
			end			
		end
		function obj = model_init(obj, tm_ver)
			
			if nargin < 2
				tm_ver = obj.model_ver;
			end

			obj.var_changed = 0;
			
			% ==============
			% C & D
			% ==============
			obj.C = zeros(obj.Nout,obj.Ns);
			fml = obj.full_model_layers;
			k=-fml;
			for i=1:obj.Nc
				k=k+fml;
				obj.C(i,k+1)=1;
			end
			% to consider the case obj.add_outputs = 0
			i = obj.Nc+1;
			cadd = sum(diag(obj.sensors_active),2) > 0;
			j=1;
			while (i <= obj.Nc+obj.add_outputs)
				while(cadd(j) ~= 1)
					j=j+1;
				end
				obj.C(i, end-obj.add_states+j) = 1;
				i=i+1;
				j=j+1;
			end
			obj.D = zeros(obj.Nout, obj.Ni);

			obj.Cc = obj.C(1:obj.Nc,:);
			
			[obj.Ac_nom, obj.Bc_nom] = obj.create_model(obj.temp_amb, 0, tm_ver);
			[obj.Ac_true, obj.Bc_true] = obj.create_model(obj.temp_amb, 1, tm_ver);
			
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
		function obj = create_model_deviation(obj)
			ddiv = 30;
			dmean = 1.01;
			obj.param_dev_per = (randn(obj.Nc, 8) / ddiv + dmean);
		end
		obj = default_floorplan_config(obj);
	end
	methods(Access=protected)
		[A, B] = create_model(obj, T, pdev, tm_ver);
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
		function value = get.add_outputs(obj)
			value = sum(obj.sensors_active>0);
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
		%{
		function obj = set.model_deviations(obj, val)
			%check modifications
			if val ~= obj.model_deviations 
				if val				
					obj = obj.create_model_deviation();
				else
					%reset: 
					obj.param_dev_per = ones(obj.Nc, obj.param_dev_dim2);					
				end
			end
			obj.model_deviations = val;
		end
		%}

		% Matrices & Vectors
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
		function obj = set.Cc(obj, val)
			cmp = [obj.Nc obj.Ns];
			if all(size(val) == cmp) && ~isempty(val)
				obj.Cc = val;
			else
				warning("[TM Error]: Cc wrong input size.");
				disp(strcat("[TM] Cc required size: [", num2str(cmp(1)), ",", num2str(cmp(2)),"]"));
			end
		end
		%
		function obj = set.param_dev_per(obj, val)
			cmp = [obj.Nc obj.param_dev_dim2];
			if all(size(val) == cmp) && ~isempty(val)
				obj.param_dev_per = val;
			else
				warning("[TM Error]: param_dev_per wrong input size.");
				disp(strcat("[TM] param_dev_per required size: [", num2str(cmp(1)), ",", num2str(cmp(2)),"]"));
			end
		end

		function obj = set.sensors_active(obj, val)
			cmp = [1 obj.add_states];
			if length(val) > size(val,2)
				val = val';
			end
			if all(size(val) == cmp) && ~isempty(val)
				obj.sensors_active = val;
			else
				warning("[TM Error]: sensor_active wrong input size.");
				disp(strcat("[TM] sensor_active required size:", num2str(cmp(2))));
			end
		end

		function obj = set.t_comp(obj, val)
			cmp = [1 obj.add_states];
			if length(val) > size(val,2)
				val = val';
			end
			if all(size(val) == cmp) && ~isempty(val)
				obj.t_comp = val;
			else
				warning("[TM Error]: t_comp wrong input size.");
				disp(strcat("[TM] t_comp required size:", num2str(cmp(2))));
			end
		end

		function obj = set.wid_comp(obj, val)
			cmp = [1 obj.add_states];
			if length(val) > size(val,2)
				val = val';
			end
			if all(size(val) == cmp) && ~isempty(val)
				obj.wid_comp = val;
			else
				warning("[TM Error]: wid_comp wrong input size.");
				disp(strcat("[TM] wid_comp required size:", num2str(cmp(2))));
			end
		end

		function obj = set.len_comp(obj, val)
			cmp = [1 obj.add_states];
			if length(val) > size(val,2)
				val = val';
			end
			if all(size(val) == cmp) && ~isempty(val)
				obj.len_comp = val;
			else
				warning("[TM Error]: len_comp wrong input size.");
				disp(strcat("[TM] len_comp required size:", num2str(cmp(2))));
			end
		end

		function obj = set.RC_fp_dim(obj, val)
			cols = obj.Nv + obj.extl_cols + obj.extr_cols;
			rows = obj.Nh + obj.extt_rows + obj.extb_rows;
			cmp = [rows, cols, 4];
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
			cmp = [rows, cols, 4];
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
			cmp = [rows, cols, obj.full_model_layers];
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
			cmp = [rows, cols, obj.full_model_layers];
			if all(size(val) == cmp) && ~isempty(val)
				obj.C_fp_material = val;
			else
				warning("[TM Error]: C_fp_material wrong input size.");
				disp(strcat("[TM] C_fp_material required size: [", num2str(cmp(1)), ",", num2str(cmp(2)), ",", num2str(cmp(3)), "]"));
			end
		end

		%%%%

		function obj = set.Nc(obj, val)
			obj.var_changed = 1;
			disp("[TM] Remember to change Nh, Nv, epos, the parameter deviation vector, and the floorplan accordingly");
			obj.Nc = val;
		end
		function obj = set.Nh(obj, val)
			obj.var_changed = 1;
			disp("[TM] Remember to change Nv, epos, and the floorplan accordingly");
			obj.Nh = val;
		end
		function obj = set.Nv(obj, val)
			obj.var_changed = 1;
			disp("[TM] Remember to change Nh, epos, and the floorplan accordingly");
			obj.Nv = val;
		end
		function obj = set.epos(obj, val)			
			un = unique(val);
			if length(val) ~= length(un)
				warning("[Thermal] epos (empty positions) was not an unique vector. It has been automatically fixed");
			end
			obj.var_changed = 1;
			disp("[TM] Remember to change Nc, and the floorplan accordingly");
			obj.epos = un;
		end
		function obj = set.Ts(obj, val)
			obj.var_changed = 1;
			obj.Ts = val;
		end

	end %method

end %class