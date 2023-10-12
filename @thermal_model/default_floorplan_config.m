function obj = default_floorplan_config(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

	% create Distance/position matrix
	obj.RC_fp_dim		= zeros(rows, cols, 4);
	obj.CPw_fp_dim		= zeros(core_rows, core_cols, 4);
	obj.R_fp_material	= zeros(rows, cols, obj.full_model_layers);
	obj.C_fp_material	= zeros(rows, cols, obj.full_model_layers);

	% Populating:
	%([1:obj.extt_rows, end-obj.extb_rows+1:end], ...
	%	[1:obj.extl_cols, end-obj.extr_cols+1:end],:)
	obj.RC_fp_dim(:,:,:)		= 0.05;
	obj.R_fp_material(:,:,:)	= 1e-16;
	obj.C_fp_material(:,:,:)	= 0.025;

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
	for r = 1:core_rows
		for c = 1:core_cols
			obj.CPw_fp_dim(r,c,North) = 1e-3;
			obj.CPw_fp_dim(r,c,South) = 1e-3;
			obj.CPw_fp_dim(r,c,West) = 0.75e-3;
			obj.CPw_fp_dim(r,c,East) = 0.75e-3;
		end
	end


end %function

