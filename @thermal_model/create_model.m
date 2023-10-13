function [A, B] = create_model(obj, T, pdev, tm_ver)

%CREATE_MODEL Summary of this function goes here
%   Detailed explanation goes here


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
	air_al_t = obj.t_comp(end-lAirPos) - (obj.t_comp(end-lAlPos)+obj.t_cu+obj.t_si+obj.t_comp(end-lPcbPos));
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
	
	%TODO
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
			R_si_hn = R_si_hn * obj.param_dev_per(ci,1);
			R_si_hs = R_si_hs * obj.param_dev_per(ci,1);
			R_si_he = R_si_he * obj.param_dev_per(ci,1);
			R_si_hw = R_si_hw * obj.param_dev_per(ci,1);

			R_cu_hn = R_cu_hn * obj.param_dev_per(ci,2);
			R_cu_hs = R_cu_hs * obj.param_dev_per(ci,2);
			R_cu_he = R_cu_he * obj.param_dev_per(ci,2);
			R_cu_hw = R_cu_hw * obj.param_dev_per(ci,2);

			R_sicu_v = R_sicu_v * obj.param_dev_per(ci,3);
			R_cual_v = R_cual_v * obj.param_dev_per(ci,4);

			C_si = Ci_si * obj.param_dev_per(ci,5);
			C_cu = Ci_cu * obj.param_dev_per(ci,6);
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
					R_si_hn = R_si_hn * obj.param_dev_per(ci,1);
					R_si_hs = R_si_hs * obj.param_dev_per(ci,1);
					R_si_he = R_si_he * obj.param_dev_per(ci,1);
					R_si_hw = R_si_hw * obj.param_dev_per(ci,1);

					R_cu_hn = R_cu_hn * obj.param_dev_per(ci,2);
					R_cu_hs = R_cu_hs * obj.param_dev_per(ci,2);
					R_cu_he = R_cu_he * obj.param_dev_per(ci,2);
					R_cu_hw = R_cu_hw * obj.param_dev_per(ci,2);

					R_sicu_v = R_sicu_v * obj.param_dev_per(ci,3);
					R_cual_v = R_cual_v * obj.param_dev_per(ci,4);

					C_si = Ci_si * obj.param_dev_per(ci,5);
					C_cu = Ci_cu * obj.param_dev_per(ci,6);
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
	A(end-lAlPos, end-lAirPos) = 1/C_al * (1/R_alair_v + 2/R_alair_len + 2/R_alair_wid) + obj.al_fan_dis/C_al;
	A(end-lAirPos, end-lAlPos) = 1/C_air * (1/R_alair_v + 2/R_alair_len + 2/R_alair_wid) + obj.al_fan_dis/C_al; %here C_al or C_air?
	A(end-lAlPos, end-lAlPos) = 0;
	A(end-lAlPos, end-lAlPos) = -sum(A(end-lAlPos,:));
	
	%PCB
	%negligible with air
	A(end-lPcbPos, end-lMbPos) =1/(C_pcb*R_pcbmb_v) * obj.pcb_mb_fact;
	A(end-lMbPos, end-lPcbPos) =1/(C_mb*R_pcbmb_v) * obj.pcb_mb_fact;
	A(end-lPcbPos, end-lPcbPos) = 0;
	A(end-lPcbPos, end-lPcbPos) = -sum(A(end-lPcbPos,:));
	
	%Motherboard
	A(end-lMbPos, end-lAirPos) = 1/(C_mb * R_mbair_v) * obj.mb_air_fact;
	A(end-lAirPos, end-lMbPos) = 1/(C_air * R_mbair_v) * obj.mb_air_fact;
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
			C_core = obj.C_fp_material(irow + obj.extt_rows, icol + obj.extl_cols, Si) * li*wi*obj.t_si;					
			if pdev==1
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

