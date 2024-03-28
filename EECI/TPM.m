
classdef TPM < HybridSystem & hpc_lab
    % A Thermal Power management modeled as a HybridSystem subclass.

    % Define variable properties that can be modified.
    properties
		F_target = 3.6;
		wl_target = [0 0 0.6 0.2 0.2];
		T_ctrl = 5e-1;
		p_budget = 40;
		%fan_speed = 1000;
		T_pid = 85 + 273.15;
		sT_case = 2.5; %steps
		sT_core = 2.5; %steps
		tT_case = 35+ 273.15; %threshold
		tT_core = 60+ 273.15;

		% PID Controller
		kp = 0.1; %1.9518;				% PID Kp
		ki = 1; %73.931;				% PID Ki
		aw_up = 0; %0.05;				% PID Anti-Windup Upper Saturation
		aw_down_c = -0.75;			% PID Anti-Windup Down Sat. Coefficient

		sat_up = 0;					% PID Upper Saturation
		
		pid_e_down = 2.5;			% PID banding, down constraint (absolute value)
		pid_e_up = 1;				% PID banding, up constraint (absolute value)
		pid_e_band_coeff = 0.7;		% PID banding coefficient

		disable_fan_boost = 0;
    end

    % Define constant properties that cannot be modified (i.e., "immutable").
    properties(SetAccess = immutable) 
		CFp = 1;
		HFp = 2;
		%qp = 1;
	end

	% Define constant properties that cannot be modified (i.e., "immutable").
    properties(SetAccess = private) 
		pid_integ = 0;
		mat_case;
		mat_al;
		mat_alpos;
		Fa;
		C_al_hs;
    end

    methods 
        function this = TPM()
            % Constructor for instances of the TPM class.

            % Call the constructor for the HybridSystem superclass and
            % pass it the state dimension. This is not strictly necessary, 
            % but it enables more error checking.
            state_dim = 22+9+3; %this.Ns+this.Nc+3;
            this = this@HybridSystem(state_dim);

			%init
			this = this.init();

			lAlPos = this.al_pos;
			lAirPos = this.air_pos;
			this.C_al_hs = this.c_al * this.wid_comp(end-lAlPos) * this.len_comp(end-lAlPos) * this.t_comp(end-lAlPos);

			this.mat_case = this.Ac_true(end) + this.case_fan_dis + this.al_fan_dis/this.C_al_hs;
			
			this.mat_al(1) = this.Ac_true(end-lAlPos, end-lAirPos) - this.al_fan_dis/this.C_al_hs;
			this.mat_alpos(1) = (this.Ns-lAlPos) + (this.Ns-lAirPos-1)*this.Ns;
			this.mat_al(2) = this.Ac_true(end-lAirPos, end-lAlPos) - this.al_fan_dis/this.C_al_hs;
			this.mat_alpos(2) = (this.Ns-lAirPos) + (this.Ns-lAlPos-1)*this.Ns;

			this.mat_al(3) = this.Ac_true(end-lAlPos, end-lAlPos) + this.al_fan_dis/this.C_al_hs;
			this.mat_alpos(3) = (this.Ns-lAlPos) + (this.Ns-lAlPos-1)*this.Ns;
		end

        % To define the data of the system, we implement 
        % the abstract functions from HybridSystem.m

		function this = init(this)
		
			this.Fa = this.F_target*ones(1,this.Nc);
			this.pid_integ = 0;
		end

        function xdot = flowMap(this, x, t, j)

			% Extract the state components.
			T = x(1:this.Ns);
			Tc = this.Cc*T;
			wl = this.wl_target;
			fan = x(end-this.CFp);
			I = x(1+this.Ns : this.Ns+this.Nc);

			%if abs(t-3.25)<1e-3 || abs(t-3.15)<1e-3
			%	aa=1;
			%end

			tt = I>this.aw_up;
			I = I + tt.*(this.aw_up-I);
			e_pid = this.T_pid - this.Cc*T;
			%this is to avoid strange errors
			e_pid = e_pid + (e_pid<-5).*(-5-e_pid );

			%banding
			e_pid = e_pid .* (~(e_pid>-1)&(e_pid<1));

			res = e_pid*this.kp + I*this.ki;
			%res = zeros(this.Nc,1);
			%fo = this.F_target + this.f_pid(Tc, this.T_pid, this.F_target);
			fo = this.F_target + min(res,0);
			fo = max(fo, 0.4);
			%if abs(this.Fa(end,:) - fo) >0.3
			%	aa = 1;
			%end
			this.Fa = [this.Fa; fo'];
			vdd_alpha = 0.3095; %0.2995
			vdd_offset = 0.07;
			Vo = vdd_offset + fo*vdd_alpha;% this.FV_table(sum(max(fo) > this.FV_table(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue
			%Vo = 1.0;

			P = this.power_compute(fo,Vo,Tc,wl,1,1);
			%this.pstore = this.pstore + sum(P);

			Tdot = this.Ac_true*T + this.Bc_true * [P; this.temp_amb*this.case_fan_nom_speed];
			
			%saturate Integral;
			e_pid = e_pid.*((~tt)|(e_pid<0));

			xdot = [Tdot; e_pid; 0; 0; 1/this.T_ctrl];
        end

        function xplus = jumpMap(this, x)

            % Extract the state components.
			T = this.Cc*x(1:this.Ns);
			%wl = this.wl_target;
			%pp = this.Nc+this.TCp-1;
			%twa = x(end-pp : end-this.TCp);
			%ctrla = x(end-this.qp);
			%fan = x(end-this.FANp);
			T_case = x(this.Ns-this.mb_pos);
			T_hs = max(T);
			prev_fan_hs = x(end-this.HFp);
			prev_fan_case = x(end-this.CFp); 

			%P = this.power_compute(this.Fa,this.Va,T,wl,1,1);
			%tt = T > this.T_case_fan;
			%dT = (T - this.T_case_fan) + 1;
			%tp = (sum(P)) > this.p_budget;

			%this setup is because control is applied at the next time step
			%fo = this.F_target;            
			%V = this.FV_table(sum(max(fo) > this.FV_table(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue

			%thermal capping
			%this.Fa = fo - twa*this.F_discretization_step.*dT;
			%this.Fa = this.Fa + (this.Fa>fo).*(fo-this.Fa);
			%this.Fa = fo;
			%this.Va = this.FV_table(sum(max(this.Fa) > this.FV_table(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue

			%diff = this.T_pid+5 - this.T_case_fan;
			%fan_case = 1 + sum(0.1*tt .* dT / diff*4*2);

			fan_hs = 1 + ((T_hs - this.tT_core) / this.sT_core)/ 2.2;
			fan_hs = floor(fan_hs*2) /2;
			fan_hs = max(fan_hs,1);
			fan_hs = min(prev_fan_hs+1, fan_hs);
			fan_case = 1 + ((T_case - this.tT_case) / this.sT_case);
			fan_case = floor(fan_case*2) /2;
			fan_case = max(fan_case,1);
			fan_case = min(prev_fan_case+1, fan_case);
			
			hsp = fan_hs * 1;
			al_fan_boost = (this.al_fan_dis/this.C_al_hs*hsp);

			if ~this.disable_fan_boost
			this.Ac_true(end) = this.mat_case - (this.case_fan_dis*fan_case) - (this.al_fan_dis/this.C_al_hs*hsp);
			this.Bc_true(end) = fan_case*this.case_fan_dis/this.case_fan_nom_speed;

			this.Ac_true(this.mat_alpos(1)) = this.mat_al(1) + al_fan_boost;
			this.Ac_true(this.mat_alpos(2)) = this.mat_al(2) + al_fan_boost;
			this.Ac_true(this.mat_alpos(3)) = this.mat_al(3) - al_fan_boost;
			end


			% Extract the state components.
			xplus = [x(1:this.Ns+this.Nc); fan_hs; fan_case; 0];            
        end
        
        function inC = flowSetIndicator(this, x)
            % Extract the state components.
			time = x(end);

            % Set 'inC' to 1 if 'x' is in the flow set and to 0 otherwise.
            inC = 1;
        end

        function inD = jumpSetIndicator(this, x)
            % Extract the state components.
            time = x(end);

			T = this.Cc*x(1:this.Ns);
			%tt = T > this.T_case_fan;
			T_case = x(this.Ns-this.mb_pos);
			T_hs = max(T);
			fan_hs = x(end-this.HFp);
			fan_case = x(end-this.CFp); 

			th_hs = 1 + (T_hs - this.tT_core) / this.sT_core / 2.2;
			th_hs = floor(th_hs*2) /2;
			th_hs = max(th_hs,1);
			th_case = 1 + (T_case - this.tT_case) / this.sT_case;
			th_case = floor(th_case*2) /2;
			th_case = max(th_case,1);

			time_cond = ((time-floor(time)) < 1e-2) & time > 1;
			if time >1
				aa=1;
			end

			
            % Set 'inD' to 1 if 'x' is in the jump set and to 0 otherwise.
            inD = time_cond && ( (th_hs~=fan_hs) || (th_case~=fan_case) );
		end


		function [upid] = f_pid(obj, T, pid_target, F)
			
			kp_l = obj.kp;
			ki_l = obj.ki;

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
			obj.pid_integ = obj.pid_integ + ki_l*e;
			eI = obj.pid_integ>obj.aw_up;
			obj.pid_integ = eI*obj.aw_up + (~eI).*obj.pid_integ;
			eI = obj.pid_integ<(F*obj.aw_down_c);
			obj.pid_integ = eI.*(F*obj.aw_down_c) + (~eI).*obj.pid_integ;

			% PID Proporitonal (& Output)
			upid = kp_l*e + obj.pid_integ;

			% PID Output Saturation
			eU = upid>obj.sat_up;
			upid = eU*obj.sat_up + (~eU).*upid;
			%TODO
			eU = upid < (-(F-ones(obj.Nc,1).*0.4));
			upid = eU.*(-(F-ones(obj.Nc,1).*0.4)) + (~eU).*upid;
		end
    end
end