classdef hpc_simulation < handle
    %HPC_SIMULATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        F;
        V;
        process;
        uop;
        xop;
        fop;
        vop;
        wlop;

        pwm;
        ctrl_pwm;
        wl;
        ctrl_wl;
        pvt;

        target_counter;
        target_index;
        nip_counter;
        nip_sign;
    end
    properties(SetAccess=protected, GetAccess=public)
        wl_index;
	 	qt_storage;

		F_cng_times;
		V_cng_times;
		F_cng_us;
		V_cng_us;
		F_cng_error;
		V_cng_error;

		V_T;
		F_T;

	 	V_s;
	 	F_s;

		delay_F_index;
		delay_V_index;

		delay_F_div;
		delay_V_div;

	 	A_s;
        B_s;
    end
    
    methods
        function obj = hpc_simulation(chip, A, B)
            %HPC_SIMULATION Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj.init_compute_model(chip, A, B);
        end
    end
        %% Simulations
	methods
		[] = sim_tm_autonomous(obj, ts, Nsample, exp_gamma)
		function [obj] = init_compute_model(obj, chip, A, B)
			obj.qt_storage = chip.quantum_instr*ones(chip.Nc, 1);
			obj.wl_index = ones(chip.Nc, 1);
			obj.F_cng_times = zeros(chip.Nc,1);
			obj.V_cng_times = zeros(chip.vd,1);
			obj.F_cng_us = zeros(chip.Nc,1);
			obj.V_cng_us = zeros(chip.vd,1);
			obj.F_cng_error = zeros(chip.Nc,1);
			obj.V_cng_error = zeros(chip.vd,1);

			obj.F_T = chip.F_min*ones(chip.Nc,1);
			obj.V_T = chip.V_min*ones(chip.vd,1);
			obj.F_s = chip.F_min*ones(chip.Nc,1);
			obj.V_s = chip.V_min*ones(chip.vd,1);
			obj.delay_F_index = zeros(chip.Nc,1);
			obj.delay_V_index = zeros(chip.vd,1);
			obj.delay_F_div = ceil(chip.delay_F_mean / chip.Ts);
			obj.delay_V_div = ceil(chip.delay_V_mean / chip.Ts);
			obj.A_s = A;
			obj.B_s = B;
			% cannot put it here, if I want it to be constant among several
			% runs
			%obj.core_pw_noise_char = ones(obj.Nc,1);
			%if (obj.model_variation)
			%obj = obj.create_core_pw_noise();
			%obj = obj.create_thermal_model_noise();
			%end
            obj.F = chip.F_min*ones(chip.Nc,1);
	        obj.V = chip.V_min*ones(chip.vd,1);
	        obj.process = ones(chip.Nc,1);
	        obj.pwm = 1;
	        obj.wl  = [ones(chip.Nc,1) zeros(chip.Nc, chip.ipl -1)];

            obj.target_counter = 0;
            obj.target_index = 0;
            obj.nip_counter = 0;
		end
        function [obj] = update_counters(obj, target_mul, nip_mul, ctrl_mul)
            obj.target_counter = (obj.target_counter-1)*(obj.target_counter>0) + (obj.target_counter<=0)*(target_mul-1);
		    obj.nip_counter = (obj.nip_counter-1)*(obj.nip_counter>0) + (obj.nip_counter<=0)*(nip_mul-1);
		    obj.target_index = obj.target_index + (obj.target_counter==target_mul-1) + obj.nip_sign*(obj.nip_counter==nip_mul-1) + ctrl_mul;

            obj.ctrl_pwm = obj.pwm;
            obj.ctrl_wl = obj.wl;
        end        
        function [uplot, xplot, d_is, pw_ms, obj] = compute_model(obj, N, x, V, F, d_p, temp_amb, wltrc, chip)
			
			d_is = zeros(1,chip.ipl);
			%delay_F_index = zeros(obj.Nc,1);
			%delay_V_index = zeros(obj.vd,1);		
			pw_ms = 0;			
			
			uplot = zeros(N, length(F));
			xplot = zeros(N, length(x));

			%F
			ttd = (obj.delay_F_index > 0)|(obj.F_T ~= obj.F_s);
			cng = (~ttd).*(obj.F_T ~= F);

			obj.F_cng_times = obj.F_cng_times + cng;

			obj.delay_F_index = obj.delay_F_index.*ttd + cng.*obj.delay_F_div;
			obj.F_T = (~cng).*obj.F_T + cng.*F;
			obj.F_cng_error = obj.F_cng_error + ttd.*(obj.F_T ~= F);

			%V
			ttd = (obj.delay_V_index > 0)|(obj.V_T ~= obj.V_s);
			cng = (~ttd).*(obj.V_T ~= V);

			obj.V_cng_times = obj.V_cng_times + cng;

			obj.delay_V_index = obj.delay_V_index.*ttd + cng.*obj.delay_V_div;
			obj.V_T = (~cng).*obj.V_T + cng.*V;
			obj.V_cng_error = obj.V_cng_error + ttd.*(obj.V_T ~= V);
			
			for sim=1:N
				%index = (ix0-1)*N + sim;

				% Delay Application
				tv = obj.delay_V_index == 0;
				obj.V_s = obj.V_s.*(~tv) + obj.V_T.*tv;
				tf = obj.delay_F_index == 0;
				obj.F_s = obj.F_s.*(~tf) + obj.F_T.*tf;

				obj.F_cng_us = obj.F_cng_us + ~tf;
				obj.V_cng_us = obj.V_cng_us + ~tv;

				obj.delay_F_index = obj.delay_F_index - (obj.F_T~=obj.F_s);
				obj.delay_V_index = obj.delay_V_index - (obj.V_T~=obj.V_s);
				%

				%noise 

				instr = obj.F_s * 1e9 * (chip.Ts);
				lwl = zeros(chip.Nc, chip.ipl);
				wld = zeros(chip.Nc, 1);
				while (sum(instr) > 0)
					vidx = max(mod(obj.wl_index, size(wltrc,3)),1);
					[m,n,l] = size(wltrc);
					idx = sub2ind([m,n,l],repelem(1:m,1,n),repmat(1:n,1,m),repelem(vidx(:).',1,n));
					wlp = reshape(wltrc(idx),[],m).';
					
					[pwl, instr] = chip.quanta2wl(wlp, instr, (chip.F_min./obj.F_s), obj.qt_storage);					

					% add to accumulator
					wld = wld + pwl;
					lwl = lwl + pwl .* wlp;

					% compute new values
					obj.qt_storage = obj.qt_storage - pwl;
					ttwl = (obj.qt_storage<=0);
					obj.qt_storage = obj.qt_storage + chip.quantum_instr .* ttwl;
					obj.wl_index = obj.wl_index + ttwl;
				end

				lwl = lwl ./ wld;

				d_is = d_is + lwl;

				pu_s = chip.power_compute(obj.F_s,chip.VDom*obj.V_s,chip.C(1:chip.Nc,:)*x,lwl,d_p);
				pw_ms = pw_ms + sum(pu_s);
				x = obj.A_s*x + obj.B_s*[pu_s;temp_amb*1000];		

				uplot(sim,:) = pu_s;
				xplot(sim,:) = x;
			end	%ssim
			d_is = d_is / N;
			pw_ms = pw_ms / N;
		end
	end
        
end

