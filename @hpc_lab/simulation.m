function asimcl = ... %[cpxplot, cpuplot, cpfplot, cpvplot, wlop] = ...
	simulation(obj, ctrl, chip, show)
%SIMULATION Summary of this function goes here
%   Detailed explanation goes here	

    lt1 = length(ctrl);
    lt2 = length(chip);
    
    %{ 
    % With cell inputs)
    if xor(lt1,lt2)
        error("[LAB] Simulation inputs have wrong dimensions")
    end
    if lt2
        N_chip = length(ichip);
        if (N_chip ~= length(ictrl))
            error("[LAB] Simulation inputs have wrong dimensions")
        end
        chip = ichip;
        ctrl = ictrl;
    else
        N_chip = 1;
        chip{1} = ichip;
        ctrl{1} = ictrl;
    end
    %}
    % With array inputs
    if (lt1 ~= lt2)
        error("[LAB] Simulation inputs have wrong dimensions")
    end
    N_chip = lt1;
    
    Nsim = zeros(N_chip,1);
    sys_mul = zeros(N_chip,1);
    target_mul = zeros(N_chip,1);
    ctrl_mul = zeros(N_chip,1);
    nip_mul = zeros(N_chip,1);

    asimcl = hpc_simulation.empty(N_chip,0);

    for i=1:N_chip
	    chip(i).anteSimCheckTM();
	    obj.anteSimCheckLab(chip(i));
	    chip(i).anteSimCheckPM();

        % INIT
        asimcl(i) = hpc_simulation(chip(i),chip(i).Ad_true, chip(i).Bd_true);
	    %obj.init_compute_model(, chip.Nc, chip.vd, chip); %called automatically

	    % To Optimize execution, istead of having several functions called and
	    %	several ifs statement, we use ctrl.Ts_ctrl as the main time step of
	    %	the simulation.
	    Nsim(i) = ceil(obj.tsim / ctrl(i).Ts_ctrl);
    
	    % System:
	    % initial offset:
	    %sys_ts_offset = 0;
	    % frequency
	    sys_mul(i) = round(ctrl(i).Ts_ctrl / chip(i).Ts);
	    % output completed
	    %ctrl_output_mul = round(obj.ctrl_Texec / obj.Ts);
    
	    % Inputs
	    % initial offset:
	    target_ts_offset = 0;
	    % frequency
	    ts_div = obj.Ts_target / ctrl(i).Ts_ctrl;
	    target_mul(i) = floor(ts_div);
	    ctrl_mul(i) = floor(1/ts_div) -1;
	    ctrl_mul(i) = (ctrl_mul(i)<0)*0 + (ctrl_mul(i)>=0)*ctrl_mul(i);
	    asimcl(i).target_counter = target_mul(i) + target_ts_offset;
	    % Manage non-integer part:
	    % TODO: This does not work well with y.xx when xx> 60 and y = [2,3,4,5] 
	    %	because it removes too much. instead with 1, and 6, 7... works well 
	    %	I should investigate the math further.
	    if target_mul(i)>ctrl_mul(i)
		    ntd = (ts_div-target_mul(i))/target_mul(i)/target_mul(i);
		    asimcl(i).nip_sign = -1;
	    else
		    ntd = (1/ts_div)-(ctrl_mul(i)+1);
		    asimcl(i).nip_sign = +1;
	    end
	    nip_mul(i) = ceil(1/ntd);
	    asimcl(i).nip_counter(i) = nip_mul(i) + target_ts_offset;
    
	    %TODO
	    x = obj.t_init{i};% + (rand(obj.Ns,1) - 0.5*ones(obj.Ns,1));
        asimcl(i).cpuplot = zeros(Nsim(i)+1,chip(i).Nc);
	    asimcl(i).cpxplot = zeros(Nsim(i)+1,chip(i).Ns);
        asimcl(i).cpfplot = zeros(Nsim(i)+1,chip(i).Nc);
	    asimcl(i).cpvplot = zeros(Nsim(i)+1,chip(i).vd);
        asimcl(i).cpxplot(1,:) = x;
	    asimcl(i).cpuplot(1,:) = NaN;
	    asimcl(i).cpfplot(1,:) = asimcl(i).F;
	    asimcl(i).cpvplot(1,:) = asimcl(i).V;
	    asimcl(i).pvt{chip.PVT_P} = asimcl(i).process;
	    asimcl(i).pvt{chip.PVT_V} = [];
	    asimcl(i).pvt{chip.PVT_T} = ctrl(i).C*x;
    
	    ctrl(i) = ctrl(i).init_fnc(obj, chip(i), Nsim(i));
    end

    % Nsim multiplicity check
    Nsim_max = max(Nsim);
    rp = 1;
    for i=1:N_chip
        rp = rp && (mod(Nsim_max,Nsim(i))==0);
    end
    if (~rp)
        %TODO: Improve this instead of throw error
        error("[LAB] simulation length multiciplity are not in sync")
    end

    % Create Nsim multiciplity
    Nsim_mul = ones(N_chip,1);
    for i=1:N_chip
        Nsim_mul(i) = fix(Nsim_max / Nsim(i));
    end

    % Start pool of parallel workers
    %pool = parpool;

	% LOOOP
	for s=1:Nsim_max

        %parfor i=1:N_chip
        for i=1:N_chip
            if (mod(s,Nsim_mul)==0)
                % Input step managing
                asimcl(i).update_counters(target_mul, nip_mul, ctrl_mul);

                %Compute model:
		        index = 1+(s-1)*sys_mul(i);
		        [asimcl(i).cpuplot(index+1:index+sys_mul(i),:), asimcl(i).cpxplot(index+1:index+sys_mul(i),:), asimcl(i).wl, asimcl(i).pwm, asimcl(i)] = ...
                    asimcl(i).compute_model(sys_mul(i), asimcl(i).cpxplot(index,:)', asimcl(i).V, asimcl(i).F, asimcl(i).process, obj.temp_amb, obj.wltrc{i}, chip(i));	
        
		        % Sim output managing
		        T = ctrl(i).C*asimcl(i).cpxplot(index,:)';
		        % Noise:
		        dim = min(length(T), chip(i).Nc);
		        T = T + [(chip(i).sensor_noise)*( (rand(dim,1) - 0.5)*2 * chip(i).sensor_noise_amplitude(chip(i).PVT_T) ); zeros(length(T)-dim,1)];	
        
		        %if (mod(s-1 + ctrl_ts_offset, ctrl_mul) == 0)
			        asimcl(i).pvt{chip(i).PVT_P} = asimcl(i).process;
			        asimcl(i).pvt{chip(i).PVT_V} = [];
			        asimcl(i).pvt{chip(i).PVT_T} = T;
                    f_ref = obj.frtrc{i}(min(asimcl(i).target_index, size(obj.frtrc,1)),:)';
			        pwbdg = obj.tot_pw_budget(min(asimcl(i).target_index, length(obj.tot_pw_budget)));
			        [asimcl(i).F, asimcl(i).V, ctrl(i)] = ...
                        ctrl(i).ctrl_fnc(f_ref, pwbdg, asimcl(i).pvt, asimcl(i).ctrl_pwm, asimcl(i).ctrl_wl);
		        %end
        
		        asimcl(i).cpfplot(s+1,:) = asimcl(i).F;
		        asimcl(i).cpvplot(s+1,:) = asimcl(i).V;
            end
        end
    end

    %delete the parallelize
    %delete(pool);

	cmp = [];
	if obj.compare_vs_baseline
		if obj.pmc_need_update
			warning("[HPC LAB] Warning! You are not comparing with the baseline! You should run base_ideal_unr() or checkante!");
        end
        for i=1:N_chip
		    cmp{i} = obj.perf_max_check;
        end
    else
        for i=1:N_chip
		    cmp{i} = size(obj.wltrc{i},3)-1;
        end
    end
    for i=1:N_chip
	    asimcl(i).wlop = asimcl(i).wl_index ./ cmp{i} * 100;
	    ctrl(i) = ctrl(i).cleanup_fnc();

	    if show
		    % Pause because it is bugged on Linux
		    pause(0.5);
		    obj.xutplot(chip(i),asimcl(i).cpxplot,asimcl(i).cpuplot);
		    pause(0.5);
		    obj.powerconstrplot(chip(i),asimcl(i).cpuplot);
		    pause(0.5);
		    obj.tempconstrplot(chip(i),asimcl(i).cpxplot);
		    pause(0.5);
		    obj.perfplot(chip(i),asimcl(i).cpfplot,asimcl(i).wl_index, cmp{i}, i);	
		    pause(0.5);
		    obj.fvplot(chip(i),asimcl(i).cpfplot,asimcl(i).cpvplot);
    
		    t1 = chip(i).Ts*[1:Nsim(i)*sys_mul(i)]';
		    t2 = chip(i).Ts*sys_mul(i)*[1:Nsim(i)]';
    
		    ctrl(i) = ctrl(i).plot_fnc(t1, t2, asimcl(i).cpxplot, asimcl(i).cpuplot, asimcl(i).cpfplot, asimcl(i).cpvplot, asimcl(i).wlop);
        end
    end

end

