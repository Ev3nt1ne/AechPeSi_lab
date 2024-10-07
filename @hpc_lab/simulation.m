function asimcl = ... %[cpxplot, cpuplot, cpfplot, cpvplot, wlop] = ...
	simulation(obj, ctrl, chip, CM, show)
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

    if (size(CM)~=[N_chip,N_chip])
        error("[LAB] Wrong Communication MAtrix (CM) dimensions");
    end
    
    Nsim = zeros(N_chip,1);
    sys_mul = zeros(N_chip,1);
    target_mul = zeros(N_chip,1);
    ctrl_mul = zeros(N_chip,1);
    nip_mul = zeros(N_chip,1);

    asimcl = hpc_simulation.empty(N_chip,0);
    comm = cell(N_chip,1);
    ctrl_comm = cell(N_chip,N_chip);

    for cid=1:N_chip
	    chip(cid).anteSimCheckTM();
	    obj.anteSimCheckLab(chip(cid));
	    chip(cid).anteSimCheckPM();

        % INIT
        asimcl(cid) = hpc_simulation(chip(cid),chip(cid).Ad_true, chip(cid).Bd_true);
	    %obj.init_compute_model(, chip.Nc, chip.vd, chip); %called automatically

	    % To Optimize execution, istead of having several functions called and
	    %	several ifs statement, we use ctrl.Ts_ctrl as the main time step of
	    %	the simulation.
	    Nsim(cid) = ceil(obj.tsim / ctrl(cid).Ts_ctrl);
    
	    % System:
	    % initial offset:
	    %sys_ts_offset = 0;
	    % frequency
	    sys_mul(cid) = round(ctrl(cid).Ts_ctrl / chip(cid).Ts);
	    % output completed
	    %ctrl_output_mul = round(obj.ctrl_Texec / obj.Ts);
    
	    % Inputs
	    % initial offset:
	    target_ts_offset = 0;
	    % frequency
	    ts_div = obj.Ts_target / ctrl(cid).Ts_ctrl;
	    target_mul(cid) = floor(ts_div);
	    ctrl_mul(cid) = floor(1/ts_div) -1;
	    ctrl_mul(cid) = (ctrl_mul(cid)<0)*0 + (ctrl_mul(cid)>=0)*ctrl_mul(cid);
	    asimcl(cid).target_counter = target_mul(cid) + target_ts_offset;
	    % Manage non-integer part:
	    % TODO: This does not work well with y.xx when xx> 60 and y = [2,3,4,5] 
	    %	because it removes too much. instead with 1, and 6, 7... works well 
	    %	I should investigate the math further.
	    if target_mul(cid)>ctrl_mul(cid)
		    ntd = (ts_div-target_mul(cid))/target_mul(cid)/target_mul(cid);
		    asimcl(cid).nip_sign = -1;
	    else
		    ntd = (1/ts_div)-(ctrl_mul(cid)+1);
		    asimcl(cid).nip_sign = +1;
	    end
	    nip_mul(cid) = ceil(1/ntd);
	    asimcl(cid).nip_counter = nip_mul(cid) + target_ts_offset;
    
	    %TODO
	    x = obj.t_init{cid};% + (rand(obj.Ns,1) - 0.5*ones(obj.Ns,1));
        asimcl(cid).cpuplot = zeros(Nsim(cid)+1,chip(cid).Nc);
	    asimcl(cid).cpxplot = zeros(Nsim(cid)+1,chip(cid).Ns);
        asimcl(cid).cpfplot = zeros(Nsim(cid)+1,chip(cid).Nc);
	    asimcl(cid).cpvplot = zeros(Nsim(cid)+1,chip(cid).vd);
        asimcl(cid).cpxplot(1,:) = x;
	    asimcl(cid).cpuplot(1,:) = NaN;
	    asimcl(cid).cpfplot(1,:) = asimcl(cid).F;
	    asimcl(cid).cpvplot(1,:) = asimcl(cid).V;
	    asimcl(cid).pvt{chip.PVT_P} = asimcl(cid).process;
	    asimcl(cid).pvt{chip.PVT_V} = [];
	    asimcl(cid).pvt{chip.PVT_T} = ctrl(cid).C*x;
    
	    [ctrl(cid), comm{cid}] = ctrl(cid).init_fnc(obj, chip(cid), cid, Nsim(cid));
    end

    % Nsim multiplicity check
    Nsim_max = max(Nsim);
    rp = 1;
    for cid=1:N_chip
        rp = rp && (mod(Nsim_max,Nsim(cid))==0);
    end
    if (~rp)
        %TODO: Improve this instead of throw error
        error("[LAB] simulation length multiciplity are not in sync")
    end

    % Create Nsim multiciplity
    Nsim_mul = ones(N_chip,1);
    for cid=1:N_chip
        Nsim_mul(cid) = fix(Nsim_max / Nsim(cid));
    end

    % Start pool of parallel workers
    %pool = parpool;

	% LOOOP
	for s=1:Nsim_max
        
        % Manage Communications
        %TODO: I did not consider queue. I'm just updating with the new
        %   infoù
        %TODO: I'm not "consuming"/resetting the communication if delivered
        for i=1:N_chip
            for j=1:N_chip
                if (CM(i,j)~=0)
                    ctrl_comm{j, i} = comm{i};
                else
                    ctrl_comm{j, i} = {};
                end
            end
        end

        %TODO I did not implement offset timers. Also Very important for
        %   communications desync and delays
        %parfor i=1:N_chip
        for cid=1:N_chip
            if (mod(s,Nsim_mul)==0)
                % Input step managing
                asimcl(cid).update_counters(target_mul(cid), nip_mul(cid), ctrl_mul(cid));

                %Compute model:
		        index = 1+(s-1)*sys_mul(cid);
		        [asimcl(cid).cpuplot(index+1:index+sys_mul(cid),:), asimcl(cid).cpxplot(index+1:index+sys_mul(cid),:), asimcl(cid).wl, asimcl(cid).pwm, asimcl(cid)] = ...
                    asimcl(cid).compute_model(sys_mul(cid), asimcl(cid).cpxplot(index,:)', asimcl(cid).V, asimcl(cid).F, asimcl(cid).process, obj.temp_amb, obj.wltrc{cid}, chip(cid));	
        
		        % Sim output managing
		        T = ctrl(cid).C*asimcl(cid).cpxplot(index,:)';
		        % Noise:
		        dim = min(length(T), chip(cid).Nc);
		        T = T + [(chip(cid).sensor_noise)*( (rand(dim,1) - 0.5)*2 * chip(cid).sensor_noise_amplitude(chip(cid).PVT_T) ); zeros(length(T)-dim,1)];	
        
		        %if (mod(s-1 + ctrl_ts_offset, ctrl_mul) == 0)
			        asimcl(cid).pvt{chip(cid).PVT_P} = asimcl(cid).process;
			        asimcl(cid).pvt{chip(cid).PVT_V} = [];
			        asimcl(cid).pvt{chip(cid).PVT_T} = T;
                    f_ref = obj.frtrc{cid}(min(asimcl(cid).target_index, size(obj.frtrc{cid},1)),:)';
			        chippwb = obj.chip_pw_budget{cid}(min(asimcl(cid).target_index, length(obj.chip_pw_budget{cid})));
                    totpwb = obj.toto_pw_budget(min(asimcl(cid).target_index, length(obj.toto_pw_budget)));
                    pwbdg = [totpwb, chippwb];
			        [asimcl(cid).F, asimcl(cid).V, comm{cid}, ctrl(cid)] = ...
                        ctrl(cid).ctrl_fnc(f_ref, pwbdg, asimcl(cid).pvt, asimcl(cid).ctrl_pwm, asimcl(cid).ctrl_wl, cid, ctrl_comm(cid,:));
		            % Consumed the information
                    %ctrl_comm(i,:) = {};
                %end

		        asimcl(cid).cpfplot(s+1,:) = asimcl(cid).F;
		        asimcl(cid).cpvplot(s+1,:) = asimcl(cid).V;
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
        for cid=1:N_chip
		    cmp{cid} = obj.perf_max_check;
        end
    else
        for cid=1:N_chip
		    cmp{cid} = size(obj.wltrc{cid},3)-1;
        end
    end
    for cid=1:N_chip
	    asimcl(cid).wlop = asimcl(cid).wl_index ./ cmp{cid} * 100;
	    ctrl(cid) = ctrl(cid).cleanup_fnc();

	    if show
		    % Pause because it is bugged on Linux
		    pause(0.5);
		    obj.xutplot(chip(cid),asimcl(cid).cpxplot,asimcl(cid).cpuplot);
		    pause(0.5);
		    obj.powerconstrplot(chip(cid),cid,asimcl(cid).cpuplot);
		    pause(0.5);
		    obj.tempconstrplot(chip(cid),asimcl(cid).cpxplot);
		    pause(0.5);
		    obj.perfplot(chip(cid),asimcl(cid).cpfplot,asimcl(cid).wl_index, cmp{cid}, cid);	
		    pause(0.5);
		    obj.fvplot(chip(cid),asimcl(cid).cpfplot,asimcl(cid).cpvplot);
    
		    t1 = chip(cid).Ts*[1:Nsim(cid)*sys_mul(cid)]';
		    t2 = chip(cid).Ts*sys_mul(cid)*[1:Nsim(cid)]';
    
		    ctrl(cid) = ctrl(cid).plot_fnc(t1, t2, asimcl(cid).cpxplot, asimcl(cid).cpuplot, asimcl(cid).cpfplot, asimcl(cid).cpvplot, asimcl(cid).wlop);
        end
    end

end

