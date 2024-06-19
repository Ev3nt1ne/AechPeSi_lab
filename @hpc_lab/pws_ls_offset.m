function [lut, Vs, Ts, M_var] = pws_ls_offset(obj, ctrl, Vslot, Tslot, show )
%PWS_LS_OFFSET Summary of this function goes here
%   Detailed explanation goes here

	if ((nargin < 5) || isempty(show))
		show = 0;
    end

    add_temp = 0;

	tmin = obj.t_outside;
	tmax = obj.core_crit_temp + add_temp;
	
    %Initially I need to craft a value table/surface
    tsloti = round(tmax-tmin);
    %vsloti = 30;

    % T:
	slot = (tmax-tmin) / (tsloti-1);
	T = [tmin:slot:tmax];
    %F
	%slot = (obj.F_max-obj.F_min) / (vsloti-1);
	%F = [obj.F_min:slot:obj.F_max];
	%V = obj.FV_table((sum(F>obj.FV_table(:,3))+1 ),1 );
    F = obj.FV_table(:,3);
    V = obj.FV_table(:,1);

    %T = T(1:end-1);
    %F = F(1:end-1);
    %V = V(1:end-1);

    h0 = ctrl.psac(1);
    h1 = ctrl.psac(2);
	h2 = ctrl.psac(3);

    lut = zeros(length(T), length(V));
	
    for i=1:length(T)
		for j=1:length(V)
			ps = obj.ps_compute(V(j), T(i),1,0);
			%TODO: evaluate WHY is the one below, and not above!
			linps = h1*(T(i)-273.15) + h2*F(j)*V(j)^2 + h0;
			%linps = h1*(T(i)-273.15) + h2*F(j)^(2*alp+1) + h0;
			lut(i,j) = ps - linps;
		end
    end

    % Min regions
    %{
    N = Vslot*Tslot;
	[table, cas] = create_regions(T, V, lut, N, show);

    vtb = zeros(size(table));
    for i=1:N
        vtb(table==i) = cas(i,3);
    end
    %}

    vtb = lut;

    %matrix of initial variance:
        % do not need it
    %repeat until I have the exact number of columns and rows:
    exc_ext = 1;
    if exc_ext
        row_dim = size(vtb, 1) -2;
        col_dim = size(vtb, 2) -2;
    else
        row_dim = size(vtb, 1);
        col_dim = size(vtb, 2);
    end
    og_row_dim = size(vtb, 1);
    og_col_dim = size(vtb, 2);
    if exc_ext
        cidx = [1 1:col_dim col_dim];
        ridx = [1 1:row_dim row_dim];
    else
        cidx = 1:col_dim;
        ridx = 1:row_dim;
    end

    sequential_merging = zeros(row_dim+col_dim-(Vslot+Tslot), 2);
    sequential_merging_i = zeros(row_dim+col_dim-(Vslot+Tslot), 2);
    sm_idx = 1;

    %test with random
    %vtb = rand(size(vtb))*20-10;

    while ((row_dim > Tslot) || (col_dim > Vslot))
        storedVar = (og_row_dim)*(og_col_dim)*10;
        mergeCR = [0 0];
    
        %B = [];
        if (col_dim > Vslot)
            ldx = ridx;
        for i=1:(col_dim-1)
            idx = cidx;
            idx(idx>i) = idx(idx>i)-1;
            %idx = 1:col_dim;
            %idx(idx==i) = [];
            %idx(idx==(i+1)) = [];
            A = 0;
            for j=1:(col_dim-1)
                for k=1:row_dim
                    A = A + ...
                        var(vtb(ldx==k,idx==j), 0, "all"); %/(og_row_dim*sum(idx==j));
                end
            end
    
            A = A/((col_dim-1)*row_dim);
    
            %B = [B A];
            
            if (A < storedVar)
                storedVar = A;
                mergeCR = [0 i];
            end        
        end
        end %if
        if (row_dim > Tslot)
            ldx = cidx;
        for i=1:(row_dim-1)
            idx = ridx;
            idx(idx>i) = idx(idx>i)-1;
            A = 0;
            for j=1:(row_dim-1)
                for k=1:col_dim
                    A = A + ...
                        var(vtb(idx==j,ldx==k)', 0, "all"); %/(og_col_dim*sum(idx==j));
                end
            end
            
            A = A/((row_dim-1)*col_dim);
    
            %B = [B A];
    
            if (A < storedVar)
                storedVar = A;
                mergeCR = [i 0];
            end        
        end
        end %if

        sequential_merging(sm_idx,:) = mergeCR;
        alk1 = find(ridx==mergeCR(1),1,"last");
        alk1(isempty(alk1)) = 0;
        alk2 = find(cidx==mergeCR(2),1,"last");
        alk2(isempty(alk2)) = 0;
        sequential_merging_i(sm_idx,:) = [alk1 alk2];
        sm_idx = sm_idx+1;
    
        col_dim = col_dim - (mergeCR(2)>0);
        row_dim = row_dim - (mergeCR(1)>0);
    
        cidx(cidx>mergeCR(2)) = cidx(cidx>mergeCR(2))-(mergeCR(2)>0);
        ridx(ridx>mergeCR(1)) = ridx(ridx>mergeCR(1))-(mergeCR(1)>0);

    end %while

    Vs = zeros(Vslot-1,1);
    Ts = zeros(Tslot-1,1);
    for i=1:(row_dim-1)
        lidx = find(ridx==i,1,"last");
        Ts(i) = T(lidx);
    end
    for i=1:(col_dim-1)
        lidx = find(cidx==i,1,"last");
        Vs(i) = V(lidx);
    end

    lut = zeros(Tslot, Vslot);
    M_var = size(lut);

    plridx = 1;
    for i=1:row_dim
        lridx = find(ridx==i,1,"last");
        plcidx = 1;
        for j=1:col_dim
            lcidx = find(cidx==j,1,"last");
            lut(i, j) = mean(vtb(plridx:lridx, plcidx:lcidx), "all");
            M_var(i,j) = var(vtb(plridx:lridx, plcidx:lcidx), 0, "all");
            plcidx = lcidx+1;
        end
        plridx = lridx+1;
    end

    if show
        figure();
        surf(T,V,vtb');

        hold on;

        Tpl = zeros(Tslot, 1);
        Tpl(1) = (T(1)+Ts(1))/2;
        Tpl(end) = (T(end)+Ts(end))/2;
        for i=2:row_dim-1
            Tpl(i) = (Ts(i-1)+Ts(i)) / 2; 
        end
        Vpl = zeros(Vslot, 1);
        Vpl(1) = (V(1)+Vs(1))/2;
        Vpl(end) = (V(end)+Vs(end))/2;
        for i=2:col_dim-1
            Vpl(i) = (Vs(i-1)+Vs(i)) / 2; 
        end

        surf(Tpl, Vpl, lut','FaceAlpha',0.6);

        for i=1:(row_dim-1)
            lidx = find(ridx==i,1,"last");
            plot3(T(lidx)*ones(size(V)), V, vtb(lidx, :), 'LineWidth', 3.5);
        end
   
        for i=1:(col_dim-1)
            lidx = find(cidx==i,1,"last");
            plot3(T, V(lidx)*ones(size(T)), vtb(:,lidx), 'LineWidth', 3.5);
        end

        %{
        figure();

        hold on;

        colorplot = [0.8941    0.1020    0.1098;
            0.2157    0.4941    0.7216;
            0.3020    0.6863    0.2902;
            0.5961    0.3059    0.6392;
            1.0000    0.4980         0;
            1.0000    1.0000    0.2000;
            0.6510    0.3373    0.1569;
            0.9686    0.5059    0.7490;
            0.6000    0.6000    0.6000;];
    
        plidx = 1;
        for i=1:row_dim
            lidx = find(ridx==i,1,"last");
            surf(T(plidx:lidx)', V, vtb(plidx:lidx, :)', 'FaceColor', colorplot(mod(i, size(colorplot,1)),:),'FaceAlpha',0.5, 'EdgeColor', 'none' );
            plidx = lidx;
        end
   
        plidx = 1;
        for i=1:col_dim
            lidx = find(cidx==i,1,"last");
            surf(T, V(plidx:lidx), vtb(:,plidx:lidx)', 'EdgeColor', colorplot(mod(i, size(colorplot,1)),:), 'FaceAlpha',0,'LineWidth', 3.5 );
            plidx = lidx;
        end
        %}
    end
        
end


    %% OLD considerations
	% T:

    %P0 = obj.ps_compute(obj.V_min, tmin,1,0);
	%P1 = obj.ps_compute(obj.V_max, tmax,1,0);
	
    %slot = (P1-P0) / (Tslot);
	%v = [P0:slot:P1];
	
	%TODO: this depends on the model I should have something like "inverse..."
	%vm = (obj.V_min + obj.V_max) / 2;
	%slot = (obj.V_max - obj.V_min) / Tslot
	%vp = [obj.V_min:slot:obj.V_max]
	%T1 = ( log(v/(obj.V_min*obj.leak_vdd_k + obj.leak_process_k)) ...
	%	- obj.leak_exp_vdd_k*obj.V_min - obj.leak_exp_k ) / obj.leak_exp_t_k
	%T2 = ( log(v/(obj.V_max*obj.leak_vdd_k + obj.leak_process_k)) ...
	%	- obj.leak_exp_vdd_k*obj.V_max - obj.leak_exp_k ) / obj.leak_exp_t_k;
	%T3 = ( log(v/(vm*obj.leak_vdd_k + obj.leak_process_k)) ...
	%	- obj.leak_exp_vdd_k*vm - obj.leak_exp_k ) / obj.leak_exp_t_k
	%T4 = ( log(v./(vp*obj.leak_vdd_k + obj.leak_process_k)) ...
	%	- obj.leak_exp_vdd_k*vp - obj.leak_exp_k ) / obj.leak_exp_t_k
	
	%Tl = (T1+T2+T3+T4)/ 4
	
	%TODO this 273.15 (parametrize, and also fix it in formulas/math)
	%T = T2 + 273.15;
	%T = [tmin T(T>tmin)];

	%This is actually worse! it is better to make the intevarl even. Problably
	%	because in the middle the difference is greater. So I should make
	%	the non-linear interval in the inverse way, with more at the center
	%	or, even better, depending on the abs(diff)


