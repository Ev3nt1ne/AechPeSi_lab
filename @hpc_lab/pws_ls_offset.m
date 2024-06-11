function [lut, F, T] = pws_ls_offset(obj, ctrl, Vslot, Tslot, show )
%PWS_LS_OFFSET Summary of this function goes here
%   Detailed explanation goes here

	if ((nargin < 5) || isempty(show))
		show = 0;
    end

    add_temp = 5;

	tmin = obj.t_outside;
	tmax = obj.core_crit_temp + add_temp;
	
    %Initially I need to craft a value table/surface
    tsloti = 100;
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

    if show
        figure()
	    surf(T,V, lut');
    end

    % Min regions
    N = Vslot*Tslot;
	[table, cas] = create_regions(T, V, lut, N, show);

    vtb = zeros(size(table));
    for i=1:N
        vtb(table==i) = cas(i,3);
    end

    vtb

    %matrix of initial variance:
        % do not need it
    %repeat until I have the exact number of columns and rows:
    row_dim = size(vtb, 1);
    col_dim = size(vtb, 2);
    og_row_dim = row_dim;
    og_col_dim = col_dim;            
    cidx = 1:col_dim;
    ridx = 1:row_dim;

    while ((row_dim > Tslot) || (col_dim > Vslot))
        storedVar = 1000;
        mergeCR = [0 0];
    
        B = [];
        if (col_dim > Vslot)
        for i=1:(col_dim-1)
            idx = cidx;
            idx(idx>i) = idx(idx>i)-1;
            %idx = 1:col_dim;
            %idx(idx==i) = [];
            %idx(idx==(i+1)) = [];
            A = 0;
            for j=1:(col_dim-1)
                A = A + ...
                    var(vtb(:,idx==j), 0, "all")/(og_row_dim*sum(idx==j));
            end
    
            A = A/(col_dim-1);
    
            B = [B A];
            
            if (A < storedVar)
                storedVar = A;
                mergeCR = [0 i];
            end        
        end
        end %if
        if (row_dim > Tslot)
        for i=1:(row_dim-1)
            idx = ridx;
            idx(idx>i) = idx(idx>i)-1;
            A = 0;
            for j=1:(row_dim-1)
                A = A + ...
                    var(vtb(idx==j,:)', 0, "all")/(og_col_dim*sum(idx==j));
            end
            
            A = A/(row_dim-1);
    
            B = [B A];
    
            if (A < storedVar)
                storedVar = A;
                mergeCR = [i 0];
            end        
        end
        end %if
    
        col_dim = col_dim - (mergeCR(2)>0);
        row_dim = row_dim - (mergeCR(1)>0);
    
        cidx(cidx>mergeCR(2)) = cidx(cidx>mergeCR(2))-(mergeCR(2)>0);
        ridx(ridx>mergeCR(1)) = ridx(ridx>mergeCR(1))-(mergeCR(1)>0);

    end %while

    hold on;

    for i=2:row_dim
        lltt = min(find(ridx==i));
        plot3(T(lltt)*ones(size(V)), V, lut(lltt, :), 'LineWidth', 2.5)
    end

    for i=2:col_dim
        lltt = min(find(cidx==i));
        plot3(T, V(lltt)*ones(size(T)), lut(:,lltt), 'LineWidth', 2.5)
    end

    aas = 1;
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


