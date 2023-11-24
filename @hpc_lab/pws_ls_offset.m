function [lut, F, T] = pws_ls_offset(obj, Fslot, Tslot, add_temp )
%PWS_LS_OFFSET Summary of this function goes here
%   Detailed explanation goes here
	tmin = 20+273.15;%obj.temp_amb;
	tmax = 90+273.15;%obj.core_crit_temp+add_temp;
	
	P0 = obj.ps_compute(obj.V_min, tmin,1,0);
	P1 = obj.ps_compute(obj.V_max, tmax,1,0);
	
	% T:
	slot = (P1-P0) / (Tslot);
	v = [P0:slot:P1];
	
	%TODO: this depends on the model I should have something like "inverse..."
	%vm = (obj.V_min + obj.V_max) / 2;
	%slot = (obj.V_max - obj.V_min) / Tslot
	%vp = [obj.V_min:slot:obj.V_max]
	%T1 = ( log(v/(obj.V_min*obj.leak_vdd_k + obj.leak_process_k)) ...
	%	- obj.leak_exp_vdd_k*obj.V_min - obj.leak_exp_k ) / obj.leak_exp_t_k
	T2 = ( log(v/(obj.V_max*obj.leak_vdd_k + obj.leak_process_k)) ...
		- obj.leak_exp_vdd_k*obj.V_max - obj.leak_exp_k ) / obj.leak_exp_t_k;
	%T3 = ( log(v/(vm*obj.leak_vdd_k + obj.leak_process_k)) ...
	%	- obj.leak_exp_vdd_k*vm - obj.leak_exp_k ) / obj.leak_exp_t_k
	%T4 = ( log(v./(vp*obj.leak_vdd_k + obj.leak_process_k)) ...
	%	- obj.leak_exp_vdd_k*vp - obj.leak_exp_k ) / obj.leak_exp_t_k
	
	%Tl = (T1+T2+T3+T4)/ 4
	
	%TODO this 273.15 (parametrize, and also fix it in formulas/math)
	T = T2 + 273.15;
	T = [tmin T(T>tmin)];

	%This is actually worse! it is better to make the intevarl even. Problably
	%	because in the middle the difference is greater. So I should make
	%	the non-linear interval in the inverse way, with more at the center
	%	or, even better, depending on the abs(diff)
	slot = (tmax-tmin) / (Tslot);
	T = [tmin:slot:tmax];
	
	%F
	slot = (obj.F_max-obj.F_min) / (Fslot);
	F = [obj.F_min:slot:obj.F_max];
	
	%
	%TODO: convert to formula
	V = obj.FV_table((sum(F>obj.FV_table(:,3))+1 ),1 );
	%TODO parametrize
	alp = 0.2995;
	lut = zeros(length(T), length(F));
	[h0, h1, h2] = obj.pws_ls_approx();
	for i=1:length(T)
		for j=1:length(F)
			ps = obj.ps_compute(V(j), T(i),1,0);
			%TODO: evaluate WHY is the one below, and not above!
			linps = h1*(T(i)-273.15) + h2*F(j)*V(j)^2 + h0;
			%linps = h1*(T(i)-273.15) + h2*F(j)^(2*alp+1) + h0;
			lut(i,j) = ps - linps;
		end
	end
	
end

