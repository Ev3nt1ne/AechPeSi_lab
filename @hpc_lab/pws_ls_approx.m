function [k0, k1, k2] = pws_ls_approx(obj, I, T, Toff, C, alp, alp_I0, using_voltage)

	if (nargin < 8) || isempty(using_voltage)
		using_voltage = 0;
    end

	a = T(1)+Toff;
	b = T(2)+Toff;
	c = I(1);
	d = I(2);

	KT = obj.leak_exp_t_k;
    Ke = obj.leak_exp_k - Toff*KT;
	KV = obj.leak_exp_vdd_k;
	Ks = obj.leak_process_k;
	Kw = obj.leak_vdd_k;

	if using_voltage

		c1 = exp(Ke)*(exp(KT*b) - exp(KT*a)) / (KT*(d-c)*(b-a));
		c2 = exp(KV)/(KV^2);
		c3 = (exp(d) - exp(c))*(KV*Ks - Kw) + (d*exp(d) - c*exp(c))*KV*Kw;

		m0 = c1*c2*c3;

		c1 = 6*exp(Ke)*( KT*(b-a)*(exp(KT*b) + exp(KT*a)) - 2*(exp(KT*b) - exp(KT*a)) ) / (KT^2*(d-c)*(b-a)^3);
		c2 = exp(KV)/(KV^2);
		c3 = (exp(d) - exp(c))*(KV*Ks - Kw) + (d*exp(d) - c*exp(c))*KV*Kw;

		m1 = c1*c2*c3;

		if length(alp) < 2
			%Exponential

			c1 = exp(Ke)*(exp(KT*b) - exp(KT*a)) / (KT*KV^(5)*(b-a));
			
			Viarr = [c, d];
			V0 = C;
			p0 = V0^(1/alp) * (alp-1)/alp;
			Kp = V0^((1-alp)/alp)/alp;
			eta = (d+c)/2;
	
			Ndc = [0;0];
	
			for i=1:2
				Vi = Viarr(i);
				l1 = exp(KV*Vi);
				l2 = KV^4 * (Kw*Vi + Ks)*(Vi^2 * (p0 + Kp*Vi) - eta);
				l3 = -KV^3 * (Kp * Vi^2 * (4*Kw*Vi + 3*Ks) + p0*Vi*(3*Kw*Vi + 2*Ks));
				l4 = 2*KV^2*(3*Kp*Vi*(2*Kw*Vi+Ks) + p0*(3*Kw*Vi+Ks));
				l5 = -6*KV*(Kp*(4*Kw*Vi+Ks)+p0*Kw) + 24*Kp*Kw;
	
				Ndc(i) = l1*(l2+l3+l4+l5);
			end
	
			c2 = Ndc(2) - Ndc(1);
	
			Ddc = [0;0];
			for i=1:2
				Vi = Viarr(i);
	
				l1 = Vi^5 /5;
				l11 = 10*alp*alp_I0*Vi^(1/alp) / (5*alp +1);
				l12 = 5*alp*Vi^(2/alp) / (5*alp+2);
				l13 = alp_I0^2;
	
				l2 = -2*eta*Vi^3 / 3;
				l21 = 3*alp*Vi^(1/alp) / (3*alp + 1);
				l22 = alp_I0;
	
				Ddc(i) = l1*(l11+l12+l13) + l2*(l21+l22);			
			end
	
			c3 = eta^2*(d-c) + (Ddc(2) - Ddc(1));

		else
			%Polynomial
			c1 = exp(Ke)*(exp(KT*b) - exp(KT*a)) / (KT*KV^(6)*(b-a));

			Viarr = [c, d];
			eta = (d+c)/2;
	
			Ndc = [0;0];
	
			for i=1:2
				Vi = Viarr(i);
				l1 = exp(KV*Vi);
				l2 = KV^5 * Vi^2 * (Kw*Vi + Ks)*(alp_I0 + Vi*(alp(1) + alp(2)*Vi));
				l3 = -KV^4 * Vi * ( alp_I0*(3*Kw*Vi + 2*Ks) + Vi*( alp(1)*(4*Kw*Vi + 3*Ks) + alp(2)*Vi*(5*Kw*Vi+4*Ks) ) );
				l34 = -KV^4 / 2 *(c+d) * (KV*(Kw*Vi + Ks) - Kw);
				l4 = 2*KV^3*(alp_I0*(3*Kw*Vi+Ks) + Vi*(3*alp(1)*(2*Kw*Vi+Ks) + 2*alp(2)*Vi*(5*Kw*Vi+2*Ks) ) );
				l5 = -6*KV^2*( alp_I0*Kw + alp(1)*(4*Kw*Vi+Ks) + 2*alp(2)*Vi*(5*Kw*Vi+2*Ks) );
				l6 = 24*KV*(alp(1)*Kw + alp(2)*(5*Kw*Vi + Ks)) - 120*alp(2)*Kw;
	
				Ndc(i) = l1*(l2+l3+l34+l4+l5+l6);
			end

			c2 = Ndc(2) - Ndc(1);

			Ddc = [0;0];
			for i=1:2
				Vi = Viarr(i);
	
				l1 = alp_I0^2 / 5;
				l2 = (2*alp_I0*alp(2) + alp(1)^2) / 7;
				l3 = alp(1)*alp(2)/4 + alp(2)^2*Vi/9;
				l4 = alp_I0*alp(1)/3;
				
				l0 = alp_I0*Vi^3/3 + alp(1)*Vi^4/4 + alp(2)*Vi^5/5;
					
				Ddc(i) = -2*eta*l0 + Vi^5*(l1 + Vi*( Vi*(l2 + Vi*l3) + l4));
			end

			c3 = eta^2*(d-c) + (Ddc(2) - Ddc(1));

		end

		m2 = c1*c2/c3;

		k1 = m1;
		k2 = m2;
		k0 = m0 - (m1*(b+a) + m2*(d+c))/2;	

	else
		%% using Frequency
		% I = F
		% V0 = C
		% V0 = C

		%Fint = [obj.F_min obj.F_max];
		%Tint = [20 90];
	
		%Fint = [2.8 obj.F_Max];
		%Tint = [60 90];
	
		%V0 = 0.9;
		V0 = C;
		KV1 = exp(KV*V0)*(1-KV*V0);
		KV2 = KV / (1 - KV*V0);
	
		%alp = 0.2995; %0.15;
	
		c1 = exp(Ke)*KV1*(exp(KT*b) - exp(KT*a)) / (KT*(d-c)*(b-a));
		c2 = (d^(alp+1) - c^(alp+1))*(Ks*KV2 + Kw) / (alp+1);
		c3 = KV2*Kw*(d^(2*alp+1) - c^(2*alp+1)) / (2*alp+1);
		c4 = Ks*(d-c);
	
		m0 = c1 * (c2 + c3 + c4);
	
		c1 = 6*exp(Ke)*KV1*( KT*(b-a)*(exp(KT*b) + exp(KT*a)) - 2*(exp(KT*b) - exp(KT*a)) ) / (KT^2*(d-c)*(b-a)^3);
		c2 = (d^(alp+1) - c^(alp+1))*(Ks*KV2 + Kw) / (alp+1);
		c3 = KV2*Kw*(d^(2*alp+1) - c^(2*alp+1)) / (2*alp+1);
		c4 = Ks*(d-c);
	
		m1 = c1 * (c2 + c3 + c4);
	
		c1 = 2*exp(Ke)*KV1*(exp(KT*b) - exp(KT*a)) / (KT*(b-a));
		cd1 = (d-c)*(d+c)^2 - 2*(d^(2*alp+2) - c^(2*alp+2))*(d+c)*(alp+1)^(-1) + 4*(d^(4*alp+3) - c^(4*alp+3))*(4*alp+3)^(-1);
		c2 = (d^(alp+1) - c^(alp+1))*(Ks*(d^(alp+1)+c^(alp+1))-(d+c)*(Ks*KV2+Kw))/(alp+1);
		c3 = KV2*Kw*(d^(2*alp+1)-c^(2*alp+1))*((d^(2*alp+1)+c^(2*alp+1)) - (d+c)) / (2*alp+1);
		c4 = 2*(d^(3*alp+2)-c^(3*alp+2))*(Ks*KV2 + Kw) / (3*alp+2);
		c5 = Ks*(d^2 - c^2);
	
		m2 = c1/cd1 * (c2 + c3 + c4 - c5);
	
		k1 = m1;
		k2 = m2;
		k0 = m0 - (m1*(b+a) + m2*(d+c))/2;	

	end

end