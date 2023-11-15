function [k0, k1, k2] = pws_ls_approx(obj)

	Fint = [obj.F_min obj.F_Max];
	Tint = [20 90];

	%Fint = [2.8 obj.F_Max];
	%Tint = [60 90];

	a = Tint(1);
	b = Tint(2);
	c = Fint(1);
	d = Fint(2);

	Ke = obj.leak_exp_k;
	KT = obj.leak_exp_t_k;
	KV = obj.leak_exp_vdd_k;
	Ks = obj.leak_process_k;
	Kw = obj.leak_vdd_k;

	V0 = 0.9;
	%V0 = 1.2;
	KV1 = exp(KV*V0)*(1-KV*V0);
	KV2 = KV / (1 - KV*V0);

	alp = 0.2995; %0.15;

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