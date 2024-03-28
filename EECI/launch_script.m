
%%
addpath ../

ms = TPM();
% Core number
ms.Nc = 9;

% Rows and Columns
ms.Nv = 3;
ms.Nh = 3;

% Number of Voltage Domains
ms.vd = 3;

% How the cores are distributed per domain
%hpc.VDom = ... ;
% Alternatively
ms.default_VDom_config();

% Default Floorplan and Parameter deviation vector
ms.default_floorplan_config();
ms.create_model_deviation();
ms.model_init();
ms.create_core_pw_noise();


%%

ms.init();
config = HybridSolverConfig();

ms.disable_fan_boost = 0;

x0 = [(45+273.15)*ones(ms.Ns-1,1); (35+273.15); zeros(ms.Nc,1); 1; 1; 0];
tspan = [0, 4];
jspan = [0, tspan(2)/ms.T_ctrl+2];
%config = HybridSolverConfig('refine', 32); % Improves plot smoothness for demo.
config.priority('jump');
sol = ms.solve(x0, tspan, jspan, config);
%plotFlows(sol); % Display solution

%%
% ********* PLOTTING
figure();

plot_color = [    0.1490    0.5490    0.8660;
    0.9600    0.4660    0.1600;
    1.0000    0.9090    0.3920;
    0.7520    0.3600    0.9840;
    0.2860    0.8580    0.2500;
    0.4230    0.9560    1.0000;
    0.9490    0.4030    0.7720;
    0.9960    0.7520    0.2980;
    0.4900    0.6620    1.0000;
    1.0000    0.4780    0.4540;
    0.1210    0.8110    0.7450;
    0.8620    0.6000    0.4230];

lsname = ["Heat-sink", "PCB", "Motherboard", "Air"];

hybrid_arc1 = sol.select(1:ms.Ns);                   % Pick component.
hybrid_arc1 = hybrid_arc1.transform(@(x) x-273.15);
%hybrid_arc = hybrid_arc.restrictT([1.5, 12]); % Truncate to t-values between 4.5 and 7.
%hybrid_arc = hybrid_arc.restrictJ([2, inf]);  % Truncate to j-values >= 2.

ee = length(sol.x0);
hybrid_arc2 = sol.select(ee-ms.HFp:ee-ms.HFp+1);
hybrid_arc2 = hybrid_arc2.transform(@(x) x.*[1000;1000]);


%hybrid_arc1 = hybrid_arc1.restrictT([0, 5]);
%hybrid_arc2 = hybrid_arc2.restrictT([0, 5]);

% Plot hybrid arcs
%figure();
tsp = 3;
for sb=1:2
	ax(sb) = subplot(tsp,1,sb);
	%pa = hybrid_arc1.select(sb:2:ms.Nc*2);
	%plotFlows(pa);
	hpb = HybridPlotBuilder();
	hold on;
	if sb==1
	name = "Core ";
	else
	name = "heat-spreader ";
	end
	for i=1:ms.Nc
		hpb.color(plot_color(i,:)).legend(strcat(name, num2str(i))).plotFlows(hybrid_arc1.select((i-1)*2+sb)).flowLineWidth(2);
		%hpb.plotFlows(hybrid_arc1)
		%plotHybrid(hybrid_arc);
	end
	ylabel("Temperature [°C]");
	xlabel("Time [s]");
	set(gca,'FontSize',16)
end
ax(3) = subplot(tsp,1,3);
hpb = HybridPlotBuilder();
hold on;
for i=1:ms.add_states
	hpb.color(plot_color(i,:)).legend(lsname(i)).plotFlows(hybrid_arc1.select(ms.Ns-ms.add_states+i)).flowLineWidth(2);
	%hpb.plotFlows(hybrid_arc1)
	%plotHybrid(hybrid_arc);
end
ylabel("Temperature [°C]");
xlabel("Time [s]");
set(gca,'FontSize',16)

yyaxis right;
hpb.color(plot_color(ms.add_states+1,:)).legend("Heat-Sink Fan").plotFlows(hybrid_arc2.select(1)).flowLineWidth(2);
hold on;
hpb.color(plot_color(ms.add_states+2,:)).legend("Case Fan").plotFlows(hybrid_arc2.select(2)).flowLineWidth(2);
ylabel("Fan Speed [RPM]");
%hpb.plotFlows(sol.select(ee-2));

yyaxis left;
if tsp>3
	ax(4) = subplot(4,1,4);
	steps = (120) / (length(ms.Fa)-1);
	plot(0:steps:120,ms.Fa);
end
ylabel("Temperature [°C]");
xlabel("Time [s]");

linkaxes(ax, 'x');
set(gca,'FontSize',16)


%%
mean(ms.Fa,"all")
%mean(ms.Fa)

mean(hybrid_arc1.x(:,1:2:ms.Nc*2),"all")
%mean(hybrid_arc1.x(:,1:2:ms.Nc*2))

mean(hybrid_arc1.x(:,ms.Ns-3:ms.Ns))
%%

[this.Ac_true(this.mat_alpos(1)) this.Ac_true(this.mat_alpos(2)) this.Ac_true(this.mat_alpos(3))]
%%

figure();
plot(ms.Fa);
mean(ms.Fa)
mean(hybrid_arc1.x(:,1:2:2*ms.Nc))
sol.xf(ms.Ns-ms.al_pos) - 273.15
sol.xf(ms.Ns-ms.air_pos) - 273.15


figure();
ee = length(sol.x0);
hybrid_arc2 = sol.select(ee-ms.HFp:ee-ms.HFp+1);                   % Pick component.
%hybrid_arc = hybrid_arc.transform(@(x) x(1:ms.Ns));
%hybrid_arc = hybrid_arc.restrictT([1.5, 12]); % Truncate to t-values between 4.5 and 7.
%hybrid_arc = hybrid_arc.restrictJ([2, inf]);  % Truncate to j-values >= 2.

% Plot hybrid arcs
hpb = HybridPlotBuilder();
%hpb.color('black').legend('Original').plotFlows(sol.select(1));
%hold on
hpb.plotFlows(hybrid_arc2)
%plotHybrid(hybrid_arc);

figure();
pos = ms.Ns+1;
hybrid_arc3 = sol.select(pos:pos+ms.Nc-1);                   % Pick component.
%hybrid_arc3 = hybrid_arc3.transform(@(x) x);
%hybrid_arc = hybrid_arc.restrictT([1.5, 12]); % Truncate to t-values between 4.5 and 7.
%hybrid_arc = hybrid_arc.restrictJ([2, inf]);  % Truncate to j-values >= 2.

% Plot hybrid arcs
hpb = HybridPlotBuilder();
%hpb.color('black').legend('Original').plotFlows(sol.select(1));
%hold on
hpb.plotFlows(hybrid_arc3)
%plotHybrid(hybrid_arc);

%%
figure();
hybrid_arc10 = sol.select(ee-1);                   % Pick component.
%hybrid_arc10 = hybrid_arc10.transform(@(x) x-273.15);
%hybrid_arc = hybrid_arc.restrictT([1.5, 12]); % Truncate to t-values between 4.5 and 7.
%hybrid_arc = hybrid_arc.restrictJ([2, inf]);  % Truncate to j-values >= 2.

% Plot hybrid arcs
clf
hpb = HybridPlotBuilder();
%hpb.color('black').legend('Original').plotFlows(sol.select(1));
%hold on
hpb.plotFlows(hybrid_arc10)
%plotHybrid(hybrid_arc);


%%
% Extract the state components.
	T = this.Cc*x(1:this.Ns);
	%wl = this.wl_target;
	%pp = this.Nc+this.TCp-1;
	%twa = x(end-pp : end-this.TCp);
	%ctrla = x(end-this.qp);
	%fan = x(end-this.FANp);

	%P = this.power_compute(this.Fa,this.Va,T,wl,1,1);
	tt = T > this.T_case_fan;
	dT = (T - this.T_case_fan) + 1;
	%tp = (sum(P)) > this.p_budget;

	%this setup is because control is applied at the next time step
	%fo = this.F_target;            
	%V = this.FV_table(sum(max(fo) > this.FV_table(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue

	%thermal capping
	%this.Fa = fo - twa*this.F_discretization_step.*dT;
	%this.Fa = this.Fa + (this.Fa>fo).*(fo-this.Fa);
	%this.Fa = fo;
	%this.Va = this.FV_table(sum(max(this.Fa) > this.FV_table(:,3)+1e-6)+1,1); %+1e-6 to fix matlab issue
	diff = this.T_pid+5 - this.T_case_fan;
	fan_case = 1 + sum(0.1*tt .* dT / diff*4*2);

	this.Ac_true(end) = this.mat_value - (this.case_fan_dis*fan_case);
	this.Bc_true(end) = fan_case*this.case_fan_dis/this.case_fan_nom_speed;