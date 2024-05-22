


V = [0.5:0.05:1.2];
d_p = 1;
T = [20:1:100];

P = [hpc.static_pw_coeff hpc.exp_leak_coeff]/ 1000;

Power_static = V*P(1) + d_p*P(3);
maxT = 125+273.15;

res = [];
plance = [];
for i=1:length(T)
	res(:,i) = (Power_static .* exp(V*P(4) + T(i)*P(5) - P(6)))'*0.95 + 0.05;
	plane(:,i) = Power_static;
end

F = hpc.FV_table( [sum(V>hpc.FV_table(:,1))+1]',3);

ci = hpc.ceff_pw_coeff / 1000;
Ceff_low = [1 zeros(1,hpc.ipl-1)] * ci';
Ceff_high = [zeros(1,hpc.ipl-1) 1] * ci';
Power_dyn_low = Ceff_low .* F .* (V' .* V');
Power_dyn_high = Ceff_high .* F .* (V' .* V');

%%
fig = figure()

subplot(3, 3, [1:2, 4:5]);

surf(T, V, res, 'EdgeAlpha', 0.6)
hold;
h = gca;
surf([T(1) T(end)], V, [Power_static(:), Power_static(:)], 'FaceColor', 'm', ... %"#A2142F", 
	'FaceAlpha', 0.6)

view([-81.6598139306647 13.4415460949037]);

xlim([T(1) T(end)]);
ylim([V(1) V(end)]);

xlb = xlabel('Temperature [°C]');
set(xlb,'rotation',49)
ylabel('Voltage [V]');
zlabel('Power [W]');
ax = gca;
ax.FontSize = 12;

title("Comparison Exponential leakage (Blue) and Linear leakage (Magenta)", 'FontSize', 20);

%view([-80.5895672336695 16.4241547850151])

subplot(3, 3, 7)

surf(T, V, res)
hold;

%surf(T, V, repelem(Power_dyn_low, 1,length(T)), 'FaceColor', 'y', ... %"#A2142F", 
%	'FaceAlpha',0.3)

surf(T, V, plane, 'FaceColor', 'm', ... %"#A2142F", 
	'FaceAlpha',0.6)

fnz= 13;

xlabel('Temperature [°C]', 'FontSize', fnz);
ylabel('Voltage [V]', 'FontSize', fnz);
zlabel('Power [W]', 'FontSize', fnz);
ax = gca;
ax.FontSize = fnz;
xlim([T(1) T(end)]);
ylim([V(1) V(end)]);

title({"Top View:"; "Comparison Exponential and Linear leakage"}, 'FontSize', 15);

view(0,90)
%view(0,0)
%view(90,0)


%%%figure()

subplot(3, 3, 3);

surf(T, V, res, 'FaceAlpha',1, 'EdgeAlpha', 0.6);%0.6)
hold;

h = gca;

for i=1:length(T)
	planel(:,i) = Power_dyn_low;
	planeh(:,i) = Power_dyn_high;
end

%plot3(h.XLim(1):(h.XLim(2) - h.XLim(1))/(length(V)-1):h.XLim(2), V, Power_dyn_high, 'r', 'LineWidth', 1.2);

for i=1:length(T)
	plot3(T(i)*ones(size(V)), V, Power_dyn_low, 'y', 'LineWidth', 1.2);
	plot3(T(i)*ones(size(V)), V, Power_dyn_high, 'r', 'LineWidth', 1.2);
	plot3(T(i)*ones(size(V)), V, Power_static, 'm', 'LineWidth', 1.2);
end
%surf(T, V, planel, 'FaceColor', "#77AC30")
%surf(T, V, planeh, 'FaceColor', "#EDB120", 'FaceAlpha',0.2)

xlabel('Temperature [°C]');
ylabel('Voltage [V]');
zlabel('Power [W]');
ax = gca;
ax.FontSize = fnz;
xlim([T(1) T(end)]);
ylim([V(1) V(end)]);

title({"Comparison Exponential leakge"; "with Max and min Dynamic Power"}, 'FontSize', 15);

view(90,0)

subplot(3, 3, [6, 9]);

surf(T, V, res, 'FaceAlpha',1, 'EdgeAlpha', 0.6);%0.6)
hold;

h = gca;
surf([20 100], V, [Power_dyn_low(:), Power_dyn_low(:)], 'FaceColor', 'y', 'FaceAlpha', 0.8);
surf([20 100], V, [Power_dyn_high(:), Power_dyn_high(:)], 'FaceColor', 'r', 'FaceAlpha', 0.8);
xlb = xlabel('Temperature [°C]');
set(xlb,'rotation',-35)
ylb = ylabel('Voltage [V]');
set(ylb,'rotation',40)
zlabel('Power [W]');
view([-138.844996600224 39.7336634094982]);
ax = gca;
ax.FontSize = fnz;
xlim([T(1) T(end)]);
ylim([V(1) V(end)]);

%title({"Comparison Exponential leakge"; "with Max and min Dynamic Power"}, 'FontSize', 15);

%%figure()

subplot(3, 3, 8);
surf([20 100], V(end-lij: end), [Power_dyn_high(end-lij:end), Power_dyn_high(end-lij:end)], 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.2);
hold on;
surf(T, V, res, 'FaceAlpha',1, 'EdgeAlpha', 0.6);%0.6)


%surf(T, V, planel, 'FaceColor', "#77AC30")
%surf(T, V, planeh, 'FaceColor', "#EDB120", 'FaceAlpha',0.2)

h = gca;

lij = 14;

surf([20 100], V, [Power_dyn_low(:), Power_dyn_low(:)], 'FaceColor', 'y', 'FaceAlpha', 0.8, 'EdgeAlpha', 0.3);

%{
for i=1:length(Power_dyn_high)
	plot3([20 100], [V(i) V(i)], [Power_dyn_low(i) Power_dyn_low(i)], 'y', 'LineWidth', 0.8);
	plot3([20 100], [V(i) V(i)], [Power_dyn_high(i) Power_dyn_high(i)], 'r', 'LineWidth', 0.8);
end
%}

xlabel('Temperature [°C]');
ylabel('Voltage [V]');
zlabel('Power [W]');
ax = gca;
ax.FontSize = fnz;
xlim([T(1) T(end)]);
ylim([V(1) V(end)]);

title({"Comparison Exponential leakge"; "with Max and min Dynamic Power"}, 'FontSize', 15);

view(0,0)

%%
%{
if isunix
	path_name = "/tmp/MATLAB-Figures/expleak";
else
	path_name = "C:\temp\MATLAB-Figures\expleak";
end

bmpres = 1;

% Mange screens
graph_dpi = 416;
file_ext = ".png";
%bmpres = 1; %1.5;

%save as fig
saveas(fig, strcat(path_name,".fig"));
%change res for bmp:
fig.Visible = 'off';
fig.Position = [1, 1, 1920*bmpres,1080*bmpres];
exportgraphics(fig, strcat(path_name,file_ext), 'Resolution', graph_dpi);
close(fig);	
%}

%%

figure();
subplot(1,2,1);

Fvect = hpc.FV_table(1,2):hpc.F_discretization_step:hpc.FV_table(end,3);
Vvect = hpc.FV_table(sum(Fvect > hpc.FV_table(:,3))+1,1);

Fvect = Fvect';

Fvect2 = ones(hpc.FV_levels,1)*Fvect';

for i=1:hpc.FV_levels
	Fvect2(i,Fvect2(i,:) > hpc.FV_table(i,3)) = 0;
end
T = ones(length(Fvect),1)*(70+273.15);
	
d_i = ones(length(Fvect),1)*[0,0,0.25,0.75,0];
d_p = ones(length(Fvect),1);

h = [];
for j=1:hpc.FV_levels
	pup = hpc.power_compute(Fvect2(j,:)', hpc.FV_table(j,1)*ones(size(Fvect2,2),1), T, d_i,d_p, 0);
	%pup = hpc.power_compute(Fvect, Vvect, T, d_i,d_p, 0);
	pup(Fvect2(j,:)'<=0)=NaN;
	hp = plot(Fvect,pup,'LineWidth', 1.5);
	h = [h hp];
	hold on;
end	
pup = hpc.power_compute(Fvect, Vvect, T, d_i,d_p, 0);
hp = plot(Fvect,pup, 'LineWidth', 3, 'Color', 'b');
h = [h hp];

legend(h(end), "Max F-V line", 'FontSize',14);

xlim([Fvect(1), Fvect(end)]);
hold on;
fnz = 16;
xlabel('Frequency [GHz]', 'FontSize',14);
ylabel('Power [W]','FontSize',14);
ax = gca;
ax.FontSize = fnz;
grid on;

title('P-F Relation @Voltage \in [0.5; 1.2]V', 'FontSize',24 );
%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2);

Fvect = hpc.FV_table(1,2):hpc.F_discretization_step*2:hpc.FV_table(end,3);
Vvect = hpc.FV_table(sum(Fvect > hpc.FV_table(:,3))+1,1);

hold on;

Fvect = Fvect';

%Fvect = Fvect3;
%Vvect = Vvect3;

T = ones(hpc.FV_levels,1)*(70+273.15);
d_i = ones(hpc.FV_levels,1)*[0,0,0.25,0.75,0];
d_p = ones(hpc.FV_levels,1);

Fvect2 = ones(hpc.FV_levels,1)*Fvect';
for i=1:hpc.FV_levels
	Fvect2(i,Fvect2(i,:) > hpc.FV_table(i,3)) = 0;
end

h =[];
for j=1:size(Fvect2,2)
	pup = hpc.power_compute(Fvect2(:,j), hpc.FV_table(:,1), T, d_i,d_p, 0);
	%pup = hpc.power_compute(Fvect, Vvect, T, d_i,d_p, 0);
	pup(Fvect2(:,j)'<=0)=NaN;
	hp = plot(hpc.FV_table(:,1),pup, 'LineWidth', 1.5);
	h = [h hp];
	hold on;
end	

%Vvect3 = hpc.FV_table(:,1);
%Fvect3 = hpc.FV_table(:,3);
Fvect3 = hpc.FV_table(1,2):hpc.F_discretization_step:hpc.FV_table(end,3);
Vvect3 = hpc.FV_table(sum(Fvect3 > hpc.FV_table(:,3))+1,1);
T = ones(length(Vvect3),1)*(70+273.15);
d_i = ones(length(Vvect3),1)*[0,0,0.25,0.75,0];
d_p = ones(length(Vvect3),1);

pup = hpc.power_compute(Fvect3', Vvect3, T, d_i,d_p, 0);
hp = plot(Vvect3,pup, 'LineWidth', 3, 'Color', 'b');
h = [h hp];

legend(h(end), "Max F-V line", 'FontSize',14)

xlim([Vvect3(1), Vvect3(end)]);

xlabel('Voltage [V]', 'FontSize',14);
ylabel('Power [W]', 'FontSize',14);
ax = gca;
ax.FontSize = fnz;
title('P-V Relation @Frequency \in [0.4; 3.6]GHz', 'FontSize',24);

grid on;
%%
figure();

Vvect3 = hpc.FV_table(:,1);
Fvect3 = hpc.FV_table(:,3);
T = ones(length(Vvect3),1)*(40+273.15);


d_i = ones(length(Vvect3),1)*[0,0,0,0,1];
d_p = ones(length(Vvect3),1);

pup = hpc.power_compute(Fvect3, Vvect3, T, d_i,d_p, 0);
plot(Fvect3,pup);

hold on;
d_i = ones(length(Vvect3),1)*[0,1,0,0,0];
d_p = ones(length(Vvect3),1);

pup = hpc.power_compute(Fvect3, Vvect3, T, d_i,d_p, 0);
plot(Fvect3,pup);

xlim([Fvect3(1), Fvect3(end)]);

xlabel('Frequency [GHz]', 'FontSize',14);
title('P-F Relation @Workload varying [INT; VECT]', 'FontSize',20 );
ylabel('Power [W]','FontSize',14);
