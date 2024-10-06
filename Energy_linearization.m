% Define the range for x and y
x = linspace(0.4, 4, 20);
y = linspace(0.4, 4, 20);
[X, Y] = meshgrid(x, y);

% Define the functions
F = Y ./ X;
G = Y - X;
G = 1 - X ./ Y;

% Create the plots
figure;

% Plot f = y/x
subplot(1, 2, 1);
surf(X, Y, F);
title('f = y/x');
xlabel('x');
ylabel('y');
grid on;
hold on;
%plot([min(x) max(x)], [0 0], 'k--'); % x-axis
%plot([0 0], [min(y) max(y)], 'k--'); % y-axis
hold off;

% Plot g = y - x
subplot(1, 2, 2);
surf(X, Y, G);
title('g = y - x');
xlabel('x');
ylabel('y');
grid on;
hold on;
%plot([min(x) max(x)], [0 0], 'k--'); % x-axis
%plot([0 0], [min(y) max(y)], 'k--'); % y-axis
hold off;

% Adjust layout
set(gcf, 'Position', [100, 100, 1200, 500]);

figure;
surf(X, Y, F);
grid on;
hold on;
surf(X, Y, G);

%%



instr = 3;

F = linspace(0.4, 4, 20);

F = F;

time = instr ./ F;
diff = F - F(1);
rat = F./F(1);
ndiff = diff;

show = [time ndiff rat rat.*(1-ndiff)];

x0 = 1; %2 %instr
k = instr;

func = k./x0 - k./(x0.^2).*(F - x0);
func1 = k./x0 - k./(x0.^2).*(F - x0).*(k-F);

show = [time func];

figure, plot(F,time), hold on, grid on, plot(F,func) %plot(F,func1) 


%%

hpc.leak_exp_t_k = 10e-3;
hpc.leak_exp_k = -4.5;


hpc.leakageplot()




%hpc.discrplot()