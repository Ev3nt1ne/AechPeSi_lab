nfunction [x,u] = SimulateMSD(sys,cntr, x_init, T, w)
if nargin < 5
    w = zeros(size(x_init,1),T);
end

% Set up variables
x = zeros(sys.n,T+1);
u = zeros(sys.m,T);
x(:,1)=x_init;

for i=1:T
    u(:,i) = cntr(x(:,i));
    x(:,i+1) = sys.step(x(:,i),u(:,i),w(:,i));
end

% ==============
% Plot System Trajectory within Constraints
% ==============

figure()
% Plot trajectory 
plot(x(1,:),x(2,:))
hold on
% Plot state constraints
plot(sys.PxP,'alpha',0.1)
hold off
xlabel('Silicon')
ylabel('Copper')

figure()
plot(x(end,:))
title('Temp Amb')

figure()
subplot(2,1,1),hold on,grid on,title('State'),xlabel('Time (s)')
plot([0:T]',x-ones(sys.n,T+1)*273.15)
subplot(2,1,2),grid on, hold on,title('Input Power'),xlabel('Time (s)')
plot([1:T]',u(1:(end-1),:))

figure()
subplot(3,1,1)
plot(x(1,1:T))
hold on
plot([1,T],[sys.xmax(1), sys.xmax(1)],'k--')
plot([1,T],[sys.xmin(1), sys.xmin(1)],'k--')
hold off
xlim([1,T])
ylabel('Silicon')

subplot(3,1,2)
plot(x(2,:))
hold on
plot([1,T],[sys.xmax(2), sys.xmax(2)],'k--')
plot([1,T],[sys.xmin(2), sys.xmin(2)],'k--')
hold off
xlim([1,T])
ylabel('Copper')

subplot(3,1,3)
plot(u(1,:))
hold on
plot([1,T],[sys.umax(1), sys.umax(1)],'k--')
plot([1,T],[sys.umin(1), sys.umin(1)],'k--')
hold off
xlim([1,T])
ylabel('Input')
xlabel('Time steps')

end

