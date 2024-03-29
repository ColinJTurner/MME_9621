%Example-2
clc; clear all;
Fode=@(t,y)[y(2);y(3);-3.*y(1).*y(3)+2.*y(2).^2-y(4);y(5);-3.*0.7.*y(1).*y(5)]; % anonymous function
tspan=[0 5]; % time interval
y0=[0 0 0.68 1 -0.5]; % initial condition 
[t y]=ode45(Fode,tspan,y0,[ ]); % anonymous function


figure(1);
plot(t,y(:,1)); % plots t vs x % time series plot
hold on;
plot(t,y(:,2)); % plots t vs x % time series plot
hold on;
plot(t,y(:,3)); % plots t vs x % time series plot
hold on;
plot(t,y(:,4)); % plots t vs x % time series plot
hold on;
plot(t,y(:,5)); % plots t vs x % time series plot
hold off;

xlabel('n'); ylabel('y');