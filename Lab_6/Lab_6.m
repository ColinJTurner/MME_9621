%Example-2
clc; clear all;
Fode=@(t,y,Pr)[y(2);y(3);-3.*y(1).*y(3)+2.*y(2).^2-y(4);y(5);-3.*Pr.*y(1).*y(5)]; % anonymous function
tspan=[0 5]; % time interval
y0=[0 0 0.68 1 -0.5]; % initial condition
Pr=0.7; % Parameter
[t y]=ode45(Fode,tspan,y0,[ ],Pr); % anonymous function


figure(1);
plot(t,y(:,1)); % Plot row 1
hold on;
plot(t,y(:,2)); % Plot row 2
hold on;
plot(t,y(:,3)); % Plot row 3
hold on;
plot(t,y(:,4)); % Plot row 4
hold on;
plot(t,y(:,5)); % Plot row 5
hold off;

xlabel('n'); ylabel('y');
