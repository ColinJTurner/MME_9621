%% Q1
clear all; close all; clc;
% Given parameters
H = 0.1;rho = 1.2;mu = 4e-5;
D = 0.1;n =9;K = 12 * D * mu / H^2;
delY2 = (H / (n + 1))^2;B = -12*D*delY2/H^2;
% Define the matrix A
A = diag(-2 * ones(n, 1)) + diag(ones(n - 1, 1), 1) + diag(ones(n - 1, 1), -1);
% Define the vector b
b = B * ones(n, 1);
% Solve the system of equations
tic;
u_numerical = A \ b;
u_numerical = [0; u_numerical; 0];
time_numerical = toc;
% Define the y values corresponding to each node
y_values = linspace(-H/2, H/2, n+2)';
% Analytical solution
tic;
u_analytical = 3 * D / 2 * (1 - (y_values / (H/2)).^2);
time_analytical = toc;
% bvp4c Solver
bvpfcn = @(y,u) [u(2); -12 * D / H^2];
bcfcn = @(ua,ub) [ua(1); ub(1)];
guess = @(y) [0; 0];
solinit = bvpinit(linspace(-H/2, H/2, n+2), guess);
tic;
sol = bvp4c(bvpfcn, bcfcn, solinit);
time_bvp4c = toc;
u_bvp4c = deval(sol, y_values);
% Plotting
figure;
plot(y_values, u_numerical, '-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName','Numerical Solution');
hold on;
plot(y_values, u_analytical, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName','Analytical Solution');
hold on;
plot(y_values, u_bvp4c(1,:), '-o','LineWidth', 2, 'MarkerSize', 15, 'DisplayName','bvp4c Solution');
xlabel('y');
ylabel('Velocity (u)');
title('Velocity Profile');
legend;
grid on; hold off;
% Evaluate the velocity at y=0 from the numerical solution
u_max_numerical = u_numerical((n+3)/2);
% Calculate the maximum velocity using the given formula
D_max = 3 * D / 2;
error_numerical = norm(u_numerical - u_analytical, 'inf'); % Maximum absolute error
error_bvp4c = norm(u_bvp4c(1,:) - u_analytical', 'inf'); % Maximum absolute error
% Display results analytical
fprintf('Analytical Method:\n');
for i = 1:1:7
 fomatSpec='Node 1 Velocity %d = %12.5f\n';
 fprintf(fomatSpec,i,u_analytical(i));
end
fprintf('Maximum velocity (analytical): %.6f m/s\n', D_max);
fprintf('Computational time: %.6f seconds\n', time_analytical);
% Display results Numerical
fprintf('Numerical Method:\n');
for i = 1:1:7
 fomatSpec='Node 1 Velocity %d = %12.5f\n';
 fprintf(fomatSpec,i,u_numerical(i));
end
fprintf('Maximum velocity (numerical): %.6f m/s\n', u_max_numerical);
fprintf('Maximum absolute error: %12.5e\n', error_numerical);
fprintf('Computational time: %.6f seconds\n', time_numerical);
% Display results bvp4c
fprintf('\nbvp4c Solver:\n');
for i = 1:1:7
 fomatSpec='Node 1 Velocity %d = %12.5f\n';
 fprintf(fomatSpec,i,u_bvp4c(1,i));
end
fprintf('Maximum velocity (bvp4c): %.6f m/s\n', u_bvp4c(1,(n+3)/2));
fprintf('Maximum absolute error: %12.5e\n', error_bvp4c);
fprintf('Computational time: %.6f seconds\n', time_bvp4c);
%% Q2
clc; clear all; close all;
% Given parameters
delt=0.01; v=0.000217; h=40e-3;
yL=0; time_final=2; my = 5;
dely=h/my; u0 = 40;
nt=time_final/delt;dely2=dely*dely;
beta=v*delt/dely2;
Tini = 4e-2;
y = linspace(0, Tini, 6);
t = linspace(0, 2, 200);
xmt = y';
uM=0;
f_initial=@(y) (0);
i=1:(my-1);
u_initial(i)=f_initial(y(i+1));
% explicit method
tic;
[usol_1]=Explicit_1D(beta,u0,uM,my,nt,u_initial);
explicit_time = toc;
% implicit method
tic;
[usol_2]=Implicit_1D(beta,u0,uM,my,nt,u_initial);
implicit_time = toc;
% crank-nicolson method
tic;
[usol_3]=CN_1D(beta,u0,uM,my,nt,u_initial);
C_N_time = toc;
% pdepe solver
tic
sol = pdepe(0, @pdefun1, @pdeIC1, @pdeBC1, y, t, [], v, Tini, u0);
time_pdepe = toc;
fprintf('Comp Time %es\n', time_pdepe);
fprintf('Explicit Method Comp Time: %fs\n', explicit_time);
fprintf('Implicit Method Comp Time: %fs\n', implicit_time);
fprintf('Crank-Nicolson Method Comp Time: %fs\n', C_N_time);
% Boundary Points
for jj=2:nt
 u_1(:,jj)=[u0;usol_1(:,jj);uM];
 u_2(:,jj)=[u0;usol_2(:,jj);uM];
 u_3(:,jj)=[u0;usol_2(:,jj);uM];
end
% plotting explicit method
figure(1);
plot(u_1(:,10),y,'-o',u_1(:,30),y,'-o',u_1(:,60),y,'-o',u_1(:,end),y,'-o');
legend('t=0.1','t=0.3','t=0.6','t=2');
xlabel('Velocity (m/s)');
ylabel('Height (m)');
title('Explicit method');
% plotting implicit method
figure(2);
plot(u_2(:,10),y,'-o',u_2(:,30),y,'-o',u_2(:,60),y,'-o',u_2(:,end),y,'-o');
legend('t=0.1','t=0.3','t=0.6','t=2');
xlabel('Velocity (m/s)');
ylabel('Height (m)');
title('Implicit method');
% plotting C-N method
figure(3);
plot(u_3(:,10),y,'-o',u_3(:,30),y,'-o',u_3(:,60),y,'-o',u_3(:,end),y,'-o');
legend('t=0.1','t=0.3','t=0.6','t=2');
xlabel('Velocity (m/s)');
ylabel('Height (m)');
title('Crank-Nicolson method');
figure(4);
plot(sol(10,:), xmt, '-o', sol(30,:), xmt, '-o', sol(60,:), xmt, '-o', sol(end,:),xmt, '-o');
xlim([0, u0]);
legend('t=0.1', 't=0.3', 't=0.6', 't=2');
xlabel('Velocity (m/s)');
ylabel('Height (m)');
title('Couette Flow pdepe solver');
%Accuracy comparison
figure(5);
timesteps = [10, 30, 60, 200];
for idx = 1:length(timesteps)
 subplot(2, 2, idx);
 hold on;
 plot(u_1(:, timesteps(idx)), y, '-o', 'LineWidth', 1.5,'MarkerSize', 5);
 plot(u_2(:, timesteps(idx)), y, '-o', 'LineWidth', 1.5,'MarkerSize', 10);
 plot(u_3(:, timesteps(idx)), y, '-o', 'LineWidth', 1.5,'MarkerSize', 15);
 hold off;
 title(['Timestep: ', num2str(timesteps(idx))]);
 legend('Explicit', 'Implicit', 'Crank-Nicolson');
 xlabel('Velocity (m/s)'); ylabel('Height (m)');
end
% error calculation at each point compared to pdepe solver
u_1error(1) = norm(u_1(:,10) - sol(30,:), 'inf'); % Maximum absolute error
u_2error(1) = norm(u_2(:,10) - sol(30,:), 'inf'); % Maximum absolute error
u_3error(1) = norm(u_3(:,10) - sol(30,:), 'inf'); % Maximum absolute error
u_1error(2) = norm(u_1(:,30) - sol(30,:), 'inf'); % Maximum absolute error
u_2error(2) = norm(u_2(:,30) - sol(30,:), 'inf'); % Maximum absolute error
u_3error(2) = norm(u_3(:,30) - sol(30,:), 'inf'); % Maximum absolute error
u_1error(3) = norm(u_1(:,end) - sol(end,:), 'inf'); % Maximum absolute error
u_2error(3) = norm(u_2(:,end) - sol(end,:), 'inf'); % Maximum absolute error
u_3error(3) = norm(u_3(:,end) - sol(end,:), 'inf'); % Maximum absolute error
fprintf('Error Explicit Method\n');
for i=1:1:3
%formatpec = '%5f\n';
fprintf('%5f\n',u_1error(i));
end
fprintf('Error Implicit Method\n');
for i=1:1:3
%formatpec = '%5f\n';
fprintf('%5f\n',u_2error(i));
end
fprintf('Crank-Nicolson Method\n');
for i=1:1:3
%formatpec = '%5f\n';
fprintf('%5f\n',u_3error(i));
end
% explicit method function
function [u]=Explicit_1D(beta,u0,uM,my,nt,u_initial)
A=zeros(my-1,my-1);
u(:,1)=u_initial';
S=[u0; zeros(my-3,1); uM]; % BC
for i=1:my-1
 A(i,i)=1-2*beta;
 if i<my-1
 A(i,i+1)=beta;
 end
 if i>=2
 A(i,i-1)=beta;
 end
end
for n=2:nt % time loop
 u(:,n)=A*u(:,n-1)+beta*S;
end
end
% implicit method function
function [u]=Implicit_1D(beta,u0,uM,my,nt,u_initial)
A=zeros(my-1,my-1);
u(:,1)=u_initial';
S=[u0; zeros(my-3,1); uM]; %BC
for i=1:my-1
 A(i,i)=-(1+2*beta);
 if i<my-1
 A(i,i+1)=beta;
 end
 if i>=2
 A(i,i-1)=beta;
 end
end
for n=2:nt % time loop
 b=-u(:,n-1)-beta*S;
 u(:,n)=A\b;
end
end
% crank nicolson method funciton defined
function [u]=CN_1D(beta,u0,uM,my,nt,u_initial)
A=zeros(my-1,my-1); B=A;
u(:,1)=u_initial';
S=[u0; zeros(my-3,1); uM]; %BC
for i=1:my-1
 A(i,i)=-2*(1+beta);
 B(i,i)=-2*(1-beta);
 if i<my-1
 A(i,i+1)=beta;
 B(i,i+1)=-beta;
 end
 if i>=2
 A(i,i-1)=beta;
 B(i,i-1)=-beta;
 end
end
for n=2:nt %time loop
 b=B*u(:,n-1)-beta*(S+S);
 u(:,n)=A\b;
end
end
function [c, f, s] = pdefun1(y, t, u, dudy, v, Tini, u0)
 c = 1/v;
 f = dudy;
 s = 0;
end
function u0 = pdeIC1(y, v, Tini, u0)
 u0 = zeros(size(y));
 u0(1) = u0;
end
function [pl, ql, pr, qr] = pdeBC1(yL, uL, yR, uR, t, v, Tini, u0)
 pl = uL - u0; ql = 0;
 pr = uR; qr = 0;
end