%% Q1 
clc;clear all;format long e;
[Aug1, Text]=xlsread('DataFile_Assn2.xlsx');
A=Aug1(:,1:21); 
b=Aug1(:,22); 

tic;
max_iter=1000;TOL=1e-5;
[xcJ, iterJ]=Jacobi(A,b,max_iter,TOL);
tJ=toc;

tic;
max_iter=1000;TOL=1e-5;
[xcG, iterG]=GaussSeidel(A,b,max_iter,TOL);
tG =toc;

tic;
max_iter=1000;TOL=1e-5;
[xcS0, iterS0]=SOR(A,b,0.5,max_iter,TOL);
tS0=toc;

tic;
max_iter=1000;TOL=1e-5;
[xcS1, iterS1]=SOR(A,b,1.5,max_iter,TOL);
tS1=toc;

tic;
xc=A\b;
tB=toc;

% Residuals
residual_J = norm(b-A*xcJ);
residual_G = norm(b-A*xcG);
residual_S0 = norm(b-A*xcS0);
residual_S1 = norm(b-A*xcS1);
residual_backslash = norm(b-A*xc);

fprintf("Jacobi Method\t\t Gauss Seidel Method\t SOR w < 1 Method\t SOR w > 1 Method\t BackSlash Method\n");
for i = 1:1:21
    fomatSpec='F%d = %12.5f\t F%d = %12.5f\t F%d = %12.5f\t F%d = %12.5f\t F%d = %12.5f\n';
    fprintf(fomatSpec,i,xcJ(i),i,xcG(i),i,xcS0(i),i,xcS1(i),i,xc(i));
end

formatSpec1='time:\nt = %fs\t\t t = %fs\t\t t = %fs\t\t t = %fs\t\t t = %fs\n';
fprintf(formatSpec1,tJ,tG,tS0,tS1,tB);
formatSpec2='Iterations:\niter = %d\t\t iter = %d\t\t iter = %d\t\t iter = %d\n';
fprintf(formatSpec2,iterJ,iterG,iterS0,iterS1);
formatSpec3='Residuals:\nRes = %e\t Res = %e\t Res = %e\t Res = %e\t Res = %e\n';
fprintf(formatSpec3,residual_J,residual_G,residual_S0,residual_S1,residual_backslash);

figure;

% Sparsity patern of matrix representing non zero elements
subplot(1, 2, 1);
spy(A);
title('Structure of the Coefficient Matrix');
xlabel('Column Index');
ylabel('Row Index');
grid on;

% 3d view of plot showing magnitude difference in coefficients
subplot(1, 2, 2);
surf(A);
title('Surface Plot of the Coefficient Matrix');
xlabel('Column Index');
ylabel('Row Index');
zlabel('Value');
%% Q2 
clc; clear all; format long e;
% Given parameters
L=10;Hp=0.05;sig=2.7*10^-9;Ti=200;Ta=300;Tb=400;dx=2;
b=Hp*dx^2;a=2+b;c=sig*dx^2;d=b*Ti+c*Ti^4;e=Ta+d;f=Tb+d;

% Initial guesses
Guess1=[2; 4; 6; 8]; Guess2=[2.1; 4.1; 6.1; 8.1];

% Tolerance and maximum iterations for the solvers
TOL = 0.5e-3; max_iter = 100;

% System of non linear eqn
Fn=@(x,a,b,c,d,f,e) [...
    a.*x(1)+c.*x(1).^4-x(2)-e;...
    -x(1)+a.*x(2)+c.*x(2).^4-x(3)-d;...
    -x(2)+a.*x(3)+c.*x(3).^4-x(4)-d;...
    -x(3)+a.*x(4)+c.*x(4).^4-f];
options=optimset('display','iter');

% Compute the symbolic Jacobian matrix
syms x1 x2 x3 x4;
Fs = [a*x1 + c*x1^4 - x2 - e;
      -x1 + a*x2 + c*x2^4 - x3 - d;
      -x2 + a*x3 + c*x3^4 - x4 - d;
      -x3 + a*x4 + c*x4^4 - f];
Xs = [x1; x2; x3; x4];
DFs = jacobian(Fs, Xs);

% Multivariate Newton
tic;
[xcN, iterN] = multivariateNewton_fs(Fn, a, b, c, d, f, e, @Jac_fs, DFs, Xs, Guess1, TOL, max_iter);
tN=toc;

% Jacobian is assumed to be Identity matrix
A= eye(4);

% Broyden Method
tic
[xcB,iterB]=BroydenMethod1(Fn,a, b, c, d, f, e,Guess1,Guess2,A,TOL,max_iter);
tB=toc;

% fSolve Method
tic
[xcF]=fsolve(Fn,Guess1,options,a,b,c,d,f,e);
tF=toc;

% Residuals
residual_N = norm(Fn(xcN, a, b, c, d, f, e));
residual_B = norm(Fn(xcB, a, b, c, d, f, e));
residual_F = norm(Fn(xcF, a, b, c, d, f, e));

%Display the final solutions, comp times, iterations & Residuals
fprintf("Multi Newton Method\t Broyden Method\t\t fSolve\n");
for i = 1:1:4
    fomatSpec='F%d = %12.5f\t F%d = %12.5f\t F%d = %12.5f\n';
    fprintf(fomatSpec,i,xcN(i),i,xcB(i),i,xcF(i));
end

formatSpec1='time:\nt = %fs\t\t t = %fs\t\t t = %fs\n';
fprintf(formatSpec1,tN,tB,tF);
formatSpec2='Iterations:\niter = %d\t\t iter = %d\t\t iter = %d\n';
fprintf(formatSpec2,iterN,iterB,12);
formatSpec3='Residuals:\nRes = %e\t Res = %e\t Res = %e\n';
fprintf(formatSpec3,residual_N,residual_B,residual_F);

% Temperature distribution plot
x = 0:2:10;
Ttot = [Ta; xcF; Tb];
plot(x, Ttot);
xlabel('Position along the rod (m)');
ylabel('Temperature (K)');
title('Temperature Distribution along the Rod');
grid on;
%% Q3
clc; clear all; close all; format long e;
% Given parameters
m=[12000 10000 8000];k=[3000 2400 1800];
a = k(1)/m(1); b = k(2)/m(1); c = k(2)/m(2);
d = k(3)/m(2); e = k(3)/m(3);

% RK4 parameters
h = 0.1; t0 = 0; tn = 20;
NStep = abs(tn - t0) / h;

% Initial conditions
Y0 = [0; 1; 0; 0; 0; 0];

% RK4 method
tic;
[t_rk4, Y_rk4] = VectorRK4(@ode_function, t0, Y0, NStep, h, a, b, c, d, e);
trk4 = toc;
% ode45 solver
tic;
[t_ode45, Y_ode45] = ode45(@(t, Y) ode_function(t, Y, a, b, c, d, e), [t0, tn], Y0);
tode45 = toc;

% Extract positions and velocities
x1_rk4 = Y_rk4(:, 1);   x1_ode45 = Y_ode45(:, 1);
v1_rk4 = Y_rk4(:, 2);   v1_ode45 = Y_ode45(:, 2);
x2_rk4 = Y_rk4(:, 3);   x2_ode45 = Y_ode45(:, 3);
v2_rk4 = Y_rk4(:, 4);   v2_ode45 = Y_ode45(:, 4);
x3_rk4 = Y_rk4(:, 5);   x3_ode45 = Y_ode45(:, 5);
v3_rk4 = Y_rk4(:, 6);   v3_ode45 = Y_ode45(:, 6);

% Plotting Displacement
figure;
subplot(2, 1, 1);
plot(t_rk4, [x1_rk4, x2_rk4, x3_rk4], 'LineWidth', 1.5);
hold on;
plot(t_ode45, [x1_ode45, x2_ode45, x3_ode45], '--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Displacement');
title('Displacement versus Time');
legend('x1 (RK4)', 'x2 (RK4)', 'x3 (RK4)', 'x1 (ode45)', 'x2 (ode45)', 'x3 (ode45)');
grid on;

% Plotting Velocity
subplot(2, 1, 2);
plot(t_rk4, [v1_rk4, v2_rk4, v3_rk4], 'LineWidth', 1.5);
hold on;
plot(t_ode45, [v1_ode45, v2_ode45, v3_ode45], '--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity');
title('Velocity versus Time');
legend('v1 (RK4)', 'v2 (RK4)', 'v3 (RK4)', 'v1 (ode45)', 'v2 (ode45)', 'v3 (ode45)');
grid on;
hold off;

% Plot the 3D phase plane of displacements
figure;
plot3(x1_rk4, x2_rk4, x3_rk4, 'LineWidth', 1.5);
hold on;
plot3(x1_ode45, x2_ode45, x3_ode45, '--', 'LineWidth', 1.5);
legend('RK4','ode45')
xlabel('x1');
ylabel('x2');
zlabel('x3');
title('3D Phase Plane of Displacements');

%Comp time for each process
fprintf("Computation time:\n");
fomatSpec='Vectorized RK4:%fs\tODE45 solver:%fs\n';
fprintf(fomatSpec,trk4,tode45);

% Calculate the absolute error for displacement components
common_length = min(length(x1_rk4), length(x1_ode45));
abs_error_x1 = abs(x1_rk4(1:common_length) - x1_ode45(1:common_length));
abs_error_x2 = abs(x2_rk4(1:common_length) - x2_ode45(1:common_length));
abs_error_x3 = abs(x3_rk4(1:common_length) - x3_ode45(1:common_length));
% Calculate the absolute error for velocity components
common_length_v = min(length(v1_rk4), length(v1_ode45));
abs_error_v1 = abs(v1_rk4(1:common_length_v) - v1_ode45(1:common_length_v));
abs_error_v2 = abs(v2_rk4(1:common_length_v) - v2_ode45(1:common_length_v));
abs_error_v3 = abs(v3_rk4(1:common_length_v) - v3_ode45(1:common_length_v));

% Display the maximum absolute errors for each displacement component
fprintf('Maximum Absolute Error for x1: %e\n', max(abs_error_x1));
fprintf('Maximum Absolute Error for x2: %e\n', max(abs_error_x2));
fprintf('Maximum Absolute Error for x3: %e\n', max(abs_error_x3));
% Display the maximum absolute errors for each velocity component
fprintf('Maximum Absolute Error for v1: %e\n', max(abs_error_v1));
fprintf('Maximum Absolute Error for v2: %e\n', max(abs_error_v2));
fprintf('Maximum Absolute Error for v3: %e\n', max(abs_error_v3));
