clear all; close all; clc;
% Define circuit parameters
R = 100; % 100 ohms
R2 = 100;
L1 = 0.1; % 0.1 H
L2 = 1; % 1 H (second inductor)
M = 1; % 1 Mutual inductance (H)
C1 = 0.04; % 0.04 F
C2 = 0.01; % 0.01 F

% Define time varying power source
V = @(t) 50 * sin(2*pi*60*t);

% Define initial conditions
i0 = 0.5; % 0.5 A
vc0 = 0; % 0 V
i20 = 0.2; % 0.2 Initial current for the second inductor
vc1 = 0; % 0 V
initial_conditions = [i0; vc0; i20; vc1];

% Define time span
tspan = [0 0.1]; % Solve from t=0 to t=0.1 seconds
% h=1e-5;
h=tspan(2)*10^-4;
NStep = abs(tspan(1) - tspan(2)) / h;
% Solve using ODE45
tic;
[t_ode45, y_ode45] = ode45(@(t, y) circuit_odes(t, y, V, R, L1, L2, M, C1, C2), tspan, initial_conditions);
comp_ode45 = toc;
tic;
[t_ode15, y_ode15] = ode15s(@(t, y) circuit_odes(t, y, V, R, L1, L2, M, C1, C2), tspan, initial_conditions);
comp_ode15 = toc;
tic;
[t_rk4, y_rk4] = Vector_RK4(@circuit_odes, tspan(1), initial_conditions, NStep, h, V, R, L1, L2, M, C1, C2);
comp_rk4 = toc;

% Interpolate ODE45 solution onto ODE15s time grid
y_ode45_interp = interp1(t_ode45, y_ode45, t_ode15);
y_rk4_interp = interp1(t_rk4, y_rk4, t_ode15);

% Calculate absolute error for ODE15s
abs_error_ode15_ode45 = abs(y_ode15 - y_ode45_interp);

% Calculate absolute error for RK4
abs_error_rk4_ode15 = abs(y_ode15 - y_rk4_interp);

% Calculate absolute error between ODE45 and RK4
abs_error_ode45_rk4 = abs(y_ode45_interp - y_rk4_interp);


disp(['Computation Time for ode45: ', num2str(comp_ode45), ' seconds']);
disp(['Computation Time for ode15s: ', num2str(comp_ode15), ' seconds']);
disp(['Computation Time for rk4: ', num2str(comp_rk4), ' seconds']);
fprintf('Mean error for ode15s & ode45 for: %e\n', max(abs_error_ode15_ode45));
fprintf('Mean error for ode15s & rk4: %e\n', max(abs_error_rk4_ode15));
fprintf('Mean error for ode45 & rk4: %e\n', max(abs_error_ode45_rk4));

% Plot results ode45
figure;
subplot(2,1,1);
plot(t_ode45, y_ode45(:,1), 'b', t_ode45, y_ode45(:,3), 'r');
title('Currents');
xlabel('Time (s)');
ylabel('Current (A)');
legend('i1', 'i2');
grid on;

subplot(2,1,2);
plot(t_ode45, y_ode45(:,2), 'b', t_ode45, y_ode45(:,4), 'r');
title('Capacitor Voltages');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Vc1', 'Vc2');
grid on;

% Plot results ode15s
figure;
subplot(2,1,1);
plot(t_ode15, y_ode15(:,1), 'b', t_ode15, y_ode15(:,3), 'r');
title('Currents (ode15s)');
xlabel('Time (s)');
ylabel('Current (A)');
legend('i1', 'i2');
grid on;

subplot(2,1,2);
plot(t_ode15, y_ode15(:,2), 'b', t_ode15, y_ode15(:,4), 'r');
title('Capacitor Voltages (ode15s)');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Vc1', 'Vc2');
grid on;

% Plot results RK4
figure;
subplot(2,1,1);
plot(t_rk4, y_rk4(:,1), 'b', t_rk4, y_rk4(:,3), 'r');
title('Currents (RK4)');
xlabel('Time (s)');
ylabel('Current (A)');
legend('i1', 'i2');
grid on;

subplot(2,1,2);
plot(t_rk4, y_rk4(:,2), 'b', t_rk4, y_rk4(:,4), 'r');
title('Capacitor Voltages (RK4)');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Vc1', 'Vc2');
grid on;

% Plot absolute error
% Plot results
figure;
subplot(2,1,1);
plot(t_ode15, abs_error_ode15_ode45(:,1), 'b', t_ode15, abs_error_rk4_ode15(:,1), 'r', t_ode15, abs_error_ode45_rk4(:,1), 'g');
title('Absolute Error Between Different Methods (Currents)');
xlabel('Time (s)');
ylabel('Absolute Error');
legend('ODE15s vs ODE45', 'RK4 vs ODE15s', 'RK4 vs ODE45');
grid on;

subplot(2,1,2);
plot(t_ode15, abs_error_ode15_ode45(:,2), 'b', t_ode15, abs_error_rk4_ode15(:,2), 'r', t_ode15, abs_error_ode45_rk4(:,2), 'g');
title('Absolute Error Between Different Methods (Capacitor Voltages)');
xlabel('Time (s)');
ylabel('Absolute Error');
legend('ODE15s vs ODE45', 'RK4 vs ODE15s', 'RK4 vs ODE45');
grid on;