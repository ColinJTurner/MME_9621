close all; clear all; clc;

% Define circuit parameters
R = 100; % ohms
L1 = 0.1; % H
L2 = 0.05; % H (second inductor)
M = 0.02; % Mutual inductance (H)
C = 0.01; % F

% Define initial conditions
i0 = 0.5; % A
vc0 = 0; % V
i20 = 0.2; % Initial current for the second inductor

% Define time span
tspan = [0 0.1]; % Solve from t=0 to t=0.1 seconds

% Define excitation voltage function
V = @(t) 10 * sin(2*pi*50*t); % V(t) = 10*sin(2*pi*50*t)

% Define the differential equation system with mutual inductance
ode_system = @(t, Y) [
    (1/(L1-M))* (V(t) - R*Y(1) - (L2-M)*Y(3) - Y(2));   % di1/dt = (1/(L1-M))*(V - Ri1 - (L2-M)i2 - vc)
    (1/C)*Y(1);                                         % dvc/dt = i1/C
    (1/(L2-M))*(Y(2) - R*Y(3) - (L1-M)*Y(1) - Y(4));    % di2/dt = (1/(L2-M))*(vc - Ri2 - (L1-M)i1 - vc2)
    (1/C)*Y(3)                                          % dvc2/dt = i2/C
    ];

% Compute the Jacobian matrix for eigenvalue analysis
Jacobian_matrix = @(t, Y) [
    -R/(L1-M),         -1/(L1-M),          -(L2-M)/(L1-M),          0;
    1/C,               0,                  0,                       0;
    -(L1-M)/(L2-M),    0,                  -R/(L2-M),               -1/(L2-M);
    0,                 0,                  1/C,                     0
    ];

% Solve the differential equations using ode45
[t_ode45, Y_ode45] = ode45(ode_system, tspan, [i0 vc0 i20 vc0]);

% Solve the differential equations using ode15s
options = odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
[t_ode15s, Y_ode15s] = ode15s(ode_system, tspan, [i0 vc0 i20 vc0], options);

h=1e-6;
% Solve the differential equations using RK4
% 1. Note the presence of eigenvalues with magnitudes differing by several orders 
% of magnitude (e.g., -3.3e+03 and -1.0) suggests stiffness in the system. 
% When eigenvalues exhibit such a wide range, it indicates that the system's 
% dynamics vary significantly, which can make it challenging for numerical 
% solvers to handle.
% 2. Stiffness can affect the performance of numerical methods. Traditional
% explicit methods like RK4 might require very small step sizes to handle
% stiff systems accurately, which can significantly increase computational
% cost. Implicit methods or stiff solvers like ode15s are better suited
% for stiff systems as they can handle the wide range of timescales more efficiently.
NStep = round((tspan(2) - tspan(1)) / h); % Make sure NStep is a scalar
[t_vecRK4, Y_vecRK4]=VectorRK4(ode_system,tspan,[i0 vc0 i20 vc0],NStep,h);

% Compute eigenvalues of the Jacobian matrix at a specific point
Y_test = [i0 vc0 i20 vc0];  % Initial conditions
eigenvalues = eig(Jacobian_matrix(0, Y_test));  % Compute eigenvalues at t=0

% Display the eigenvalues
disp('Eigenvalues of the Jacobian matrix:');
disp(eigenvalues);

% Extract current (i1, i2) and voltage across the capacitor (vc)
i1_ode45 = Y_ode45(:,1);
vc_ode45 = Y_ode45(:,2);
i2_ode45 = Y_ode45(:,3);
vc2_ode45 = Y_ode45(:,4);

i1_ode15s = Y_ode15s(:,1);
vc_ode15s = Y_ode15s(:,2);
i2_ode15s = Y_ode15s(:,3);
vc2_ode15s = Y_ode15s(:,4);

i1_vecRK4 = Y_vecRK4(:, 1);
vc_vecRK4 = Y_vecRK4(:, 2);
i2_vecRK4 = Y_vecRK4(:, 3);
vc2_vecRK4 = Y_vecRK4(:, 4);

% Plot the results
figure;
subplot(2,1,1);
plot(t_ode45, i1_ode45, t_ode45, i2_ode45);
hold on;
plot(t_ode15s, i1_ode15s, '--', t_ode15s, i2_ode15s, '--');
% hold on;
% plot(t_vecRK4, i1_vecRK4, '-o', t_vecRK4, i2_vecRK4, '-o');
hold off;
% legend('Current i1 (ode45)', 'Current i2 (ode45)', 'Current i1 (ode15s)',...
%     'Current i2 (ode15s)', 'Current i1 (vecRK4)', 'Current i2 (vecRK4)');
legend('Current i1 (ode45)', 'Current i2 (ode45)', 'Current i1 (ode15s)',...
    'Current i2 (ode15s)');
xlabel('Time (s)');
ylabel('Current (A)');
title('Currents vs Time');

subplot(2,1,2);
plot(t_ode45, vc_ode45, t_ode45, vc2_ode45);
hold on;
plot(t_ode15s, vc_ode15s, '--', t_ode15s, vc2_ode15s, '--');
% hold on;
% plot(t_vecRK4, vc_vecRK4, '-o', t_vecRK4, vc_vecRK4, '-o');
hold off;
% legend('Voltage across Capacitor (ode45)', 'Voltage across Capacitor 2 (ode45)',...
%     'Voltage across Capacitor (ode15s)', 'Voltage across Capacitor 2 (ode15s)',...
%     'Voltage across Capacitor (vecRK4)', 'Voltage across Capacitor 2 (vecRK4)');
legend('Voltage across Capacitor (ode45)', 'Voltage across Capacitor 2 (ode45)',...
    'Voltage across Capacitor (ode15s)', 'Voltage across Capacitor 2 (ode15s)');
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Voltages across Capacitors vs Time');


% Plot the results for VectorRK4
figure;
subplot(2,1,1);
plot(t_vecRK4, i1_vecRK4, t_vecRK4, i2_vecRK4);
xlabel('Time (s)');
ylabel('Current (A)');
title('Currents vs Time (VectorRK4)');
legend('Current i1', 'Current i2');

subplot(2,1,2);
plot(t_vecRK4, vc_vecRK4, t_vecRK4, vc2_vecRK4);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Voltages across Capacitors vs Time (VectorRK4)');
legend('Voltage across Capacitor', 'Voltage across Capacitor 2');
