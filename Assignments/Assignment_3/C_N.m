
clear all; close all; clc;
% Given parameters
h = 0.04; % Plate spacing in meters
v = 0.000217; % Kinematic viscosity in m^2/s
u0 = 40; % Initial velocity in m/s
M = 5; % Number of spatial nodes
T = 10; % Total simulation time
dt = 0.01; % Time step
dy = h / (M + 1); % Spatial step

% Calculate alpha and beta
alpha = v * dt / (dy^2);
beta = v * dt / (2 * dy^2);

% Initialize arrays
u = zeros(M+2, round(T/dt)+1); % Add two ghost nodes for boundary conditions
u(:,1) = u0; % Initial condition

% Construct the coefficient matrix A
A = diag(ones(M,1)*(1+2*alpha)) + diag(-alpha*ones(M-1,1),1) + diag(-alpha*ones(M-1,1),-1);

% Apply boundary conditions to A
A(1,1) = 1 + alpha;
A(end,end) = 1 + alpha;

% Time integration loop
for n = 1:round(T/dt)
    % Construct the right-hand side vector b
    b = zeros(M, 1);
    b(1) = beta * u(1,n) + (1 - 2*beta) * u(2,n) + beta * u(3,n);
    b(end) = beta * u(M+1,n) + (1 - 2*beta) * u(M,n+1) + beta * u(M+2,n);
    
    % Solve the linear system
    u(2:M+1,n+1) = A \ b;
end

% Plot the results
figure;
[X,Y] = meshgrid(0:dy:h,0:dt:T);
surf(X,Y,u');
xlabel('Position (m)');
ylabel('Time (s)');
zlabel('Velocity (m/s)');
title('Velocity Profile over Time');
