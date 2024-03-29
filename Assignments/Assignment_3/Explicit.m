
clear all; close all; clc;
% Given parameters
h = 0.04; % plate spacing in meters (40 mm)
v = 0.000217; % kinematic viscosity in m^2/s
u0 = 40; % initial velocity in m/s
M = 5; % number of space nodes
T = 1; % total simulation time
dt = 0.001; % time step size
dy = h / (M - 1); % spatial grid size

% Initialize velocity matrix
u = zeros(M, 1);
u(1) = u0; % initial condition at y = 0

% Time stepping
for t = dt:dt:T
    % Calculate alpha
    alpha = v * dt / dy^2;
    
    % Create coefficient matrix A
    A = diag(1 - 2 * alpha * ones(M, 1)) + diag(alpha * ones(M-1, 1), 1) + diag(alpha * ones(M-1, 1), -1);
    A(1, 1) = 1; % Boundary condition at y = 0
    A(end, end) = 1; % Boundary condition at y = h
    
    % Create right-hand side vector b
    b = u;
    
    % Solve for u at the next time step
    u = A \ b;
end

% Display the final velocity profile
disp('Final velocity profile:');
disp(u);