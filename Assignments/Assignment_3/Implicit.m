clear all; close all; clc;
% Given parameters
h = 0.04; % plate spacing in meters
v = 0.000217; % kinematic viscosity in m^2/s
u0 = 40; % initial velocity in m/s
M = 5; % number of space nodes
delta_y = h / M; % space step size

% Time parameters
delta_t = 0.1; % time step size
total_time = 10; % total simulation time

% Calculate alpha
alpha = v * delta_t / (delta_y^2);

% Initialize coefficient matrix A and right-hand side vector b
A = zeros(M, M);
b = zeros(M, 1);

% Set up coefficient matrix A
for i = 1:M
    if i == 1
        A(i, i) = 1 + 2 * alpha;
        A(i, i+1) = -1;
    elseif i == M
        A(i, i) = 1 + 2 * alpha;
        A(i, i-1) = -1;
    else
        A(i, i) = 1 + 2 * alpha;
        A(i, i-1) = -1;
        A(i, i+1) = -1;
    end
end

% Set up right-hand side vector b
b(1) = u0;
b(M) = 0;

% Initialize solution vector
u = zeros(M, 1);
u_new = zeros(M, 1);

% Time loop
num_steps = total_time / delta_t;
for t = 1:num_steps
    % Solve the tridiagonal system of equations
    u_new = A \ b;
    
    % Update boundary conditions
    u_new(1) = u0;
    u_new(M) = u0 * (1 - delta_y / h);
    
    % Update u for the next time step
    u = u_new;
end

% Display the final solution
disp('Final solution:');
disp(u);