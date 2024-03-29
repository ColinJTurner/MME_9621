% clc; 
% clear all; 
% close all;
% 
% % Input Parameters
% H = 4e-2;   % Spacing between the plates (m)
% v = 2.17e-4;% Kinematic viscosity (m^2/s)
% u0 = 40;    % Velocity of the moving plate (m/s)
% uR = 0;     % Velocity of the stationary plate (m/s)
% 
% % Grid
% xL = 0;                 % Left boundary
% xR = H;                 % Right boundary
% xm = linspace(xL, xR, 6);   % Spatial grid
% t0 = 0;                 % Initial time
% tend = 2;               % Final time
% tm = linspace(t0, tend, 200); % Time grid
% xmt = xm';
% 
% % Solve the PDE using pdepe
% tic
% sol = pdepe(0, @pdefun1, @pdeIC1, @pdeBC1, xm, tm, [], v, H, u0);
% time_pdepe = toc;
% fprintf('Time required to solve the problem is %e seconds.\n\n', time_pdepe);
% 
% % Post-processing and Plotting
% figure(1); % Plot (x-y plot)
% plot(xmt, sol(10,:), '-o', xmt, sol(30,:), '-o', xmt, sol(60,:), '-o', xmt, sol(end,:), '-o');
% ylim([0, 40]); 
% legend('t=0.1', 't=0.3', 't=0.6', 't=2');
% xlabel('Height (mm)');
% ylabel('Velocity (m/s)');
% title('Velocity Profile for Couette Flow');
% 
% % Define the PDE for Couette flow
% function [c, f, s] = pdefun1(x, t, u, dudx, v, H, u0)
%     c = 1;
%     f = v * dudx;
%     s = 0;
% end
% 
% % Define the initial condition for Couette flow
% function u0 = pdeIC1(x, v, H, u0)
%     u0 = zeros(size(x));   % Initialize all velocities to zero
%     u0(1) = u0;            % Set the velocity at the first node to u0
% end
% 
% % Define the boundary conditions for Couette flow
% function [pl, ql, pr, qr] = pdeBC1(xL, uL, xR, uR, t, v, H, u0)
%     pl = uL - u0;   % Left boundary condition: moving plate
%     ql = 0;         % Left boundary condition
%     pr = uR;        % Right boundary condition: stationary plate
%     qr = 0;         % Right boundary condition
% end

clc; clear all; close all;

%---Input Parameters
H = 4e-2; %Spacing between the plates
v = 2.17e-4; %kinematic viscosity
u0 = 40; % velocity of the moving plate
uR = 0;   % velocity of the stationary plate

% Grid
xL = 0;             % Left boundary
xR = H;             % Right boundary
xm = linspace(xL, xR, 6);   % Spatial grid
t0 = 0;             % Initial time
tend = 2;           % Final time
t = linspace(t0, tend, 200); % Time grid
xmt = xm';

% Solve the PDE using pdepe
tic
sol = pdepe(0, @pdefun1, @pdeIC1, @pdeBC1, xm, t, [], v, H, u0);
time_pdepe = toc;
fprintf('Comp Time: %e \n\n',time_pdepe);

%-----Post processing---------------
figure(1); % plot(x-y plot)
 plot(xmt,sol(10,:),'-o',xmt,sol(30,:),'-o',xmt,sol(60,:),'-o',xmt,sol(200,:),'-o');
 ylim([0, 40]); 
 legend('t=0.1','t=0.3','t=0.6','t=2');
 xlabel('Velocity (m/s)'); ylabel('Height (mm)');
 title('Plot using MATLAB pdepe function');

% Define the PDE for Couette flow
function [c, f, s] = pdefun1(x, t, u, dudx, v, H, u0)
    c = 1;
    f = v*dudx;
    s = 0;
end

% Define the initial condition for Couette flow
function u0 = pdeIC1(x, v, H, u0)
    u0 = zeros(size(x));   % Initialize all velocities to zero
    u0(1) = u0;            % Set the velocity at the first node to u0
end

% Define the boundary conditions for Couette flow
function [pl, ql, pr, qr] = pdeBC1(xL, uL, xR, uR, t, v, H, u0)
    pl = uL -u0;    % Left boundary condition: moving plate
    ql = 0;     % Left boundary condition
    pr = uR;    % Right boundary condition: stationary plate
    qr = 0;     % Right boundary condition
end

