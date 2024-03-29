clc;close all; clear all;

Q = quadl(@f,0,10.5e-6);

disp(['Energy flux absorbed by the first plate: ', num2str(Q) , ' W/m^2']);
% Define wavelength range
lambda_range = linspace(0, 11e-6, 1000); 

% Calculate energy absorbed for each wavelength
energy_absorbed = arrayfun(@f, lambda_range);

% Plot the energy absorbed
plot(lambda_range, energy_absorbed);
xlabel('\lambda (um)');
ylabel('Energy absorbed (W/m^2)');
title('Energy absorbed by the first plate');
grid on;

function f = f(lambda)
% Constants
c0 = 2.9979e8;  % Speed of light (m/s)
h = 6.626e-34;  % Planck's constant (Js)
kB = 1.3806e-23;  % Boltzmann constant (J/K)
T = 925;  % Temperature (K)
% Function for E(𝜆, T)
C1 = h .* c0.^2;
C2 = h .* c0 ./ kB;

if lambda <= 10.5e-6
    epsilon = 0.850 .* (1 - lambda.*10.^-6 ./ 10.5);
else
    epsilon =0;
end
if lambda <= 4.5e-6
    rho = 0.35;
else
    rho = 0.82;
end
E = 2 .* pi .* C1 ./ (lambda.^5.*10.^-6.*(exp(C2./lambda.*10.^-6.*T)-1));
f = rho.*epsilon.^2.*E;
end


