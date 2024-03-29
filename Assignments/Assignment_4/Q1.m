clc; close all; clear all;

for j =1:4

n = 1e-7; % Step size for Simpson's 1/3 rule
a = [0 1e-6 5e-7 1e-9]; b = 4.5e-6; c = 10.5e-6; % Bounds
tic;
Q_1 = 0.35*quadl(@f, a(j), b); 
Q_2 = 0.82*quadl(@f, b, c); 
Qt = Q_1 + Q_2;
tQ = toc;
tic;
Q_1S = 0.35*simpsons_rule(@f, a(j), b, n);
Q_2S = 0.82*simpsons_rule(@f, b, c, n);
QtS = Q_1S + Q_2S;
tQS = toc;

formatspec = 'when bound "a" is: %5.1d\n';
fprintf(formatspec,a(j));
fprintf('1. Energy flux absorbed by the first plate using quadl: %5.5e W/m^2\n', Qt);
fprintf("From 0 < lambda <= 4.5: %5.5e W/m^2\n", Q_1);
fprintf("From 4.5 < lambda <= 10.5: %5.5e W/m^2\n", Q_2);
fprintf("The computation time for quadl is: %5.3e s\n", tQ);
fprintf("2. Energy flux absorbed by the first plate using quadl: %5.5e W/m^2\n", QtS);
fprintf("From 0 < lambda <= 4.5: %5.5e W/m^2\n", Q_1S);
fprintf("From 4.5 < lambda <= 10.5: %5.5e W/m^2\n", Q_2S);
fprintf("The computation time for Simpsons 1/3 rule is: %5.3e s\n\n", tQS);
end

function f = f(lambda)
% Constants
c0 = 2.9979e8;  % Speed of light (m/s)
h = 6.626e-34;  % Planck's constant (Js)
kB = 1.3806e-23;  % Boltzmann constant (J/K)
T = 925;  % Temperature (K)

C1 = h * c0^2;
C2 = h * c0 / kB;

% Functions
epsilon = 0.850 * (1 - lambda / 10.5e-6);
E = 2 .* pi .* C1 ./ ((lambda).^5.*(exp(C2./(lambda.*T))-1));
f = epsilon.^2.*E;
end

function Q = simpsons_rule(func, a, b, h)

    % Number of intervals
    n = (b - a) / h;

    % Initialize sum
    sum = func(a) + func(b);

    % Odd terms
    odd_sum = 0;
    for i = 1:2:n-1
        x = a + i * h;
        odd_sum = odd_sum + func(x);
    end

    % Even terms
    even_sum = 0;
    for i = 2:2:n-2
        x = a + i * h;
        even_sum = even_sum + func(x);
    end

    % Final result
    Q = h / 3 * (sum + 4 * odd_sum + 2 * even_sum);
end