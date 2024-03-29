clc; close all; clear all;

for j = 1:3
    a = [0 1e-7 1e-8];  % Corrected bounds
    b = 10.5e-6; % Upper bound, remains constant for all cases
    n = 1e-11; % Step size for Simpson's 1/3 rule (scalar)

    tic;
    Qt(j) = quadl(@f, a(j), b);
    tQ = toc;

    tic;
    QS(j) = simpsons_rule(@f, a(j), b, n);
    tQS = toc;

    accuracy = norm(QS(j) - Qt(j)) / norm(Qt(j));  % Corrected accuracy calculation

    fprintf('When bound "a" is: %5.1d\n', a(j));
    fprintf('1. Energy flux absorbed by the first plate using quadl: %5.5e W/m^2\n', Qt(j));
    fprintf("   Computation time for quadl is: %5.3e s\n", tQ);
    fprintf('2. Energy flux absorbed by the first plate using Simpsons 1/3 rule: %5.5e W/m^2\n', QS(j));
    fprintf("   Computation time for Simpsons 1/3 rule is: %5.3e s\n", tQS);
    fprintf("   Accuracy compared to quadl: %5.3e\n\n", accuracy);
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
epsilon = (lambda <= 10.5e-6) .* 0.850 .* (1 - lambda./ 10.5e-6) + (lambda > 10.5e-6) .* 0; % Corrected condition
rho = (lambda <= 4.5e-6) * 0.35 + (lambda > 4.5e-6) * 0.82;
E = 2 .* pi .* C1 ./ ((lambda).^5.*(exp(C2./(lambda.*T))-1));
f = rho.*epsilon.^2.*E;
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
