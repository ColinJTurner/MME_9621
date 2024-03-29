%% Q1
clc; close all; clear all;

for j =1:3
    n = 1e-11; % Step size for Simpson's 1/3 rule
    a = [0 1e-7 1e-8]; b = 10.5e-6; % Bounds
    tic;
    Qt(j) = quadl(@f, a(j), b);
    tQ = toc;
    tic;
    QS(j) = simpsons_rule(@f, a(j), b, n);
    tQS = toc;
    accuracy = norm(QS(j) - Qt(1)) / norm(Qt(1));
    formatspec = 'when bound "a" is: %5.1d\n';
    fprintf(formatspec,a(j));
    fprintf('1. Energy flux absorbed by the first plate using quadl: %5.5e W/m^2\n', Qt(j));
    fprintf("The computation time for quadl is: %5.3e s\n", tQ);
    fprintf("2. Energy flux absorbed by the first plate using quadl: %5.5e W/m^2\n", QS(j));
    fprintf("The computation time for Simpsons 1/3 rule is: %5.3e s\n", tQS);
    fprintf("The accuracy compared to quadl bound 0: %5.3e s\n\n", accuracy);
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
epsilon = (lambda <= 10.5e-6) .* 0.850 .* (1 - lambda./ 10.5e-6) + (lambda > 10.5) .* 0;
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
%% Q2
clc; close all; clear all;
% Given parameters
k1 = 18000; k2 = 20000; k3 = 20000;
l1 = 1.0; l2 = 1.5;
M = 1000; m1 = 100; m2 = 200;
r = 0.9; J = M * r^2;
a = k1 + k2;
b = k1 + k3;
c = k2 + k3;
d = k2*l1^2 + k3*l2^2;
x = [1,1,1,1];
max_iter = 100;
TOL = 10e-7;
A = [a/m1, 0, -k2/m1, k2*l1/m1;
    0, b/m2, -k3/m2, -k3*l2/m2;
    -k2/M, -k3/M, c/M, (k3*l2-k2*l1)/M;
    k2*l1/J, -k3*l2/J, (k3*l2-k2*l1)/J, d/J];

tic;
lamda_old = diag(A);
[n, m] = size(A);
Qbar = eye(n);
for k = 1:max_iter
    [Q, R] = qr(A);
    A = R * Q;
    Qbar = Qbar * Q;
    errornorm = norm(lamda_old - diag(A), inf);
    if errornorm < TOL
        break;
    end
    lamda_old = diag(A);
end
lamdaQR = diag(A);
iter_number = k;
xout = Qbar;
tQ = toc;


% Solve eigenvalue problem using eig function
tic;
[eigenvectors, lambda] = eig(A);
tE = toc;

% Extract natural frequencies
frequenciesE = sqrt(diag(lambda));
frequenciesQ = sqrt(lamdaQR);

% Calculate accuracy (relative error)
accuracy = norm(frequenciesQ - frequenciesE) / norm(frequenciesE);

disp('Natural Frequencies for Eig Function:');
disp(num2str(frequenciesE));
disp('Natural Frequencies for Q-R Method:');
disp(num2str(frequenciesQ));

disp(['Computation Time for Eig Method: ', num2str(tE), ' seconds']);
disp(['Computation Time for Q-R Function: ', num2str(tQ), ' seconds']);

disp(['Accuracy (Relative Error): ', num2str(accuracy)]);
