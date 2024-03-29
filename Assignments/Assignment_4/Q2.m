clc; close all; clear all;
% Given parameters
k1 = 18000; k2 = 20000; k3 = 20000;
l1 = 1.0; l2 = 1.5;
M = 1000; m1 = 100; m2 = 200;
r = 0.9; J = M * r^2;
a = k1 + k2;
b = k1 + k2;
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
disp(frequenciesE);
disp('Natural Frequencies for Q-R Method:');
disp(frequenciesQ);

disp(['Computation Time for Eig Method: ', num2str(tE), ' seconds']);
disp(['Computation Time for Q-R Function: ', num2str(tQ), ' seconds']);

disp(['Accuracy (Relative Error): ', num2str(accuracy)]);