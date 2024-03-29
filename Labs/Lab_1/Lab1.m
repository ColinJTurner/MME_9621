format long e;
clc; % clear the command window
clear all;
%definitions
d1 = 38e-3;
v1 = 0.3;    d2 = 70e-3;
v2 = v1;     F = 450;
E1 = 206e+9; L = 50e-3;
E2 = E1;     z = 0.025e-3;
%optimization
F_2_L = 2*F/pi/L;
%b & Pmax EQN
b= sqrt((F_2_L)*(((1-v1^2)/E1+(1-v2^2)/E2)/(1/d1+1/d2)));
Pm = F_2_L/b;
%optimization
z_b = z/b;   z_b_1 = sqrt(1+z_b^2);

shearx = -2*v2*Pm*(z_b_1-z_b);
fprintf ('shear x = %e \v', shearx);
sheary = -Pm*(((2-z_b_1^(-2))*z_b_1)-2*z_b);
fprintf ('shear y = %e \n', sheary);
shearz = -Pm/z_b_1;
fprintf ('shear z = %e \v', shearz);
tensileyz = 0.5*(sheary-shearz);
fprintf ('tensileyz = %e \n', tensileyz);
