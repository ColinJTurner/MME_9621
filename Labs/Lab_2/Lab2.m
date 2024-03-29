format long e;
clc; % clear the command window
clear all;
%definitions
h = 7;     a = 1;
r = 2;     b = 10;
H = 10;    n = 12;
           s = 0;
Hr = H-r;  b1 = 1/b;
            
%Problem 1
if (h<r)
    V = 1/3*pi*h^2*3*r-h;
elseif (h<Hr)
    V = 2/3*pi*r^3+pi*r^2*(h-r);
else
    V = 4/3*pi*r^3+pi*r^2*(H-2*r)-1/3*pi*(H-h)^2*(3*r-H+h);
end;

fprintf ('V = %e \n', V);

% Problem 2
a = 1;
b = 10;
n = 12;
s = 0;
b1 = 1/b;
y1 = [];
y2 = [];
    % linsapce(x1,x2,n)
    %n is generated points
    s=  linspace(0,360,600);
    J = a+b1.*tanh(b.*sind(n.*s));
    y1 = J.*cosd(s);
    y2 = J.*sind(s);
   
    plot (y1,y2);