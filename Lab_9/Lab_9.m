clc;close all; clear all;
I = dblquad(@f, 0, pi/2, 0, pi/2, 1e-7,'quadl');

fid = fopen('Lab9_Integration_result_Turner.txt','w');
fprintf(fid,'surface area of the tumor (cm2) = %f',I);
fclose(fid);


function [S]=f(theta,phi)
a = 9.5/2; b = 8/2; c = 4.2/2;
del=1-c^2/a^2;ep=1-c^2/b^2;
p=del.*sin(phi).^2+ep.*cos(phi).^2;
S=8.*a.*b.*sin(theta).*sqrt(1-p.*sin(theta).^2);
end

