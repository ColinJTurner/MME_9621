clc;close all; clear all;
s=[40.0,37.5,43.0,52.0,60.0,55.0];
e=[0.02,0.05,0.10,0.15,0.20,0.25];

I = trapz(s,e);
fprintf('Modulus of Toughness = %f',I);

