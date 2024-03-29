function [c,f,s]=pdefunBioreactor(y,x,u,dudx,beta,Vavg,hD,Da,Cin,Cmm)
V=Vavg*(1-(2*y-1).^2)*3/2;
K=V*hD;
c=K;
f=dudx*beta;
s=0;