function [pL,qL,pR,qR]=pdeBC_Bioreactor(xL,uL,xR,uR,t,beta,Vavg,hD,Da,Cin,Cmm)
pL=-Da*Cin*uL/(uL+Cmm);
qL=1/beta;
pR=0;
qR=1;
