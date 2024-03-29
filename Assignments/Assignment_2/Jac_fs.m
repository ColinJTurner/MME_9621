function [fJac]=Jac_fs(Xn,DFs,Xs)
fJac=subs(DFs,Xs,Xn);
fJac=double(fJac);
