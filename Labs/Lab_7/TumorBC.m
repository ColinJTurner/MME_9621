function res=TumorBC(ya,yb,S,KOb)
res=[ya(2)-0; % 2 is for derivative(Neumann) BC: dy/dx(x=a)=0
 yb(1)-1]; % 1 is for Dirichlet BC: y(x=b)=1
% Note: since TumorODE fun have parameters (S, KOb), TumorBC also must have these
 % parameters, though there is No use of them in BC.
