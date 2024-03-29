function [xout,k]=newton(x,f,fprime,tol,max_iteration)
k=0;
xnew=x;
xold=x+5*tol; % dummy, used only to start the while loop
fprintf(1,'%d %15.10f \n',k,xnew);
while abs(xnew-xold)>tol
 xold=xnew;
 xnew=xold-f(xold)/fprime(xold);
 k=k+1;
 fprintf(1,'%d %15.10f \n',k,xnew);
 if k>=max_iteration
 break;
 end
end
xout=xnew;