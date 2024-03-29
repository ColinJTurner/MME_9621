function [xout,k]=multivariateNewton_fs(F,a, b, c, d, f, e,J,DFs,Xs,x0,TOL,max_iter)
k=0;
enorm=10; % dummy to start the iteration process
x=x0;xold=x0;
while enorm>TOL
    dx=J(x,DFs,Xs)\F(x,a, b, c, d, f, e); % J-Jacobian (computed symbolically), % F-function vector
    x=x-dx;
    k=k+1;
    enorm=norm(x-xold,inf);
    xold=x;
    %fprintf(1,'%d %15.10f %15.10f %15.10f\n',k,enorm,x(1),x(2));
    if k>=max_iter
        break;
    end
end
xout=x;