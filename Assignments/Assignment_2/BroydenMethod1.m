function [xout,k]=BroydenMethod1(F,a, b, c, d, f, e,x0,x1,A,TOL,max_iter)
    for k=1:max_iter
        deltax=x1-x0; deltaF=F(x1,a, b, c, d, f, e)-F(x0,a, b, c, d, f, e);
        A=A+(deltaF-A*deltax)*deltax'/(deltax'*deltax);
        dx=A\F(x1,a, b, c, d, f, e);
        x=x1-dx;
        enorm=norm(x-x1,inf); x0=x1;x1=x;
        %fprintf(1,'%d %15.10f %15.10f %15.10f\n',k,enorm,x(1),x(2));
        if enorm <TOL || k >= max_iter
            break;
        end
    end
xout=x;
