function [xout iter_number]=Jacobi(A,b,max_iter,TOL)
    N=length(b);
    d=diag(A); % takes the main diagonal elements only
    x=zeros(N,1); % initial guess vector
        for k=1:max_iter
            x=(b-(A-diag(d))*x)./d; % Note: A=L+D+U, so L+U=A-D
            maxresidual=norm(b-A*x,inf);
            if maxresidual<TOL
                break;
            end
        end
iter_number=k;
xout=x;