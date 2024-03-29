function [xout iter_number]=GaussSeidel(A,b,max_iter,TOL)
    N=length(b);
    d=diag(diag(A)); % takes the main diagonal elements only
    x=zeros(N,1); % initial guess vector
    U=triu(A,1); % above the main diagonal
    L=tril(A,-1); % below the main diagonal
    for k=1:max_iter
        bL=b-U*x;
        for j=1:N
            x(j)=(bL(j)-L(j,:)*x)./d(j,j);
        end
        maxresidual=norm(b-A*x,inf);
        if maxresidual<TOL
            break;
        end
    end
iter_number=k;
xout=x;