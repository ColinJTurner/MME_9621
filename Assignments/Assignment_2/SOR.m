function [xout iter_number]=SOR(A,b,relax,max_iter,TOL)
    N=length(b); d=diag(diag(A));
    x=zeros(N,1); % initial guess vector
    U=triu(A,1); % above the main diagonal
    L=tril(A,-1); % below the main diagonal
    for k=1:max_iter
        bL=relax*(b-U*x)+(1-relax)*d*x;
        for j=1:N
            x(j)=(bL(j)-relax*L(j,:)*x)./d(j,j);
        end
        maxresidual=norm(b-A*x,inf);
        if maxresidual<TOL
            break;
        end
    end
iter_number=k;
xout=x;