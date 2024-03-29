function [xout]=tridiagonal(a,b,c,d)
N=length(b);
%---forward elimination------------
for k=2:N
multiplier=a(k)/b(k-1);
b(k)=b(k)-multiplier*c(k-1);
d(k)=d(k)-multiplier*d(k-1);
end
%-----back substitution----------
x(N)=d(N)/b(N);
for k=N-1:-1:1
x(k)=(d(k)-c(k)*x(k+1))/b(k);
end
xout=x.';
