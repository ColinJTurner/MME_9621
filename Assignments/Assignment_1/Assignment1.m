%% Q2 report # of correct digits and comment
clc;clear all;format long e; % done
F1=@(x)((1-sec(10^-x))/tan(10^-x)^2);
F2=@(x)(-1)/(1+sec(10^-x));
fprintf("No Error Corrections\t");
fprintf("\t\tError Correction\n");
x1=1;
while x1<11
    formatSpec='10^-%d is: %1.20f\t';
    fprintf(formatSpec,x1,F1(x1));
    formatSpec='10^-%d is: %1.20f\n';
    fprintf(formatSpec,x1,F2(x1));
    x1=x1+1;
end
%% Q3
clc;clear all;format long e;
L=25;d=0.1;ep=0.8;sig=5.67*10^-8;Q=18405;h=10;Ta=298;Tsurr=298;
a=pi*d*L*h;b=pi*d*L*ep*sig;
c=a*Ta;d=b*Tsurr^4;
F=@(x)(a*x-c+b*x^4-d-Q);

r1=bisectionMethod(F,400,501,1e-10);
formatSpec='Using Bisection Method Root 1 found at:%e\n';
fprintf(formatSpec,r1);

options=optimset('display','iter','tolX',1e-10); 
[xc1]=fzero(F,[400 501],options)
formatSpec='Using fZero Method Root 1 found at:%e\n';
fprintf(formatSpec,xc1);

x=[-9000:0.1:9000];
for jj=1:length(x);
    y(jj)=F(x(jj));
end
plot(x,y);  grid on;
%ylim([-30000 30000])
xlabel('x'); ylabel('y');
%% Q4
clc;clear all;format long e;
L=3;E=70*10^9;I=52.9*10^-6;w=15000;
z=w./(120.*E.*I);
a=3.*z.*L.^2;b=7.*z.*L;c=5.*z;d=z./L;
F=@(x)(a.*x.^2-b.*x.^3+c.*x.^4-d.*x.^5);
dF=@(x)(2.*a.*x-3.*b.*x.^2+4.*c.*x.^3-5.*d.*x.^4);
ddF=@(x)(2.*a-6.*b.*x+12.*c.*x.^2-20.*d.*x.^3);
r1=newton(1,dF,ddF,1e-10,100);
formatSpec='Using Newton Method Root 1 found at:%1.10e\nThe Deflection at this root is:%fm';
fprintf(formatSpec,r1,F(r1));

options=optimset('display','iter','tolX',1e-10); 
[xc]=fzero(dF,1,options)
formatSpec='Using fZero Method Root 1 found at:%1.10e\nThe Deflection at this root is:%fm';
fprintf(formatSpec,xc,F(xc));

x=[0:0.001:3.5];
for jj=1:length(x);
    y(jj)=F(x(jj));
end

plot(x,y);  grid on;
%ylim([-0.5 0.5])
xlabel('Distance (m)'); ylabel('Deflection(m)');
title('Deflection Curve');
%% Q5
clc;clear all;format long e;
L=0.02;k=0.5;q=10^6;Ta=100;Tb=200;dx=0.004;A=1;
a=k*A/dx;b=2*a;c=q*A*dx;d=a+b;
e=c+b*Ta;f=c+b*Tb;
% / Method
A=[d -a 0 0 0; 
   -a b -a 0 0; 
   0 -a b -a 0;
   0 0 -a b -a;
   0 0 0 -a d;
   ]; 
B=[e;c;c;c;f];
tic;
x1=A\B;
formatSpec='Using Backslash Method:\nT1=%3.0f\nT2=%3.0f\nT3=%3.0f\nT4=%3.0f\nT5=%3.0f\nComp Time:%fs\n';
fprintf(formatSpec,x1(1),x1(2),x1(3),x1(4),x1(5),toc);
% Inversion Method
tic;
xc=inv(A)*B;
formatSpec='Using Inversion Method:\nT1=%3.0f\nT2=%3.0f\nT3=%3.0f\nT4=%3.0f\nT5=%3.0f\nComp Time:%fs\n';
fprintf(formatSpec,xc(1),xc(2),xc(3),xc(4),xc(5),toc);
% Tridiagonal Method
tic;
a1=[0 -a -a -a -a];
b1=[d b b b d];
c1=[-a -a -a -a 0];
d1=[e;c;c;c;f];
[xout]=tridiagonal(a1,b1,c1,d1);
formatSpec='Using Tridiagonal Method:\nT1=%3.0f\nT2=%3.0f\nT3=%3.0f\nT4=%3.0f\nT5=%3.0f\nComp Time:%fs\n';
fprintf(formatSpec,xout(1),xout(2),xout(3),xout(4),xout(5),toc);

F=@(x)(Ta+((Tb-Ta)/L+q*L/2/k)*x-q*x^2/2/k);

x=[0:0.0001:0.02];
for jj=1:length(x); 
    y(jj)=F(x(jj));
end
plot(x,y);  grid on;

xlabel('Distance (m)'); ylabel('Temperature (C)');
title('Temperature Distribution');
hold on
x=0.002:0.004:0.018;
stem(x,x1);
ylim([100 300])
hold off