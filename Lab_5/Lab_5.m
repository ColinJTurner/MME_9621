%Example - 1a: fsolve function: for numerical solution
clc; clear all; format long e;
%---Anonymous function----------------
Fn=@(x,c1,c2) [x(1)-x(1).^2-c1.*x(1).*x(2)-c2.*x(1).*x(3);x(2)-x(2).^2-c2.*x(2).*x(1)-c1.*x(2).*x(3);x(3)-x(3).^2-c1.*x(3).*x(1)-c2.*x(3).*x(2)];
options=optimset('display','iter');
c1=0.3; c2=0.6;
x0=[0.5 0 0.5]; % initial guess vector
[xc]=fsolve(Fn,x0,options,c1,c2)
[xc]=fsolve(Fn,x0,[],c1,c2); % if no options is used
[xc fxc]=fsolve(Fn,x0,options,c1,c2)

