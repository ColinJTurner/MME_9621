clc;close all; clear all;
R=0.02; D=2e-5; Rho=1;Smax=15e-9; K=25e-10; Ob=5.234e-8;
%-------------------------------
S=Smax*R^2/(Rho*D*Ob); KOb=K/Ob;
%------------Boundary -------------------
chiL=0.0005; chiR=1.0; % left(a) and right(b) boundary
 % since at chi=0, there is a Singularity, we change chiL to a small number near zero.
%---Initial guess---- %-------assume y=psi------
N=5; % number of grid point (x-point)
x=linspace(chiL,chiR,N); % initial mesh/grid
yinit=[0 0]; % y=0, dy/dx=0 for initial guess, (x=chi)
solinit=bvpinit(x,yinit); % initial solution compatible for ‘bvp4c’
%---call the solver-------------
sol=bvp4c(@TumorODE,@TumorBC,solinit,[],S,KOb); % S and KOb are parameters
%=== Post processing ============
chifinal=sol.x;
y_sol=sol.y(1,:); %-Row 1 of sol.y represents y
dy_dx_sol=sol.y(2,:); %-Row 2 of sol.y represents dy/dx, (x=chi)
figure(1);
plot(chifinal,y_sol); %-the plot is not smooth
xlabel('\chi'); ylabel('\psi');
%---take more points----
N=500; % number of x-points to make the plot smoother
chi=linspace(chiL,chiR,N); % for interpolating the result to get
y=deval(sol,chi); % smooth graph
y=y'; % transpose
figure(2);
plot(chi,y(:,1),chi,y(:,2)); % now, column 1: y(:,1) represents y
xlabel('\chi'); ylabel('y'); % column 2: y(:,2) represents dy/dx, (x=chi)
legend('\psi','d\psi/d\chi');