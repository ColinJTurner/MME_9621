%---Input data-------------------------
L=5.5; h=0.01; Vavg=0.3; D=2e-5;Da=0.45; Cin=150e-6;Cmm=1e-5;
beta=L/h; hD=h/D;
%---Generating x and t mesh -----------
yL=0; yR=h; % here x is y
xmesh=linspace(yL/h,yR/h,101);
t0=0; tend=L; % here t is x
tmesh=linspace(t0/L,tend/L,101);
sol=pdepe(0,@pdefunBioreactor,@pdeIC_Bioreactor,@pdeBC_Bioreactor,xmesh,tmesh,[],beta,Vavg,hD,Da,Cin,Cmm);
sol=sol*1e6; % original magitudes of sol are very small
%-----Post processing---------------
figure(1); % plot(x-y plot)
 plot(tmesh,sol(:,1)); hold on;
 plot(tmesh,sol(:,10)); xlabel('x'); ylabel('C'); hold off;
figure(2); % contour plot (x-y-z plot) % contour levels are set automatically
 t=tmesh;x=xmesh;
 [mx,mt]=meshgrid(x,t);
 [CC,hh]=contour(mt,mx,sol,'LineWidth',2);
 xlabel('x'); ylabel('y'); set(hh,'ShowText','on');
figure(3); % contour plot (x-y-z plot) % contour levels are set by Clevel supplied by user
 Clevel=[0 5 10 20 30 40 50 60 80 100 110 130 140];
 [CC,hh]=contour(mt,mx,sol,Clevel,'LineWidth',2);
 xlabel('x'); ylabel('y');
 clabel(CC,hh,'manual'); % contour level texts are shown manually.
figure(4); % 3 D plot
 plot3(x,t,sol);
