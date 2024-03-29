%% 
clc;clear all;format long e;
n=5;    t=35;

F=@(x)(tand(x).*tand(x-t)+tand(x).^2-2.*n^2.*sind(x).^2.*tand(x-t).^2);
x=[35:0.5:90];
for jj=1:length(x);
    y(jj)=F(x(jj));
end

plot(x,y);  grid on;
ylim([-30 30])

%% 

clc;clear all; format long e;
n=5; t=35;
F=@(x)(tand(x).*tand(x-t)+tand(x).^2-2.*n^2.*sind(x).^2.*tand(x-t).^2);
x=[t:0.5:90];
plot(x,F(x));
grid on;
ylim([-30 30]);
xlabel('x'); ylabel('y');
options=optimset('display','iter','tolX',1e-10); 
[xc1]=fzero(F,[t 60],options)
[xc2]=fzero(F,[60 85],options)
