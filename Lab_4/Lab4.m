clc;clear all;format long e;
[Aug1 Text]=xlsread('Datafile_Lab4.xlsx');
A=Aug1(:,1:14); 
b=Aug1(:,15); 
x=A\b;
x1=inv(A)*b;
fprintf("Using Backslash Method\t Using Inv Method\n");
for i = 1:1:14
    fomatSpec='F%d = %12.5f\t F%d = %12.5f\n';
    fprintf(fomatSpec,i,x(i),i,x1(i));
end