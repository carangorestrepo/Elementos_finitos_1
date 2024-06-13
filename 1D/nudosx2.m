%manera  1
syms x1 x2 x L u1 u2

N1=simplify((0-1)/(L)*(x-x2));
N2=simplify((1-0)/(L)*(x-x1));
%% super posicion
u=N1*u1+N2*u2

u=(u2-u1)/(x2-x1)*(x-x1)+u1;

u=collect(u,[u1,u2]);
disp('u=');pretty(u)



