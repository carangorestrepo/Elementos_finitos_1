
clc
clear
syms x1 y1 z1 x2 y2 z2 x3 y3 z3 x y z  lcx lcy cx0 cy0 pm ex ey
%lcx = 0.35;   % Largo columna [m]
%lcy = 0.65;   % Ancho columna [m]

%cx0 = 1.1;   % Posición x columna [m]
%cy0 = 1;     % Posición y columna [m]

%Mx=200;
%My=150;
%P=-1000;
%ex=Mx/P;
%ey=My/P;
% pm=P/(lcx*lcy);
p1=pm*(1+6*abs(ex)/lcx+6*abs(ey)/lcy);
p2=pm*(1-6*abs(ex)/lcx+6*abs(ey)/lcy);
p3=pm*(1-6*abs(ex)/lcx-6*abs(ey)/lcy);
p4=pm*(1+6*abs(ex)/lcx-6*abs(ey)/lcy);

x1=cx0;
y1=cy0;
z1=p1;

x2=cx0+lcx;
y2=cy0;
z2=p2;

x3=cx0+lcx;
y3=cy0+lcy;
z3=p3;

A=[x1,y1,z1];
B=[x2,y2,z2];
C=[x3,y3,z3];
V=[x,y,z];
AB=B-A;
BC=C-B;
AV=V-A;

sol=det([AB',BC',AV']);

Z=simplify(solve(sol,z));

xx=linspace(cx0,cx0+lcx,30);
yy=linspace(cy0,cy0+lcy,30);

[X,Y] = meshgrid(xx,yy);

surf(X,Y,double(subs(Z,{x,y},{X,Y})));





