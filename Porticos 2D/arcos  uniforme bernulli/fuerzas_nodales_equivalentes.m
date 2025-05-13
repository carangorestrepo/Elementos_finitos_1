function R=fuerzas_nodales_equivalentes(qa,qb,xa,xb,xc,ya,yb,yc,EI)
% Definicion de simbolos
syms  x 

%% Funci�n de la par�bola

%y= coe(1)*(x+xa(1))^2+coe(2)*(x+xa(1))+coe(3)-ya(1); %% funcion de viga
y=-(x*(xa^2*yb + xb^2*ya - xa^2*yc - xc^2*ya - xb^2*yc + xc^2*yb + x*xa*yb - x*xb*ya - x*xa*yc + x*xc*ya + x*xb*yc - x*xc*yb - 2*xa*xb*ya + 2*xa*xc*ya + 2*xa*xb*yc - 2*xa*xc*yb))/((xa - xb)*(xa - xc)*(xb - xc));

%% Funci�n de momento para una carga trapezoidal

L=xc-xa;
DX=0;
DY=0;
GM=0;
a=0;
b=xc-xa;
h=yc-ya;
%V = qa*x - (x^2*(qa - qb))/(2*L);
M = (-((qa - qb) * x^3) / (6 * L) + (qa * x^2) / 2);
%% funcion para calcular  fuerzas nodales equivalentes en cordenaadas positivas
R=calculos_parabolico(y,EI,M,a,b,L,h,qa,qb,DX,DY,GM);