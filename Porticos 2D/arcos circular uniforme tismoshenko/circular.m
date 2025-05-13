

clc
clear
% Definicion de simbolos
syms f L x qa qb EI ryy rxx mm phi

% Parámetros de entrada
qa = 10;% 5;        % Carga distribuida (N/m)
qb = 13; %15;         % Carga distribuida (N/m)
bv = 25/1000;   % Ancho de la viga (m)
hv = 50/1000;   % Altura de la viga (m)
I = bv * hv^3 / 12; % Inercia de la viga (m^4)
E = 210000 * 1000; % Módulo de elasticidad del material (Pa)
v = 0.3;        % Coeficiente de Poisson
G = E / (2 * (1 + v)); % Módulo de rigidez en corte (Pa)
EI = E * I;     % Producto de módulo de elasticidad y momento de inercia (N·m²)
r=4;
L=2*r;
Ac=bv*hv*G*5/6;
DX=0;
DY=0;
GM=0;
%% Función circulo

%y= coe(1)*(x+xa(1))^2+coe(2)*(x+xa(1))+coe(3)-ya(1); %% funcion de viga
y=r*sin(phi);
x=r*cos(phi)+r;
dydx=diff(x,phi);
%% Función de momento para una carga trapezoidal
% Momento flector en la viga debido a la carga trapezoidal
%q=(qb-qa)/L*x+qa
%V=int(q,x)
%M=int(V,x)
a=0;
b=pi;
h=0;%yn(3)-yn(1);
V = qa*x - (x^2*(qa - qb))/(2*L);

M = (-((qa - qb) * x^3) / (6 * L) + (qa * x^2) / 2);
%% funcion para calcular  fuerzas nodales equivalentes en cordenaadas positivas

R=Arc_circular(r,y,x,dydx,EI,Ac,M,V,a,b,L,h,qa,qb,DX,DY,GM); %apoyo derecho

%% funcion para calcular  fuerzas nodales equivalentes en cordenaadas negativas
%y1=coe(1)*(x+xa(3))^2+coe(2)*(x+xa(3))+coe(3)-ya(3); %% funcion de viga
y1=(x*(xa^2*yb - xb^2*ya - xa^2*yc - xc^2*ya + xb^2*yc + xc^2*yb - x*xa*yb + x*xb*ya + x*xa*yc - x*xc*ya - x*xb*yc + x*xc*yb - 2*xa*xc*yb + 2*xb*xc*ya + 2*xa*xc*yc - 2*xb*xc*yc))/((xa - xb)*(xa - xc)*(xb - xc));
dy1dx=(xa^2*yb - xb^2*ya - xa^2*yc - xc^2*ya + xb^2*yc + xc^2*yb - 2*x*xa*yb + 2*x*xb*ya + 2*x*xa*yc - 2*x*xc*ya - 2*x*xb*yc + 2*x*xc*yb - 2*xa*xc*yb + 2*xb*xc*ya + 2*xa*xc*yc - 2*xb*xc*yc)/((xa - xb)*(xa - xc)*(xb - xc));
 
a=-L;
b=0;
qa=10;
qb=13;
M = (-((qa - qb) * x^3) / (6 * -L) + (qa * x^2) / 2);
R3=calculos_parabolico(y1,EI,M,a,b,-L,h,qa,qb,DX,DY,GM);
R3t=calculos_parabolico_tismoshenko(y1,dy1dx,EI,Ac,M,V,a,b,-L,h,qa,qb,DX,DY,GM);
%% calaculo de matriz de rigidez
a=0;
b=L;
qa=0;
qb=0;
M =0;
DX=1;
DY=0;
GM=0;

K1=calculos_parabolico(y,EI,M,a,b,L,h,qa,qb,DX,DY,GM);

K1t=calculos_parabolico_tismoshenko(y,dydx,EI,Ac,M,V,a,b,L,h,qa,qb,DX,DY,GM);

DX=0;
DY=1;
GM=0;

K2=calculos_parabolico(y,EI,M,a,b,L,h,qa,qb,DX,DY,GM);
K2t=calculos_parabolico_tismoshenko(y,dydx,EI,Ac,M,V,a,b,L,h,qa,qb,DX,DY,GM);

DX=0;
DY=0;
GM=1;
K3=calculos_parabolico(y,EI,M,a,b,L,h,qa,qb,DX,DY,GM);
K3t=calculos_parabolico_tismoshenko(y,dydx,EI,Ac,M,V,a,b,L,h,qa,qb,DX,DY,GM);

DX=1;
DY=0;
GM=0;
a=-L;
b=0;

K4=calculos_parabolico(y1,EI,M,a,b,-L,h,qa,qb,DX,DY,GM);
K4t=calculos_parabolico_tismoshenko(y1,dy1dx,EI,Ac,M,V,a,b,-L,h,qa,qb,DX,DY,GM);

DX=0;
DY=1;
GM=0;

K5=calculos_parabolico(y1,EI,M,a,b,-L,h,qa,qb,DX,DY,GM);
K5t=calculos_parabolico_tismoshenko(y1,dy1dx,EI,Ac,M,V,a,b,-L,h,qa,qb,DX,DY,GM);

DX=0;
DY=0;
GM=1;

K6=calculos_parabolico(y1,EI,M,a,b,-L,h,qa,qb,DX,DY,GM);
K6t=calculos_parabolico_tismoshenko(y1,dy1dx,EI,Ac,M,V,a,b,-L,h,qa,qb,DX,DY,GM);

