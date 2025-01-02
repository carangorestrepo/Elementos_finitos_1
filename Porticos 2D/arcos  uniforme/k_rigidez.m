function K=k_rigidez(xa,xb,xc,ya,yb,yc,EI)
syms x
h=yc-ya;
L=xc-xa;
%% funcion para calcular  fuerzas nodales equivalentes en cordenaadas positivas
%y= coe(1)*(x+xa(1))^2+coe(2)*(x+xa(1))+coe(3)-ya(1); %% funcion de viga
y=-(x*(xa^2*yb + xb^2*ya - xa^2*yc - xc^2*ya - xb^2*yc + xc^2*yb + x*xa*yb - x*xb*ya - x*xa*yc + x*xc*ya + x*xb*yc - x*xc*yb - 2*xa*xb*ya + 2*xa*xc*ya + 2*xa*xb*yc - 2*xa*xc*yb))/((xa - xb)*(xa - xc)*(xb - xc));

%% funcion para calcular  fuerzas nodales equivalentes en cordenaadas negativas
%y1=coe(1)*(x+xa(3))^2+coe(2)*(x+xa(3))+coe(3)-ya(3); %% funcion de viga
y1=(x*(xa^2*yb - xb^2*ya - xa^2*yc - xc^2*ya + xb^2*yc + xc^2*yb - x*xa*yb + x*xb*ya + x*xa*yc - x*xc*ya - x*xb*yc + x*xc*yb - 2*xa*xc*yb + 2*xb*xc*ya + 2*xa*xc*yc - 2*xb*xc*yc))/((xa - xb)*(xa - xc)*(xb - xc));
%% calaculo de matriz de rigidez
a=0;
b=L;
qa=0;
qb=0;
M =0;
DX=1;
DY=0;
GM=0;
%% grado de libertad Hozontal inicial
K1=calculos_parabolico(y,EI,M,a,b,L,h,qa,qb,DX,DY,GM);
DX=0;
DY=1;
GM=0;
%% grado de libertad Vertical inicial
K2=calculos_parabolico(y,EI,M,a,b,L,h,qa,qb,DX,DY,GM);
DX=0;
DY=0;
GM=1;
%% grado de libertad giro inicial
K3=calculos_parabolico(y,EI,M,a,b,L,h,qa,qb,DX,DY,GM);
DX=1;
DY=0;
GM=0;
a=-L;
b=0;
%% grado de libertad Hozontal final
K4=calculos_parabolico(y1,EI,M,a,b,-L,h,qa,qb,DX,DY,GM);
DX=0;
DY=1;
GM=0;
%% grado de libertad Vertical final
K5=calculos_parabolico(y1,EI,M,a,b,-L,h,qa,qb,DX,DY,GM);
DX=0;
DY=0;
GM=1;
%% grado de libertad giro final
K6=calculos_parabolico(y1,EI,M,a,b,-L,h,qa,qb,DX,DY,GM);
%% signos de matriz de rigidez
%sing=[1,1,-1,-1,-1,1;1,1,1,-1,-1,1;-1,1,1,1,-1,-1;-1,-1,1,1,1,-1;-1,-1,-1,1,1,-1;1,1,-1,-1,-1,1];
%% Armado de matriz de rigidez 
K=[-K1,K2,-K3,K4,-K5,K6];
%K=[abs(K1).*sing(:,1),abs(K2).*sing(:,2),abs(K3).*sing(:,3),abs(K4).*sing(:,4),abs(K5).*sing(:,5),abs(K6).*sing(:,6)];


  

