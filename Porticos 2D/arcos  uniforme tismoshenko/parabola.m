clc
clear
syms x1 x2 x3 y1 y2 y3 x
%calculo de ecuacion parabola
x1=0;y1=0;
x2=2.5;y2=1.5;
x3=4;y3=0.6;
%{
X=[x1^2,x1,1;
   x2^2,x2,1; 
   x3^2,x3,1]; 
Y=[y1;y2;y3];
coe=X^(-1)*Y;
yp=coe(1)*x^2+coe(2)*x+coe(2);%funcion de arco positivo
%}
hx =(x1^2*y2 - x2^2*y1 - x1^2*y3 + x3^2*y1 + x2^2*y3 - x3^2*y2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
ky =(x1^4*y2^2 - 2*x1^4*y2*y3 + x1^4*y3^2 + 4*x1^3*x2*y2*y3 - 4*x1^3*x2*y3^2 - 4*x1^3*x3*y2^2 + 4*x1^3*x3*y2*y3 - 2*x1^2*x2^2*y1*y2 - 2*x1^2*x2^2*y1*y3 - 2*x1^2*x2^2*y2*y3 + 6*x1^2*x2^2*y3^2 + 4*x1^2*x2*x3*y1*y2 + 4*x1^2*x2*x3*y1*y3 - 8*x1^2*x2*x3*y2*y3 - 2*x1^2*x3^2*y1*y2 - 2*x1^2*x3^2*y1*y3 + 6*x1^2*x3^2*y2^2 - 2*x1^2*x3^2*y2*y3 + 4*x1*x2^3*y1*y3 - 4*x1*x2^3*y3^2 + 4*x1*x2^2*x3*y1*y2 - 8*x1*x2^2*x3*y1*y3 + 4*x1*x2^2*x3*y2*y3 - 8*x1*x2*x3^2*y1*y2 + 4*x1*x2*x3^2*y1*y3 + 4*x1*x2*x3^2*y2*y3 + 4*x1*x3^3*y1*y2 - 4*x1*x3^3*y2^2 + x2^4*y1^2 - 2*x2^4*y1*y3 + x2^4*y3^2 - 4*x2^3*x3*y1^2 + 4*x2^3*x3*y1*y3 + 6*x2^2*x3^2*y1^2 - 2*x2^2*x3^2*y1*y2 - 2*x2^2*x3^2*y1*y3 - 2*x2^2*x3^2*y2*y3 - 4*x2*x3^3*y1^2 + 4*x2*x3^3*y1*y2 + x3^4*y1^2 - 2*x3^4*y1*y2 + x3^4*y2^2)/(4*(x1 - x2)*(x1 - x3)*(x2 - x3)*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
p =((x1 - x2)*(x1 - x3)*(x2 - x3))/(4*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
 
y=(x-hx)^2/(-4*p)+ky;
focoy=ky-p;
focox=hx;

dydxp=diff(y,x,1);%derivada 
ds=sqrt(1+dydxp^2);% longitud d arco

%componentes
%cos=1 ./ sqrt(1 + dydx.^2)
Vu=1/ds;%% cortante unitario    % cos(theta)
%sin=dydx ./ sqrt(1 + dydx.^2); % sin(theta)
Pu=dydxp/ds;%% axial unitario

bf=25/1000;
d=50/1000;
E=210000000;
Ae=bf*d;
I=bf*d^3/12;
EI=E*I;
EA=E*Ae;
G=81018518.3;
Ac=Ae*5/6*G;
qa=5;
qb=15;
a=0;
b=x3;
L=x3;
h=y3;
M = (-((qa - qb) * x^3) / (6 * L) + (qa * x^2) / 2);
V = qa*x - (x^2*(qa - qb))/(2*L);
Vz=Vu*V;
Pz=Pu*V;

DX=0;
DY=0;
GM=0;

%% 1. DEFORMACIÓN VERTICAL (Qy = 1)
AMxx=area_cuadraturas(a,b,matlabFunction((x*x/EI+Vu/Ac)*ds,"Vars", x)); % Deformación vertical por carga unitaria vetical
AMxy=area_cuadraturas(a,b,matlabFunction((x*y/EI+Pu/Ac)*ds,"Vars", x)); % Deformación horizontal por carga unitaria vetical
AMxg=area_cuadraturas(a,b,matlabFunction((x/EI)*ds,"Vars", x));         % giro por carga unitaria vetical

%% 2. DEFORMACIÓN HORIZONTAL (Qx = 1)
%%NOta los diagramas de cortante y carga axial se invierten 
AMyx=area_cuadraturas(a,b,matlabFunction((y*x/EI+Pu/Ac)*ds,"Vars", x)); % Deformación vertical por carga unitaria horizonal
AMyy=area_cuadraturas(a,b,matlabFunction((y*y/EI+Vu/Ac)*ds,"Vars", x)); % Deformación horizontal por carga unitaria horizonal
AMyg=area_cuadraturas(a,b,matlabFunction((y/EI)*ds,"Vars", x));         % giro por carga por carga unitaria horizonal

%% 3. DEFORMACIÓN GIRO (QG = 1)
AMx =area_cuadraturas(a,b,matlabFunction((x/EI)*ds,"Vars", x));               % Deformación vertical por carga unitaria momento
AMy =area_cuadraturas(a,b,matlabFunction((y/EI)*ds,"Vars", x));               % Deformación horizontal por carga unitaria momento
AMg=area_cuadraturas(a,b,matlabFunction((1/EI)*ds,"Vars", x));                % giro por carga  por carga unitaria momento

%% 4. DEFORMACIÓN CARGA EXTERNA Y
AMMx=area_cuadraturas(a,b,matlabFunction((M*x/EI+ Vz/Ac)*ds,"Vars", x));     % Deformación vertical por carga unitaria vetical
AMMy=area_cuadraturas(a,b,matlabFunction((M*y/EI+ Pz/Ac)*ds,"Vars", x));     % Deformación horizontal por carga unitaria vetical
AMMg=area_cuadraturas(a,b,matlabFunction((M/EI)*ds,"Vars", x));              % giro por carga unitaria vetical

%% ==============================================================
%% 6. SISTEMA DE ECUACIONES (compatibilidad)
%% ==============================================================
% Matriz de coeficientes
ec=[AMxx ,AMxy ,AMxg;  % Suma de fuerzas en Y
    AMyx ,AMyy ,AMyg;  % Suma de fuerzas en X
    AMx  ,AMy  ,AMg];  % Suma de momentos
% Vector de términos conocidos
cein=[-AMMx-DY;
      -AMMy-DX;
      -AMMg-GM];
% Resolución
sol=ec\cein;
Ma=sol(3);
Ray=-sol(1);
Rax=sol(2);
%% ==============================================================
%% 7. REACCIONES EN APOYO B (equilibrio global)
%% ==============================================================
Rby=qa*L - (L^2*(qa - qb))/(2*L)-Ray;
Rbx=-Rax;
Mb=-Rax*h-Ma-(-((qa - qb) * L^3) / (6 * L) + (qa * L^2) / 2)+Ray*L;
%% Vector final de reacciones
R=[Rax;Ray;Ma;Rbx;Rby;Mb];
%% ==============================================================
%% 8. AJUSTE DE SIGNOS SI EL DOMINIO ESTÁ INVERTIDO
%% ==============================================================
if a<0
    Ma=-Ma;
    Rax=-Rax;
    Ray=-Ray;
    Mb=-(-Rax*h+Ma-(-((qa - qb) * L^3) / (6 * L) + (qa * L^2) / 2)-Ray*L);
    R=[-Rbx;-Rby;Mb;Rax;Ray;Ma];
end
R1=calculos_parabolico(y,dydxp,EI,Ac,M,V,a,b,L,h,qa,qb,DX,DY,GM)


K=k_rigidez(x1,x2,x3,y1,y2,y3,EI,Ac)
