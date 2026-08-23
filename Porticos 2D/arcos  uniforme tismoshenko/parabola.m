clc
clear
syms x1 x2 x3 y1 y2 y3 x C1 C2 xa
%calculo de ecuacion parabola

%% 2. COORDENADAS DE LA PARÁBOLA (tres puntos)
% Puntos de paso del arco (pueden modificarse)
x1 = 0;   y1 = 0;      % punto inicial (apoyo A)
x2 = 2.5; y2 = 1.5;    % punto intermedio
x3 = 4;   y3 = 0.6;     % punto final (apoyo B)
%%{

X = [x1^2,x1,1;
     x2^2,x2,1; 
     x3^2,x3,1]; 
Y = [y1;y2;y3];
coe=X^(-1)*Y;
%}
y=coe(1)*x^2+coe(2)*x+coe(3);%funcion de arco 

%% 4. FUNCIONES GEOMÉTRICAS DEL ARCO (derivada, longitud diferencial)
dydx   = diff(y, x, 1);      % pendiente dy/dx
ds     = sqrt(1 + dydx^2);   % diferencial de longitud de arco

tanphi = dydx;                % tan(?)
% Componentes unitarias: Vu = cos(?), Pu = sen(?)
% (el signo sigue la convención de que el arco crece en x)

Vu = 1 / sqrt(1 + tanphi^2);     % = cos(?)
Pu = tanphi / sqrt(1 + tanphi^2); % = sin(?)

%% propieades de seccion tranversal seccion viga
bf = 25/1000; % ancho de seccion viga (m)
d = 50/1000;  % altura de seccinn viga (m)
E=210000000;  % modulo de elaticidad elemento (kPa)
Ae=bf*d;      % Area elemento (m2)
I=bf*d^3/12;  % inercia elemento (m4)
EI=E*I;
EA=E*Ae;
G=0.4*E;
alfaGA = (5/6)*Ae*G; % rigidez a cortante (factor de forma 5/6 para sección rectangular)
qa = 20;          % valor inicial de la carga (en x=0)
qb = 25;         % valor final (en x=L)
a  = 0;          % límite inferior de integración (x inicial)
b  = x3;         % límite superior (x final)
L  = x3;         % longitud proyectada en x
h  = y3;         % altura del arco en el apoyo B (para brazo de palanca)

% Diagramas de momento flector (M) y cortante (V) de la carga trapezoidal
% sobre una viga recta simplemente apoyada (NO sobre el arco). Estos
% diagramas se proyectan luego sobre las direcciones del arco.
M = -((qa - qb) * x^3) / (6 * L) + (qa * x^2) / 2;
V = qa*x - (x^2 * (qa - qb)) / (2*L);

Vz = Vu*V;% componente  cortante 
Pz = Pu*V;% componente  axial

%% 6. DESPLAZAMIENTOS CONOCIDOS EN EL APOYO A (condiciones de contorno)
DX = 0;   % desplazamiento horizontal conocido (cero si es empotramiento)
DY = 0;   % desplazamiento vertical conocido
GM = 0;   % giro conocido

%% 1. DEFORMACIÓN VERTICAL (Qy = 1)
AMxx=area_cuadraturas(a,b,matlabFunction((x*x/EI+Vu/alfaGA)*ds,"Vars", x)); % Deformación vertical por carga unitaria vetical
AMxy=area_cuadraturas(a,b,matlabFunction((x*y/EI+Pu/alfaGA)*ds,"Vars", x)); % Deformación horizontal por carga unitaria vetical
AMxg=area_cuadraturas(a,b,matlabFunction((x/EI)*ds,"Vars", x));             % giro por carga unitaria vetical

%% 2. DEFORMACIÓN HORIZONTAL (Qx = 1)
%%NOta los diagramas de cortante y carga axial se invierten 
AMyx=area_cuadraturas(a,b,matlabFunction((y*x/EI+Pu/alfaGA)*ds,"Vars", x));   % Deformación vertical por carga unitaria horizonal
AMyy = area_cuadraturas(a,b,matlabFunction((y*y/EI+Vu/alfaGA)*ds,"Vars", x)); % Deformación horizontal por carga unitaria horizonal
AMyg = area_cuadraturas(a,b,matlabFunction((y/EI)*ds,"Vars", x));             % giro por carga por carga unitaria horizonal
 
%% 3. DEFORMACIÓN GIRO (QG = 1)
AMx = area_cuadraturas(a,b,matlabFunction((x/EI)*ds,"Vars", x));             % Deformación vertical por carga unitaria momento
AMy = area_cuadraturas(a,b,matlabFunction((y/EI)*ds,"Vars", x));             % Deformación horizontal por carga unitaria momento
AMg = area_cuadraturas(a,b,matlabFunction((1/EI)*ds,"Vars", x));             % giro por carga  por carga unitaria momento

%% 4. DEFORMACIÓN CARGA EXTERNA Y
AMMx = area_cuadraturas(a,b,matlabFunction((M*x/EI+ Vz/alfaGA)*ds,"Vars", x)); % Deformación vertical por carga unitaria vetical
AMMy = area_cuadraturas(a,b,matlabFunction((M*y/EI+ Pz/alfaGA)*ds,"Vars", x)); % Deformación horizontal por carga unitaria vetical
AMMg = area_cuadraturas(a,b,matlabFunction((M/EI)*ds,"Vars", x));              % giro por carga unitaria vetical

%% ==============================================================
%% 6. SISTEMA DE ECUACIONES (compatibilidad)
%% ==============================================================
% Matriz de coeficientes
%% 8. SISTEMA DE ECUACIONES DE COMPATIBILIDAD (método de las fuerzas)
% Matriz de flexibilidad [F] (3x3)

ec=[AMxx ,AMxy ,AMxg;  % Suma de fuerzas en Y
    AMyx ,AMyy ,AMyg;  % Suma de fuerzas en X
    AMx  ,AMy  ,AMg];  % Suma de momentos
% Vector de términos conocidos
cein=[-AMMx-DY;
      -AMMy-DX;
      -AMMg-GM];
% Condiciones de compatibilidad: [F]{X} + {D_ext} = {desplazamientos conocidos}
% Los desplazamientos conocidos en A son [DY; DX; GM] (observar signos según
% la dirección positiva de las redundantes). En este código, las redundantes
% son: Ray (vertical), Rax (horizontal), Ma (momento). La ecuación resulta:
%   [F] * [Ray; Rax; Ma] + D_ext = [DY; DX; GM]
% Por tanto:
%   [F] * [Ray; Rax; Ma] = [DY; DX; GM] - D_ext
% En nuestro caso DY=DX=GM=0, así que:
%   [F] * [Ray; Rax; Ma] = - D_ext
% Resolviendo:

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




