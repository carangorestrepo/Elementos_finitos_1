clc
clear
close
% Ejemplo 11.23 Uribe Escamilla

%% Unidades en toneladas y metros
% se definen algunas constantes que hacen el codigo mas legible
% se definen algunas constantes que hacen el codigo mas legible
G = 3;
YG=12;
YGG=123;
NL1 = 1; NL2 = 2;
X = 1;   Y = 2;
bf=0.3;
d=0.4;
Ae=bf*d;
E=4700*sqrt(28)*1000;
Ge=0.4*E;
I=bf*d^3/12;
J=d/2*(bf/2)^3*(16/3-3.36*(bf/2)/(d/2)*(1-(bf/2)^4/(12*(d/2)^4)));
Ac=Ae*Ge*5/6;
EI=E*I;
GJ=Ge*J;
puntos_graficas=50;
xnod=[0,0;1.2,0;2.4,1;0,3];
LaG=[1,2;2,3;1,4];

nb = size(LaG,1);   % numero de barras (numero de filas de LaG)
%nn = size(gdl,1);   % numero de nodos  (numero de filas de gdl)
ngdl=size(xnod,1)*3;
%% separo la memoria
qe=[30,30;
    30,30;
    20,20];
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[YGG,3;
            YGG,4];         
%% escalas de Dibujo la estructura y su deformada
esc_def    = 500;          % escalamiento de la deformada
esc_faxial = 0.001;      % escalamiento del diagrama de axiales
esc_V      = 0.007;       % escalamiento del diagrama de cortantes
esc_M      = 0.01;      % escalamiento del diagrama de momentos

%% numero nudos
nudos = size(xnod,1);
%% a
napoyo=size(tipo_apoyo,1);
fe = cell(nb,1);      %% fuerzas nodales equivalentes para las diferentes barras
K   = zeros(ngdl);    % matriz de rigidez global
Le   = zeros(nb,1);   % Longitud elementos
cc   = zeros(nb,1);   % Longitud elementos
ss   = zeros(nb,1);   % Longitud elementos
f   = zeros(ngdl,1);  % vector de fuerzas nodales equivalentes global
Ke  = cell(nb,1);     % matriz de rigidez local en coordenadas globales
T   = cell(nb,1);     % matriz de transformacion de coordenadas
idx = cell(nb,1);     % almacena los 6 gdls de las barras
qe_glob = cell(nb,1);
qe_loc  = cell(nb,1);
Cey = cell(nb,1); 
ang = zeros(nb,1); 
ca = zeros(nudos*3,1);
%% ensamblo la matriz de rigidez global (K) y vector de fuerzas global (f)
for e = 1:nb  % para cada barra
   % saco los 6 gdls de la barra e
   idx{e} = [LaG(e,NL1)*3-2,LaG(e,NL1)*3-1,LaG(e,NL1)*3,LaG(e,NL2)*3-2,LaG(e,NL2)*3-1,LaG(e,NL2)*3];
   x1 = xnod(LaG(e,NL1), X);
   y1 = xnod(LaG(e,NL1), Y);
   x2 = xnod(LaG(e,NL2), X);
   y2 = xnod(LaG(e,NL2), Y);
   Le(e)  = hypot(x2-x1, y2-y1);
   L=Le(e);
   c=(x2-x1)/Le(e);
   s=(y2-y1)/Le(e);
   cc(e)=c;
   ss(e)=s;
   ang(e) = atan2(y2-y1,x2-x1);
   % matriz de transformacion de coordenadas para la barra e
   %c = cosd(theta(e)); s = sind(theta(e));
   T{e} = [c,-s,0,0, 0,0;
           s, c,0,0, 0,0;
           0, 0,1,0, 0,0;
           0, 0,0,c,-s,0;
           0, 0,0,s, c,0;
           0, 0,0,0, 0,1];
   %GJx = Jx(e)*G(e);    EIz = E(e)*Iz(e);    L=Le(e);  Acz = As2(e)*G(e);
   %M=qa(e)*L^2/12;
   %V=qa(e)*L/2;
   
   % matriz de rigidez local expresada en el sistema de coordenadas locales
   % para la barra e
   % T                                                  Mx              Vz
   Kloc =[[                                (4*EI*(Ac*L^2 + 3*EI))/(L*(Ac*L^2 + 12*EI)),     0,       (6*Ac*EI)/(Ac*L^2 + 12*EI),                               -(2*EI*(- Ac*L^2 + 6*EI))/(L*(Ac*L^2 + 12*EI)),     0,      -(6*Ac*EI)/(Ac*L^2 + 12*EI)];
          [                                                                          0,  GJ/L,                                0,                                                                            0, -GJ/L,                                0];
          [                                                 (6*Ac*EI)/(Ac*L^2 + 12*EI),     0,  (12*Ac*EI)/(L*(Ac*L^2 + 12*EI)),                                                   (6*Ac*EI)/(Ac*L^2 + 12*EI),     0, -(12*Ac*EI)/(L*(Ac*L^2 + 12*EI))];
          [ (6*Ac*EI*L)/(Ac*L^2 + 12*EI) - (4*EI*(Ac*L^2 + 3*EI))/(L*(Ac*L^2 + 12*EI)),     0,       (6*Ac*EI)/(Ac*L^2 + 12*EI), (2*EI*(- Ac*L^2 + 6*EI))/(L*(Ac*L^2 + 12*EI)) + (6*Ac*EI*L)/(Ac*L^2 + 12*EI),     0,      -(6*Ac*EI)/(Ac*L^2 + 12*EI)];
          [                                                                          0, -GJ/L,                                0,                                                                            0,  GJ/L,                                0];
          [                                                -(6*Ac*EI)/(Ac*L^2 + 12*EI),     0, -(12*Ac*EI)/(L*(Ac*L^2 + 12*EI)),                                                  -(6*Ac*EI)/(Ac*L^2 + 12*EI),     0,  (12*Ac*EI)/(L*(Ac*L^2 + 12*EI))]];
    q1=qe(e,1);
    q2=qe(e,2);
     M1=((L^2*q1)/12 - (L^4*q1)/30 + (L^4*q2)/30);
     M2=((L^4*q2)/30 - (L^4*q1)/30 - (5*L^2*q1)/12 + L^3*((L*q1)/3 - (L*q2)/3) + L*((L*q1)/2 + (3*L^3*q1)/20 - (3*L^3*q2)/20 - L^2*((L*q1)/2 - (L*q2)/2)));
     V1=-((3*L^3*q1)/20 - (L*q1)/2 - (3*L^3*q2)/20);
     V2=((L*q1)/2 + (3*L^3*q1)/20 - (3*L^3*q2)/20 - L^2*((L*q1)/2 - (L*q2)/2)); 

   Mq=[0,0,-M1,0,0,0;
       0,0  0,0,0,0;
       0,0,-V1,0,0,0;
       0,0, 0,0,0,M2;
       0,0, 0,0,0,0;
       0,0, 0,0,0,-V2];
   vq=[0;0;1;0;0;1];
   fe{e}=T{e}'*Mq*vq;
   Cey{e}=[0;0;qe(e);0;0;qe(e)];
   % matriz de rigidez local en coordenadas globales
   Ke{e} = T{e}'*Kloc*T{e};            
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke{e}; % sumo Ke{e} a K global
   f(idx{e})        = f(idx{e})        + fe{e}; % sumo a f global
end
%% grados de libertad que no tienen ceros
smy=sum(abs(K),1);
[~,fy]=find(smy~=0);
% extraigo las submatrices y especifico las cantidades conocidas
for ap=1:napoyo
    if tipo_apoyo(ap,1)==YGG
        ca(tipo_apoyo(ap,2)*3-2:tipo_apoyo(ap,2)*3,1)=tipo_apoyo(ap,2)*3-2:tipo_apoyo(ap,2)*3;
    elseif tipo_apoyo(ap,1)==YGX
        ca(tipo_apoyo(ap,2)*3-2:tipo_apoyo(ap,2)*3-1,1)=tipo_apoyo(ap,2)*3-2:tipo_apoyo(ap,2)*3-1;
    elseif tipo_apoyo(ap,1)==YGY
        ca(tipo_apoyo(ap,2)*3-1,1)=tipo_apoyo(ap,2)*3-1;
    end
end
Kel=K(fy,fy);
fel=f(fy,1);


%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
    ci = setdiff(ca,0); % desplazamientos conocidos
d = setdiff(1:ngdl, ca); % d = [4 5 6 7 8 9];

%% extraigo las submatrices y especifico las cantidades conocidas
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

Kcc = K(ci,ci); Kcd = K(ci,d); fd = f(ci);
Kdc = K(d,ci); Kdd = K(d,d); fc = f(d);

% desplazamientos para los gdls c = [1 2 3 10 11 12]
ac = zeros(length(ci),1); % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc - Kdc*ac);    % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd; % calculo fuerzas de equilibrio desconocidas

% armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(ngdl,1);  q = zeros(ngdl,1);   % separo la memoria
a(ci) = ac;   a(d) = ad; % desplazamientos 
q(ci) = qd;  %q(d) = qc; % fuerzas nodales de equilibrio

%fprintf('Desplazamientos de los nodos en coord. globales = \n'); a

%% imprimo las fuerzas internas en cada barra referidas a las coordenadas
%% globales
for e = 1:nb % para cada barra
    qe_glob{e} = Ke{e}*a(idx{e}) - fe{e};
    qe_loc{e}=T{e}*qe_glob{e};
    x1 = xnod(LaG(e,NL1), X);
    y1 = xnod(LaG(e,NL1), Y);
    x2 = xnod(LaG(e,NL2), X);
    y2 = xnod(LaG(e,NL2), Y);
    %GJx = Jx(e)*G(e);
    %EI = E(e)*Iz(e);
    %Ac = As2(e)*G(e);
    %L=Le(e); 
    deformada(GJ,EI,Ac,Cey{e},T{e}*a(idx{e}),Le(e),ang(e),esc_def,esc_faxial,esc_V,esc_M,x1,y1,x2,y2,puntos_graficas)     
end
figure(2); hold on; title('Deformada exagerada');    xlabel('x, m'); ylabel('y, m'); axis equal
figure(3); hold on; title('Fuerza cortante [kN]');   xlabel('x, m'); ylabel('y, m'); axis equal
figure(4); hold on; title('Momento flector [kN-m]'); xlabel('x, m'); ylabel('y, m'); axis equal