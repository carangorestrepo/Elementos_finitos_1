% Ejemplo 11.23 Uribe Escamilla


%% Unidades en toneladas y metros

% se definen algunas constantes que hacen el codigo mas legible
NL1 = 1; NL2 = 2;
X = 1;   Y = 2;

nb = size(LaG,1);   % numero de barras (numero de filas de LaG)
%nn = size(gdl,1);   % numero de nodos  (numero de filas de gdl)

%% separo la memoria
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
%% ensamblo la matriz de rigidez global (K) y vector de fuerzas global (f)
for e = 1:nb  % para cada barra
   % saco los 6 gdls de la barra e
   
   idx{e} = [LaG(e,NL1)*3-2,LaG(e,NL1)*3-1,LaG(e,NL1)*3,LaG(e,NL2)*3-2,LaG(e,NL2)*3-1,LaG(e,NL2)*3];
   
   x1 = xnod(LaG(e,NL1), X);
   y1 = xnod(LaG(e,NL1), Y);
   x2 = xnod(LaG(e,NL2), X);
   y2 = xnod(LaG(e,NL2), Y);
   Le(e)  = hypot(x2-x1, y2-y1);
   c=(x2-x1)/Le(e);
   s=(y2-y1)/Le(e);
   cc(e)=c;
   ss(e)=s;
   ang(e) = atan2(y2-y1,x2-x1);
   % matriz de transformacion de coordenadas para la barra e
   %c = cosd(theta(e)); s = sind(theta(e));
   T{e} = [ c  s  0  0  0  0        
           -s  c  0  0  0  0        
            0  0  1  0  0  0
            0  0  0  c  s  0
            0  0  0 -s  c  0
            0  0  0  0  0  1 ];
   GJx = Jx(e)*G(e);    EIz = E(e)*Iz(e);    L=Le(e);  Acz = As2(e)*G(e);
   M=qa(e)*L^2/12;
   V=qa(e)*L/2;
   
   % matriz de rigidez local expresada en el sistema de coordenadas locales
   % para la barra e
             % T                                                  Mx              Vz
   Kloc =[[  GJx/L,                                                   0,                                    0, -GJx/L,                                                   0,                                    0]
          [      0,    (4*EIz*(Acz*L^2 + 3*EIz))/(L*(Acz*L^2 + 12*EIz)),      -(6*Acz*EIz)/(Acz*L^2 + 12*EIz),      0, -(2*EIz*(- Acz*L^2 + 6*EIz))/(L*(Acz*L^2 + 12*EIz)),       (6*Acz*EIz)/(Acz*L^2 + 12*EIz)];
          [      0,                     - (6*Acz*EIz)/(Acz*L^2 + 12*EIz),  (12*Acz*EIz)/(L*(Acz*L^2 + 12*EIz)),      0,                     -(6*Acz*EIz)/(Acz*L^2 + 12*EIz), -(12*Acz*EIz)/(L*(Acz*L^2 + 12*EIz))];
          [ -GJx/L,                                                   0,                                    0,  GJx/L,                                                   0,                                    0];
          [      0, -(2*EIz*(- Acz*L^2 + 6*EIz))/(L*(Acz*L^2 + 12*EIz)),      -(6*Acz*EIz)/(Acz*L^2 + 12*EIz),      0,    (4*EIz*(Acz*L^2 + 3*EIz))/(L*(Acz*L^2 + 12*EIz)),       (6*Acz*EIz)/(Acz*L^2 + 12*EIz)];
          [      0,                      (6*Acz*EIz)/(Acz*L^2 + 12*EIz), -(12*Acz*EIz)/(L*(Acz*L^2 + 12*EIz)),      0,                      (6*Acz*EIz)/(Acz*L^2 + 12*EIz),  (12*Acz*EIz)/(L*(Acz*L^2 + 12*EIz))]];
  if c==1
   Mq=[0,0, 0,0,0,0;
       0,0  M,0,0,0;
       0,0,-V,0,0,0;
       0,0, 0,0,0, 0;
       0,0, 0,0,0,-M;
       0,0, 0,0,0,-V];
  elseif s==1
    Mq=[0,0,0,0,0,0;
       0,0,-M,0,0,0;
       0,0,-V,0,0,0;
       0,0, 0,0,0, 0;
       0,0, 0,0,0, M;
       0,0, 0,0,0,-V];
      
  end
   vq=[0;0;1;0;0;1];
   
   fe{e}=T{e}*Mq*vq;
   Cey{e}=[0;qa(e);0;0;qa(e);0];
   % matriz de rigidez local en coordenadas globales
   Ke{e} = T{e}'*Kloc*T{e};            
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke{e}; % sumo Ke{e} a K global
   f(idx{e})        = f(idx{e})        + fe{e}; % sumo a f global
end


%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)

d = setdiff(1:ngdl, C); % d = [4 5 6 7 8 9];

%% extraigo las submatrices y especifico las cantidades conocidas
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

Kcc = K(C,C); Kcd = K(C,d); fd = f(C);
Kdc = K(d,C); Kdd = K(d,d); fc = f(d);

% desplazamientos para los gdls c = [1 2 3 10 11 12]
ac = zeros(length(C),1); % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc - Kdc*ac);    % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd; % calculo fuerzas de equilibrio desconocidas

% armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(ngdl,1);  q = zeros(ngdl,1);   % separo la memoria
a(C) = ac;   a(d) = ad; % desplazamientos 
q(C) = qd;  %q(d) = qc; % fuerzas nodales de equilibrio

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
    GJx = Jx(e)*G(e);
    EIz = E(e)*Iz(e);
    Acz = As2(e)*G(e);
    L=Le(e);
    
    [xa,V,M,w]=deformada(T{e}*a(idx{e}),qa(e),EIz,Acz,L,puntos_graficas);
    subplot(1,3,1)
    hold on
    if cc(e)==1
        y=ones(1,puntos_graficas)*y1;
        plot3(xa+x1,y,w,'b'), hold on,grid on
    elseif ss(e)==1       
        x=ones(1,puntos_graficas)*x1;
        plot3(x,xa+y1,w,'r'), hold on,grid on
    end
    view([-37 67])
    subplot(1,3,2)

    hold on
    if cc(e)==1
        y=ones(1,puntos_graficas)*y1;
        stem3(xa+x1,y,M,'fill','b','markersize',1), hold on,grid on
    elseif ss(e)==1       
        x=ones(1,puntos_graficas)*x1;
        stem3(x,xa+y1,M,'fill','r','markersize',1), hold on,grid on
    end
    view([-37 67])
    subplot(1,3,3)
    hold on
    if cc(e)==1
        y=ones(1,puntos_graficas)*y1;
        stem3(xa+x1,y,V,'fill','b','markersize',1), hold on,grid on
    elseif ss(e)==1       
        x=ones(1,puntos_graficas)*x1;
        stem3(x,xa+y1,V,'fill','r','markersize',1), hold on,grid on
    end
    view([-37 67])
        
        
end