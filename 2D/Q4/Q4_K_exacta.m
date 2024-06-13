clear,clc        % borra la memoria y la pantalla
syms a b E nu t
%% Programa para calcular los modos de energia nula del solido
%a  = 1;       % [m]   ancho elemento
%b  = 1;       % [m]   altura elemento
%E  = 200;     % [GPa] modulo de elasticidad del elemento
%nu = 0.33;    %       coeficiente de Poisson
%t  = 0.10;    % [m]   espesor del elemento

%% Funciones de forma del rectangulo
syms r s
N1 = (1-r/a)*(1-s/b)/4;
N2 = (1+r/a)*(1-s/b)/4;
N3 = (1+r/a)*(1+s/b)/4;
N4 = (1-r/a)*(1+s/b)/4;

%% matriz constitutiva del elemento para TENSION PLANA
D = [ E/(1-nu^2)     E*nu/(1-nu^2)  0
      E*nu/(1-nu^2)  E/(1-nu^2)     0
      0              0              E/(2*(1+nu)) ];

  %% matriz de deformaciones
B1 = [diff(N1,r)   0         
      0            diff(N1,s)
      diff(N1,s)   diff(N1,r)];
      
B2 = [diff(N2,r)   0
      0            diff(N2,s)
      diff(N2,s)   diff(N2,r)];
      
B3 = [diff(N3,r)   0
      0            diff(N3,s)
      diff(N3,s)   diff(N3,r)];
      
B4 = [diff(N4,r)   0
      0            diff(N4,s)
      diff(N4,s)   diff(N4,r)];

% se ensambla la matriz de deformaciones
B = [B1 B2 B3 B4];

%% realizo la integracion para calcular la matriz K exactamente
K = simplify(int(int(B.'*D*B*t, r, -a,a), s,-b,b));

%% Se definen las funciones de forma unidimensionales y sus derivadas
N1 = (1-r/a);
N2 = (1+r/a);
idx = [ 1 2 3 ];
Nijk = [ N(1)    0    N(2)     0    N(3)       0    N(4)      0   
           0  N(1)       0  N(2)       0    N(3)       0    N(4)];



