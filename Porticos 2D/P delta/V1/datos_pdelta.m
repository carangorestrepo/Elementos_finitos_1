%% Seleccion de desplazamientos
%% Se describen las propiedades de los materiales
clc
clear
syms C1 C2 C3 C4 C5 C6 x  alfa beta L qx1 qx2 qy1 qy2 EI Ac k EA;
%%{
g=9.8066502;      % aceleracion de la gravedad
E = 24870062.3;       % [kPa] modulo de elasticidad del concreto (NSR 10  - art C.8.5)
G = 0.4*E;                    % [kPa] modulo de cortante (NSR 10  - art CR8.8.2)
L = 4;                        % [m]   longitud
Ae = 0.4^2;                   % [m^2] area  
I = 0.4^4/12;                 % [m^4] inercia
Ac=Ae*G*5/6;                  % coeficiente de correccion del cortante para seccion rectangular
k=500;
qy1=25;
qy2=45;
qx1=20;
qx2=25;
EI=E*I;
AE=Ae*E;
P=1000;
la=(P/EI)^(1/2);
EA=E*Ae;
%}
%% ========================================================================
% CARGAS DISTRIBUIDAS
% ========================================================================
% Carga axial trapezoidal
qx = qx1 + (qx2-qx1)*x/L;
% Carga transversal trapezoidal
qy = qy1 + (qy2-qy1)*x/L;
%% ========================================================================
% SOLUCION HOMOGENEA TRANSVERSAL
% ========================================================================
vh = C1 + C2*x ...
  + C3*cos(la*x) ...
  + C4*sin(la*x); 
%% ========================================================================
% MOMENTO HOMOGENEO
% Para la solucion homogenea:
% ========================================================================
Mh = simplify(EI*(diff(vh,x,2)));
%% ========================================================================
% CORTANTE HOMOGENEO
% ========================================================================
Vh = simplify(diff(Mh,x) + P*diff(vh,x));
%% ========================================================================
% GIRO HOMOGENEO

% ========================================================================
th = simplify(diff(vh,x) + Vh/Ac);
%% ========================================================================
% SOLUCION PARTICULAR TRANSVERSAL
% Para carga trapezoidal:
% qy = qy1 + (qy2-qy1)*x/L

% ========================================================================
r = (qy2-qy1)/L;
q = qy1 + r*x;
vp = - qy1*x^2/(2*P) ...
  - r*x^3/(6*P);
%
% ========================================================================
Mp = simplify(EI*(diff(vp,x,2) - q/Ac));

%% ========================================================================
% CORTANTE PARTICULAR
% ========================================================================
Vp = simplify(diff(Mp,x) + P*diff(vp,x));

%% ========================================================================
% GIRO PARTICULAR
% ========================================================================
tp = simplify(diff(vp,x) + Vp/Ac);
%% ========================================================================
% SOLUCION AXIAL HOMOGENEA
% Para ausencia de carga axial:
% N' = 0
% entonces:
% Nh = C5
% u' = N/EA
% uh = C5*x/EA + C6
% ========================================================================
Nh = C5;
uh = C5*x/EA + C6;
%% ========================================================================
% SOLUCION AXIAL PARTICULAR
% N' = -qx
% entonces:
% Np = -integral(qx dx)
% ========================================================================
Np = -int(qx,x);
Np = simplify(Np);
%% ========================================================================
% DESPLAZAMIENTO AXIAL PARTICULAR
%
% u' = N/EA
% ========================================================================
up = int(Np,x)/EA;
up = simplify(up);
%% ========================================================================
% VECTOR DE CONSTANTES DE INTEGRACION
% ========================================================================
Un = [C1;
      C2;
      C3;
      C4;
      C5;
      C6];
%% ========================================================================
% MATRIZ FUNDAMENTAL Afun
% Solucion homogenea:
% Yh = Afun*Un
% donde:
% Yh = [uh;
%       vh;
%       th;
%       Nh;
%       Vh;
%       Mh]
% ========================================================================
Afun = equationsToMatrix( ...
    [uh;
     vh;
     th;
     Nh;
     Vh;
     Mh], ...
     Un);
Afun = simplify(Afun);
%% ========================================================================
% VECTOR PARTICULAR Yp
% ========================================================================
Yp = [up;
      vp;
      tp;
      Np;
      Vp;
      Mp];
Yp = simplify(Yp);
%% ========================================================================
% SOLUCION TOTAL
% Y = Afun*Un + Yp
% ========================================================================
Y = Afun*Un + Yp;
Y = simplify(Y);
%% Variables totales
u = Y(1);
v = Y(2);
t = Y(3);
N = Y(4);
V = Y(5);
M = Y(6);
D = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];
%% Seleccion de fuerzas internas
R = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];
%% Signos en extremo inicial
Si = [-1  0  0;
       0  1  0;
       0  0 -1];
%% Signos en extremo final
Sj = [ 1  0  0;
       0 -1  0;
       0  0  1]; 
[Ke,feq,B,S,dp,fp] = matriz_elemento_EDO(Afun,Yp,x,L,D,R,Si,Sj);
Ke=double(Ke);
feq=double(feq);
%% ========================================================================
% FUNCIONES DE FORMA EXACTAS
%
% La solucion general es:
%
% Y(x) = Afun(x)*Un + Yp(x)
%
% y los desplazamientos nodales cumplen:
%
% d = B*Un + dp
%
% donde:
%
% d = [u1;
%      v1;
%      t1;
%      u2;
%      v2;
%      t2]
%
% Despejando las constantes:
%
% Un = inv(B)*(d-dp)
%
% Sustituyendo:
%
% Y(x) = Afun*inv(B)*d ...
%      + Yp ...
%      - Afun*inv(B)*dp
%
% Por tanto se define:
%
% H(x) = Afun*inv(B)
%
% Esta matriz interpola TODO el vector de estado.
%
% Numericamente/simbolicamente es preferible utilizar:
%
% H = Afun/B
%
% en lugar de inv(B).
% ========================================================================
H = simplify(Afun/B);
%% ========================================================================
% FUNCIONES DE FORMA DE LOS DESPLAZAMIENTOS
%
% El operador D selecciona:
%
% [u;
%  v;
%  t]
%
% del vector completo:
%
% Y = [u;
%      v;
%      t;
%      N;
%      V;
%      M]
%
% Entonces:
%
% [u(x);
%  v(x);
%  t(x)] =
%
% Nfun(x)*d + dq(x)
%
% donde:
%
% Nfun(x) = D*H(x)
% ========================================================================

Nfun = simplify(D*H);


%% ========================================================================
% SOLUCION PARTICULAR CORREGIDA
%
% La solucion particular Yp por si sola generalmente NO satisface
% desplazamientos nodales nulos.
%
% Por eso debe corregirse mediante:
%
% Yq(x) = Yp(x) - H(x)*dp
%
% Esta solucion cumple:
%
% D*Yq(0) = 0
% D*Yq(L) = 0
%
% y representa la deformada producida exclusivamente por las cargas
% cuando todos los grados de libertad nodales son cero.
% ========================================================================

Yq = simplify(Yp - H*dp);


%% Particular de desplazamientos

dq = simplify(D*Yq);


%% ========================================================================
% FUNCIONES DE FORMA INDIVIDUALES
%
% Nfun es una matriz 3x6:
%
%           u1   v1   t1   u2   v2   t2
%
% u(x) = [ Nu1  Nv1  Nt1  Nu2  Nv2  Nt2 ]*d
%
% v(x) = [ ...                         ]*d
%
% t(x) = [ ...                         ]*d
%
% Para mayor claridad se separan las tres filas.
% ========================================================================
Nu = simplify(Nfun(1,:));
Nv = simplify(Nfun(2,:));
Nt = simplify(Nfun(3,:));
%% ========================================================================
% INTERPOLACION DE LOS DESPLAZAMIENTOS
% ========================================================================
syms u1 v1 t1 u2 v2 t2 real
de = [u1;
      v1;
      t1;
      u2;
      v2;
      t2];
%% Desplazamiento axial
ux = simplify(Nu*de + dq(1));

%% Desplazamiento transversal
vx = simplify(Nv*de + dq(2));
%% Giro
tx = simplify(Nt*de + dq(3));
%% ========================================================================
% FUNCIONES DE FORMA DEL VECTOR DE ESTADO COMPLETO
% Como:
% Y(x) = H(x)*d + Yq(x)
% H contiene directamente las funciones exactas para:
% u
% v
% t
% N
% V
% M
% ========================================================================
HY = H;
%% Funciones para fuerza axial
HN = simplify(HY(4,:));

%% Funciones para cortante
HV = simplify(HY(5,:));
%% Funciones para momento
HM = simplify(HY(6,:));
%% ========================================================================
% RECUPERACION DEL ESTADO COMPLETO
% ========================================================================
Ytotal = simplify(HY*de + Yq);
u_forma = Ytotal(1);
v_forma = Ytotal(2);
t_forma = Ytotal(3);
N_forma = Ytotal(4);
V_forma = Ytotal(5);
M_forma = Ytotal(6);
%% ========================================================================
% COMPROBACION DE LAS FUNCIONES DE FORMA EN LOS NODOS
%
% Las funciones de desplazamiento deben satisfacer:
%
% Nfun(0) =
% [1 0 0 0 0 0
%  0 1 0 0 0 0
%  0 0 1 0 0 0]
% Nfun(L) =
% [0 0 0 1 0 0
%  0 0 0 0 1 0
%  0 0 0 0 0 1]
% ========================================================================
N0 = simplify(subs(Nfun,x,0));
NL = simplify(subs(Nfun,x,L));
disp('Funciones de forma en x=0')
disp(N0)
disp('Funciones de forma en x=L')
disp(NL)
%% ========================================================================
% COMPROBACION DE LA SOLUCION PARTICULAR CORREGIDA
% Debe cumplirse:
% dq(0) = 0
% dq(L) = 0
% ========================================================================
dq0 = simplify(subs(dq,x,0));
dqL = simplify(subs(dq,x,L));
disp('Solucion particular corregida en x=0')
disp(dq0)
disp('Solucion particular corregida en x=L')
disp(dqL)