clear
clc
format long g
syms C1 C2 C3 C4 C5 C6 alfa beta k Ac qy1 qy2 qx1 qx2 x L EI EA real
%% ========================================================================
% PROPIEDADES DEL MATERIAL Y DE LA SECCION
% ========================================================================
g = 9.8066502;           % [m/s^2] aceleracion de la gravedad
E = 24870062.3;          % [kPa] modulo de elasticidad
G = 0.4*E;               % [kPa] modulo de cortante
L  = 4;                  % [m] longitud
Ae = 0.4^2;              % [m^2] area
I  = 0.4^4/12;           % [m^4] inercia
Ac = Ae*G*5/6;           % rigidez efectiva a cortante kappa*G*A
k = 500;                 % coeficiente de balasto
EI = E*I;                % rigidez a flexion
EA = E*Ae;               % rigidez axial
%% ========================================================================
% PARAMETROS alfa Y beta
% Ecuacion diferencial homogenea:
% v'''' - k/Ac*v'' + k/EI*v = 0
% a = k/Ac
% b = k/EI
% alfa = 1/2*sqrt(a + 2*sqrt(b))
% beta = 1/2*sqrt(2*sqrt(b) - a)
% ========================================================================
a = k/Ac;
b = k/EI;
alfa = 0.5*sqrt(a + 2*sqrt(b));
beta = 0.5*sqrt(2*sqrt(b) - a);
%% ========================================================================
% CARGAS DISTRIBUIDAS
% ========================================================================
qx1 = 0;
qx2 = 0;
qy1 = 0;
qy2 = 0;
% Carga axial trapezoidal
qx = qx1 + (qx2-qx1)*x/L;
% Carga transversal trapezoidal
qy = qy1 + (qy2-qy1)*x/L;
%% ========================================================================
% SOLUCION TRANSVERSAL HOMOGENEA
% ========================================================================
f1 = cosh(alfa*x)*cos(beta*x);
f2 = cosh(alfa*x)*sin(beta*x);
f3 = sinh(alfa*x)*cos(beta*x);
f4 = sinh(alfa*x)*sin(beta*x);
vh = C1*f1 + ...
     C2*f2 + ...
     C3*f3 + ...
     C4*f4;
%% ========================================================================
% SOLUCION PARTICULAR TRANSVERSAL
% Para carga trapezoidal:
% vp = qy/k
% y como qy'' = 0:
% Mp = 0
% Vp = 0
% ========================================================================
vp = qy/k;
v = vh + vp;
%% ========================================================================
% MOMENTO FLECTOR
% M = EI*(v'' - k/Ac*v + qy/Ac)
% Para una carga trapezoidal la parte particular se cancela:
% M = EI*(vh'' - k/Ac*vh)
% ========================================================================
Mm = EI*(diff(vh,x,2) - k/Ac*vh);
%% ========================================================================
% CORTANTE
%
% V = M'
% ========================================================================
V = diff(Mm,x);
%% ========================================================================
% GIRO
% Ecuacion cinematica:
% v' = t - V/Ac
% luego:
% t = v' + V/Ac
% ========================================================================
t = diff(v,x) + V/Ac;
%% ========================================================================
% SOLUCION AXIAL
% Se adopta la convencion:
% N' = qx
% N = C5 + Np
% u' = N/EA
% ========================================================================
Np = qx1*x + (qx2-qx1)*x^2/(2*L);
N = C5 + Np;
%% Desplazamiento axial particular
up = qx1*x^2/(2*EA) + (qx2-qx1)*x^3/(6*EA*L);
%% Desplazamiento axial total
u = C6 + C5*x/EA + up;
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
% MATRIZ QUE RELACIONA CONSTANTES CON DESPLAZAMIENTOS
% Yd = [u;
%       v;
%       t]
% Solo se obtiene el coeficiente asociado a las constantes C1...C6.
% Las cargas particulares quedan fuera de esta matriz.
% ========================================================================
Bd = equationsToMatrix([u;
                        v;
                        t],Un);
Bd = simplify(Bd);
%% ========================================================================
% MATRIZ QUE RELACIONA CONSTANTES CON FUERZAS INTERNAS
% ========================================================================
Sf = equationsToMatrix([N;
                        V;
                        Mm],Un);
Sf = simplify(Sf);
%% ========================================================================
% EVALUACION EN LOS DOS EXTREMOS
% Vector de desplazamientos:
% d = [u1;
%      v1;
%      t1;
%      u2;
%      v2;
%      t2]
% ========================================================================
B0 = double(subs(Bd,x,0));
BL = double(subs(Bd,x,L));
B = [B0;
     BL];
%% ========================================================================
% FUERZAS NODALES
% Se adopta:
% f = [-N(0);
%      -V(0);
%      -M(0);
%       N(L);
%       V(L);
%       M(L)]
% Por ello se introducen los signos del extremo inicial.
% ========================================================================
S0 = double(subs(Sf,x,0));
SL = double(subs(Sf,x,L));
S = [-S0;
      SL];
%% ========================================================================
% MATRIZ DE RIGIDEZ
% d = B*C
% f = S*C
% C = inv(B)*d
% entonces:
% f = S*inv(B)*d
% K = S*inv(B)
% Numericamente es preferible NO usar inv(B).
% ========================================================================
K = S/B;
%% ========================================================================
% ELIMINACION DE PEQUENOS ERRORES NUMERICOS DE ASIMETRIA
% ========================================================================
K = 0.5*(K + K.');
%% Mostrar matriz
disp('Matriz de rigidez K =')
disp(K)