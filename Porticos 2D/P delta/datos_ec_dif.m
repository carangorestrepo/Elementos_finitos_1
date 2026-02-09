%% PROPIEDADES DEL MATERIAL Y GEOMETRÍA
% Valores típicos para hormigón armado (ejemplo)
E = 24870062;      % Módulo de elasticidad [kN/m²]
G = 0.4*E;         % Módulo de cortante [kN/m²]
I = 0.4^4/12;      % Inercia de la sección [m?]
Ae = 0.4^2;        % Área efectiva [m²]
EI = E*I;          % Rigidez a flexión [kN·m²]
Ac = Ae*G*5/6;     % Rigidez a cortante reducida [kN]
g = 9.8066502;     % Aceleración de gravedad
rho = 2.4;         % Densidad del concreto [Mg/m^3]
L = 4;             % Longitud del elemento [m].
AE = Ae*E;         % Rigidez axial [kN]
k=0;
P=1000;
b1loc=20;
b2loc=20; 
q1loc=20;
q2loc=20;

puntos_graficas=500;
Ke = Ke_ec_dif(L,EI, AE,Ac,k,P,puntos_graficas);
[X1,Y1,M1,X2,Y2,M2] = fe_ec_dif(L, EI, Ac, AE, b1loc,b2loc, q1loc,q2loc,k,P)
