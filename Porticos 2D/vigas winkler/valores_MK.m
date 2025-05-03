
E = 24870062.3;                 % [kPa] modulo de elasticidad del concreto (NSR 10  - art C.8.5)
G = 0.4*E;                    % [kPa] modulo de cortante (NSR 10  - art CR8.8.2)
L = 4;               % [m]   longitud
Ae = 0.4^2;                   % [m^2] area  
I = 0.4^4/12;                 % [m^4] inercia
Ac=Ae*G*5/6;                  % coeficiente de correccion del cortante para seccion rectangular
k=0;%k=500;
EI=E*I;
puntos_graficas=2000;
wa = 20;              % [kN/m] carga uniforme distribuida vertical inicial
wb = 20;              % [kN/m] carga uniforme distribuida vertical final
P=1000;
AE=Ae*E;
[X1,Y1,M1,X2,Y2,M2] = fe_ec_dif(L, EI, Ac, AE, wa,wb, wa,wb,k,P,puntos_graficas);
Ke = Ke_ec_dif(L,EI, AE,Ac,k,P,puntos_graficas);

Mv=[X1,Y1,M1,X2,Y2,M2]'