clc
clear
close all
% se definen algunas constantes que hacen el codigo mas legible

X = 1;
Y = 2;
G = 3;
XY = 12;
XYG=123;
E=24855578; %Modulo elaticidad viga
G=0.4*E;    %Modulo cortante

b1=0.3;     % ancho inicial viga
b2=0.5;     % ancho final viga
nb=1;       % exponente ancho final viga

h1=0.4;     % altura inicial viga
h2=0.8;     % altura final viga
nh=1;       % exponente altura viga

L=4;        % longitud viga
q1=25;      % carga vertical inicial viga
q2=60;      % carga vertical final viga
nq=1;       % exponente carga vertical final viga

b1a=25;     % carga axial inicial viga
b2a=45;     % carga axial final viga
naxi=2;     % exponenete carga axial final viga

sec=2;
%seccion C ==1
%Seccion rectangular  ==2
%Circular solida ==3 
%Seccion I ==4
%Seccion circular hueca ==5
%Tubular rectangular hueca ==6
%Perfil T==7

tf=0.09;
tw=0.06;
nudos_intemedios=100;

%% Se describen las propiedades de los materiales
%%% E G
Es=200000*1000;
v=0.3;
Gs=Es/(2*(1+v));
EG=[E,G];%[2,#materiales]];
puntos_graficas=100;
%[KE,Me]=KE_FE_seccion_variable(tf,tw,E,G,sec,b1,b2,nb,h1,h2,nh,L,q1,q2,nq,b1a,b2a,naxi);
[X1,Y1,M1,X2,Y2,M2] = fe_ec_dif(tf,tw,naxi,nq,nb,nh,L,E,G,b2,b1,h2,h1,sec,b2a,b1a, q2,q1,0,0,puntos_graficas);
FE=[X1,Y1,M1,X2,Y2,M2]';
KE = Ke_ec_dif(tf,tw,nb,nh,L,E,G,b2,b1,h2,h1,sec,0,0,puntos_graficas);
%% Matriz de rigidez y momentos de empotramiento Viga EE
[Ke,Me,Fe]=KE_FE_seccion_variable(tf,tw,E,G,b1,b2,nb,h1,h2,nh,L,q1,q2,nq,b1a,b2a,naxi,sec);

