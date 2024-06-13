
function [KE,Me,FE]=KE_FE_seccion_variable(tf,tw,E,G,b1,b2,nb,h1,h2,nh,L,q1loc,q2loc,nq,b1loc,b2loc,naxi,sec)

%clc
%clear
%syms q1loc q2loc nw L nq x E G As2
%syms b1 b2 nb L 
%syms h1 h2 nh L 
%{
E=24855578; %Modulo elaticidad viga
G=0.4*E;    %Modulo cortante
As2=5/6;    % area cortante 

b1=0.4;     % ancho inicial viga
b2=0.5;     % ancho final viga
nb=2;       % exponente ancho final viga

h1=0.4;     % altura inicial viga
h2=0.8;     % altura final viga
nh=2;       % exponente altura viga

L=4;        % longitud viga
q1loc=25;   % carga vertical inicial viga
q2loc=45;   % carga vertical final viga
nq=3;       % exponente carga vertical final viga

b1loc=25;   % carga axial inicial viga
b2loc=45;   % carga axial final viga
naxi=3;     % exponenete carga axial final viga
%}

g=0;        %deformacion vertical o axial
v=0;        %angulo de giro

[k2] =funcion(tf,tw,b1,b2,nb,h1,h2,nh,0 ,0 ,nq,L,E ,G,sec, 0,1 ,0,L);%% giroi==1 apoyada empotrada 
[k3] =funcion(tf,tw,b1,b2,nb,h1,h2,nh,0 ,0 ,nq,L,E ,G,sec,-1,0 ,0,L);%%deformadai==1   empotrada empotrada
[k5] =funcion(tf,tw,b2,b1,nb,h2,h1,nh,0 ,0 ,nq,L,E ,G,sec, 0,-1,L,0);%% girof==1 apoyada empotrada 
[k6] =funcion(tf,tw,b2,b1,nb,h2,h1,nh,0 ,0 ,nq,L,E ,G,sec,-1,0 ,L,0);%%deformadaf==1   empotrada empotrada
 
[k1]= funcion1(tf,tw,b1,b2,nb,h1,h2,nh,0,0,naxi,L,E,G,1,0,L,sec);    % axial
 
%fuezas nodales y equivalentes  flexion y cortante
[MEE]=funcion(tf,tw,b1,b2,nb,h1,h2,nh,q1loc,q2loc,nq,L,E,G,sec,g,v,0,L);% Flexion y coetante

%%fuerzas nodales equivalentes axiales
[R]=  funcion1(tf,tw,b1,b2,nb,h1,h2,nh,b1loc,b2loc,naxi,L,E,G,0,0,L,sec);     % axial
%% Mareiz de rigidez
KE=[[k1(1);0;0;k1(2);0;0],...
   [0;k2(1:2);0;k2(3:4)],...
   [0;k3(1:2);0;k3(3:4)],...
   [k1(2);0;0;k1(1);0;0],...
   [0;k5(1);k5(4);0;k5(3);k5(2)],...
   [0;k6(1);k6(4);0;k6(3);k6(2)]];
%% Fuerzas nodales equivalentes
FE=[R(1);MEE(1:2,1);R(2);MEE(3:4,1)];

%comprobacion por matlab
%puntos_graficas=101;
%Ke = Ke_ec_dif(nb,nh,L,E,G,b2,b1,h2,h1,As2,0,0,puntos_graficas);
%[X1,Y1,M1,X2,Y2,M2] = fe_ec_dif(naxi,nq,nb,nh,L,E,G,b2,b1,h2,h1, As2, b2loc,b1loc, q2loc,q1loc,0,0,puntos_graficas);
%Me=[X1,Y1,M1,X2,Y2,M2]';
X1=FE(1);
Y1=FE(2);
M1=FE(3);
X2=FE(4);
Y2=FE(5);
M2=FE(6);
    % matriz de mometos de empotramiento
Me=[-X1 ,0  ,0,0  ,0  ,0;
     0  ,Y1,0,0   ,0  ,0;
     0  ,-M1,0,0  ,0  ,0;
     0  ,0  ,0,X2 ,0  ,0;
     0  ,0  ,0,0  ,-Y2,0;
     0  ,0  ,0,0  ,M2 ,0];
 
 FE=[X1,Y1,M1,X2,Y2,M2]';
