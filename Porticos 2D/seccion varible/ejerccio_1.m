clc
clear
close all
% se definen algunas constantes que hacen el codigo mas legible
X = 1;
Y = 2;
G = 3;
XY = 12;
XYG=123;

%% Se describen las propiedades de los materiales
%%% E G
Es=2*10^7;
v=0.3;
Gs=Es/(2*(1+v));
EG=[Es,Gs;%[2,#materiales]
    Es,Gs];
%% tipo de material por elemento 
mat_elem = ones(1,5); %[1,#elemetos]
%% ingreso coordenadas de nudos de estructura % [#coordenadas,2]

xyc= [0,0;
      1,1.333;
      3,4;
      6,4;
      7,2.667;
      9,0];

%% Nudos que conectan elementos
% x y %% [#nuros,2]

ninf=[1,2;
      2,3;
      3,4;
      4,5;
      5,6];

%%  Propiedades de la geometria secciones 
%b1=0.3;     % ancho inicial viga
%b2=0.2;     % ancho final viga
nb=1;       % exponente ancho final viga

%h1=0.5;     % altura inicial viga
%h2=0.8;     % altura final viga
nh=1;       % exponente altura viga

%q1=25;      % carga vertical inicial viga
%q2=30;      % carga vertical final viga
nq=1;       % exponente carga vertical final viga

%b1a=25;     % carga axial inicial viga
%b2a=30;     % carga axial final viga
naxi=1;     % exponenete carga axial final viga

sec=2;
%seccion C ==1
%Seccion rectangular  ==2
%Circular solida ==3 
%Seccion I ==4
%Seccion circular hueca ==5
%Tubular rectangular hueca ==6
%Perfil T==7
%Perfil L==8
tf=0.09;
tw=0.06;

  %% b1,b2,h1 h2,nb,bh
bh=[1,1,1,0.9667,nb,nh;
    1,1,0.9667,0.9,nb,nh;
    1,1,0.9,0.8,nb,nh;
    1,1,0.8,0.7667,nb,nh;
    1,1,0.7667,0.70,nb,nh;];
%% cargas aplicadas distribuida (gdl carga)
% wix   , wfx , direccion   [#cargas,4 columnas]%
%[kN/m ,kN/m,elemeto cargado,direccion X   ] 
%  para cargas proyectadas cos cargas verticales
%  para cargas proyectadas sen cargas horizontales
%qx=-(q2-q1)/L^nq*x.^nq-q1;
%bax=-(b2a-b1a)/L^naxi*x.^naxi-b1a;

elemento_carga=[0 ,-50,2,Y,nq,naxi;
                -50,0  ,3,Y,nq,naxi];
%elemento_carga=[0,0,1,Y];            
%% carga puntual 
% P,nudo,direccion X Y G

nudo_carga=[ -40,5,Y;
             20,5,G];
%nudo_carga=[ 0,1,Y];
                          
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XYG,1;
            XYG,6];

[qe_glob,ad,Kdd,Mdd,d,xye,elementos,Le,GLe,nudos]=porticos(EG,mat_elem,xyc,ninf,sec,elemento_carga,tipo_apoyo,tf,tw,bh,nudo_carga);


