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
Es=24870062;
v=0.25;
Gs=Es/(2*(1+v));
EG=[Es,Gs;%[2,#materiales]
    Es,Gs];
%% tipo de material por elemento 
mat_elem = ones(1,5); %[1,#elemetos]
%% ingreso coordenadas de nudos de estructura % [#coordenadas,2]

xyc= [0,0;
      4,3;
      6.5,6;
      10,0];

%% Nudos que conectan elementos
% x y %% [#nuros,2]

ninf=[1,2;
      2,3;
      3,4];

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
bh=[0.3,0.2,0.5,0.8,nb,nh;
    0.2,0.3,0.8,0.5,nb,nh;
    0.3,0.2,0.5,0.8,nb,nh];
%% cargas aplicadas distribuida (gdl carga)
% wix   , wfx , direccion   [#cargas,4 columnas]%
%[kN/m ,kN/m,elemeto cargado,direccion X   ] 
%  para cargas proyectadas cos cargas verticales
%  para cargas proyectadas sen cargas horizontales

%qx=-(q2-q1)/L^nq*x.^nq-q1;
%bax=-(b2a-b1a)/L^naxi*x.^naxi-b1a;

elemento_carga=[-25,-30,1,Y,nq,naxi;
                -28,-35,2,X,nq,naxi;
                -35,-22,3,Y,nq,naxi];
%elemento_carga=[0,0,1,Y];            
%% carga puntual 
% P,nudo,direccion X Y G

nudo_carga=[ -30,2,Y;
             -20,3,X
             35,3,G];
%nudo_carga=[ 0,1,Y];
                          
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XYG,1;
            XYG,4];

[qe_glob,ad]=porticos(EG,mat_elem,xyc,ninf,sec,elemento_carga,tipo_apoyo,tf,tw,bh,nudo_carga);


