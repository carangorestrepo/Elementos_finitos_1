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
Es=24870062;%%4700*sqrt(28)*1000;
Es1=200000*1000;
v=0.25;
v1=0.3;
Gs=Es/(2*(1+v));
Gs1=Es/(2*(1+v1));
EG=[Es,Gs; % marerial 1 concreto %[2,#materiales]
    Es1,Gs1];% marerial 2 metalico
%% tipo de material por elemento 
mat_elem = [1,1,1,1,1,1,1,1,1,1];
    %[1,#elemetos]
%% ingreso coordenadas de nudos de estructura % [#coordenadas,2]
xyc= [0,0;0,3.5;
      0,6.5;
      6.5,0;
      6.5,3.5;
      6.5,6.5;
      13.25,0;
      13.25,3.5;
      14.75,3.5];

%% Nudos que conectan elementos
% x y %% [#nuros,2]
ninf=[1,2;
     2,3;
     4,5;
     5,6;
     7,8;
     2,5;
     5,8;
     3,6;
     8,6;
     8,9];
%%  Propiedades de la geometria secciones 
%1Seccion C
%2Seccion rectangular 
%3Circular solida 
%4Seccion I 
%5Seccion circular hueca
%6Tubular hueca 
%7Perfil L 
%8Perfil T
%% Columnas
dc = 0.6;
bfc = 0.4;
twc =0;
tfc= 0;
[Aic,Ixic,As2ic] = propiedades_geometrica_perfiles_v1(dc,bfc,twc,tfc);

%% vigas
dv = 0.45;
bfv = 0.4;
twv =0;
tfv= 0;
[Aiv,Ixiv,As2iv] = propiedades_geometrica_perfiles_v1(dv,bfv,twv,tfv);
sec1=2;
sec2=2;

% [m^2 area,% m^4 inercia_y,% m^2 area cortante];
%A  I AC  [3 columnas,numero de secciones]
sec = [Aic(sec1),Ixic(sec1),As2ic(sec1); % seccion 1 columnas
       Aiv(sec2),Ixiv(sec2),As2iv(sec2)];% seccion 2 vigas
   
%% tipo de sección para cada elemento
sec_ele = [1;1;1;1;1;2;2;2;1;2];
%% cargas aplicadas distribuida (gdl carga)
% wix   , wfx , direccion   [#cargas,4 columnas]%
%[kN/m ,kN/m,elemeto cargado,direccion X   ] 
%  para cargas proyectadas cos cargas verticales
%  para cargas proyectadas sen cargas horizontales
elemento_carga=[-27.4,-27.4,6 ,Y;
                -27.4,-27.4,7 ,Y;
                -27.4,-27.4,8 ,Y;
                -27.4,-27.4,10,Y];
%% carga puntual 
% P,nudo,direccion X Y G
nudo_carga=[ 0,1,X];
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XYG,1;
            XY,4;
            Y,7]; 
%% 
%% X resorte, resorte Y  resorte G
% reigidez resoporte,nudo,direccion X Y G
nudo_resorte=[0,1,Y];
              %0,3,Y];          
%% tipo de continuidad
%E empotrado
%R Rotula
% EE , ER, RE, RR
%% CON RESORTE
%kw
%% con nudos rigidos
% NR
tipo_conti=['EE';%1
            'EE';%2
            'EE';
            'EE';
            'EE';
            'EE';
            'EE';
            'EE';
            'EE';
            'EE'];%3
%% coeficioente de balastro para elemetos apoyado sobre lecho elastico        
kw= zeros(10,1);
PD='no';
%% longitud de nudo rigido
lalb=zeros(10,2);
%lalb=[0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2];
%% escalas de Dibujo la estructura y su deformada
esc_def    = 500;          % escalamiento de la deformada
esc_faxial = 0.003;      % escalamiento del diagrama de axiales
esc_V      = 0.02;       % escalamiento del diagrama de cortantes
esc_M      = 0.01;      % escalamiento del diagrama de momentos
fesct=max(max(xyc))*11/16*3;%% escala textos
fesc=max(max(xyc))*0.4/16;%% escala apoyos
[qe_glob,ad,Kdd,Mdd,GLe,d,elementos,nudos,xye,MMexdd,T,Ke,ang,AEe,EIe,Ace,Le,P,iteraciones,ite,puntos_graficas,gama,Ce]=porticos(EG,mat_elem,xyc,ninf,sec,sec_ele,elemento_carga,nudo_carga,tipo_apoyo,tipo_conti,nudo_resorte,kw,lalb,PD,esc_def,esc_faxial,esc_V,esc_M,fesct,fesc);