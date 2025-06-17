clc
clear
close all
datos=[1,4,1.10000000000000,1.10000000000000,5.50000000000000,1.10000000000000,0.300000000000000,0.410000000000000;2,4,1.15000000000000,1.15000000000000,5.60000000000000,1.15000000000000,0.300000000000000,0.350000000000000;1,4,1.20000000000000,1.20000000000000,5.70000000000000,1.20000000000000,0.300000000000000,0.430000000000000;1,5,1.25000000000000,1.25000000000000,5.80000000000000,1.25000000000000,0.400000000000000,0.520000000000000;2,5,1.30000000000000,1.30000000000000,5.90000000000000,1.30000000000000,0.400000000000000,0.440000000000000;1,5,1.35000000000000,1.35000000000000,6,1.35000000000000,0.400000000000000,0.540000000000000;2,4,1.40000000000000,1.40000000000000,6.10000000000000,1.40000000000000,0.300000000000000,0.330000000000000;1,5,1.45000000000000,1.45000000000000,6.20000000000000,1.45000000000000,0.400000000000000,0.570000000000000;1,4,1.50000000000000,1.50000000000000,6.30000000000000,1.50000000000000,0.350000000000000,0.490000000000000;2,5,1.55000000000000,1.55000000000000,6.40000000000000,1.55000000000000,0.400000000000000,0.500000000000000;2,4,1.60000000000000,1.60000000000000,6.50000000000000,1.60000000000000,0.400000000000000,0.410000000000000;1,4,1.65000000000000,1.65000000000000,6.60000000000000,1.65000000000000,0.400000000000000,0.520000000000000;2,5,1.70000000000000,1.70000000000000,6.70000000000000,1.70000000000000,0.400000000000000,0.530000000000000;1,4,1.75000000000000,1.75000000000000,6.80000000000000,1.75000000000000,0.400000000000000,0.550000000000000;2,5,1.80000000000000,1.80000000000000,6.90000000000000,1.80000000000000,0.500000000000000,0.550000000000000];

No=2;
n=datos(No,3);%
L1=datos(No,4);%
L2=datos(No,5);%
L3=datos(No,6);%
b=datos(No,7);
h=datos(No,8);

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
mat_elem = [1,1,1,1];
    %[1,#elemetos]
%% ingreso coordenadas de nudos de estructura % [#coordenadas,2]
xyc= [0,0;
      L1,0;
      L1+L2/2,0;
      L1+L2,0;
      L1+L2+L3,0];

%% Nudos que conectan elementos
% x y %% [#nuros,2]
ninf=[1,2;
      2,3;
      3,4;
      4,5];
%%  Propiedades de la geometria secciones 
%1Seccion C
%2Seccion rectangular 
%3Circular solida 
%4Seccion I 
%5Seccion circular hueca
%6Tubular hueca 
%7Perfil L 
%8Perfil T

%% vigas
dv = h;
bfv = b;
twv =0;
tfv= 0;
[Aiv,Ixiv,As2iv] = propiedades_geometrica_perfiles_v1(dv,bfv,twv,tfv);
sec1=2;
sec2=2;

% [m^2 area,% m^4 inercia_y,% m^2 area cortante];
%A  I AC  [3 columnas,numero de secciones]
sec = [Aiv(sec2),Ixiv(sec2),As2iv(sec2)];% seccion 2 vigas
   
%% tipo de sección para cada elemento
sec_ele = [1;1;1;1];
%% cargas aplicadas distribuida (gdl carga)
% wix   , wfx , direccion   [#cargas,4 columnas]%
%[kN/m ,kN/m,elemeto cargado,direccion X   ] 
%  para cargas proyectadas cos cargas verticales
%  para cargas proyectadas sen cargas horizontales
elemento_carga=[-38*n,-65*n,1,Y;
                -65*n,-95*n,2,Y;
                -95*n,-80*n,3,Y;
                -80*n,-88*n,4,Y];
%% carga puntual 
% P,nudo,direccion X Y G
nudo_carga=[ -95*n,1,Y;
             -150*n,2,Y;
             -250*n,3,Y;
             -280*n,4,Y;
             -275*n,5,Y
              255*n,1,G;
              235*n,2,G;
             -285*n,3,G;
             -270*n,4,G;
              228*n,5,G;];
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XY,datos(No,1);
            XY,datos(No,2)]; 
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
[qe_glob,ad,Kdd,Mdd,GLe,d,elementos,nudos,xye,MMexdd,T,Ke,ang,AEe,EIe,Ace,Le,P,iteraciones,ite,puntos_graficas]=porticos(EG,mat_elem,xyc,ninf,sec,sec_ele,elemento_carga,nudo_carga,tipo_apoyo,tipo_conti,nudo_resorte,kw,lalb,PD,esc_def,esc_faxial,esc_V,esc_M,fesct,fesc);