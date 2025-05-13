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
Es=25*10^6;
v=0.25;
Gs=Es/(2*(1+v));
EG=[Es,Gs;%[2,#materiales]
    Es,Gs];
%% tipo de material por elemento 
mat_elem = ones(1,24); %[1,#elemetos]
%% ingreso coordenadas de nudos de estructura % [#coordenadas,2]

xyc= [0,0;
      0,4;
      0,4+3;
      0,4+6;
      0,4+9;
      4,0;
      4,4;
      4,4+3;
      4,4+6;
      4,4+9;
      7,0;
      7,4;
      7,4+3;
      7,4+6;
      11,0;
      11,4;
      11,7];

%% Nudos que conectan elementos
% x y %% [#nuros,2]
ninf=[1,2;
      2,3;
      3,4;
      4,5;
      6,7;
      7,8;
      8,9;
      9,10;
      11,12;
      12,13;
      13,14;
      15,16;
      16,17;
      2,7;
      7,12;
      12,16
      3,8;
      8,13;
      13,17;
      4,9;
      9,14;
      5,10;
      15,12;
      12,17];
%%  Propiedades de la geometria secciones 
%1Seccion C
%2Seccion rectangular 
%3Circular solida 
%4Seccion I 
%5Seccion circular hueca
%6Tubular hueca 
%7Perfil L 
%8Perfil T
%% columnas
dc = 0.4;
bfc = 0.4;
tw = 9/1000;
tf=  4/1000;
[Aic,Ixic,As2ic] = propiedades_geometrica_perfiles_v1(dc,bfc,tw,tf);
%% vigas
dv = 0.6;
bfv = 0.4;
tw = 9/1000;
tf=  4/1000;
[Aiv,Ixiv,As2iv] = propiedades_geometrica_perfiles_v1(dv,bfv,tw,tf);

%% muros
dm = 2;
bfm = 0.2;
tw = 9/1000;
tf=  4/1000;
[Aim,Ixim,As2im] = propiedades_geometrica_perfiles_v1(dm,bfm,tw,tf);

%diagonales
dd = 0.1;
bfd = 0.1;
tw = 4/1000;
tf=  4/1000;
[Aid,Ixid,As2id] = propiedades_geometrica_perfiles_v1(dd,bfd,tw,tf);

sec=2;
sect=6;
% [m^2 area,% m^4 inercia_y,% m^2 area cortante];
%A  I AC  [3 columnas,numero de secciones]
sec = [Aic(sec),Ixic(sec),As2ic(sec); % seccion 1  columnas
       Aiv(sec),Ixiv(sec),As2iv(sec);% seccion 2    vigas
       Aim(sec),Ixim(sec),As2im(sec);% seccion 3   muros
       Aid(sect),Ixid(sect),As2id(sect)];% seccion 4  diagonal
%% tipo de sección para cada elemento
sec_ele = [3;3;3;3;
           1;1;1;1;
           1;1;1;
           1;1;
           2;2;2;
           2;2;2;
           2;2;
           2
           4;4];
           
%% cargas aplicadas distribuida (gdl carga)
% wix   , wfx , direccion   [#cargas,4 columnas]%
%[kN/m ,kN/m,elemeto cargado,direccion X   ] 
%  para cargas proyectadas cos cargas verticales
%  para cargas proyectadas sen cargas horizontales

elemento_carga=[-25  , -25,14 ,Y;
                -25  , -25,15 ,Y;
                -25  , -25,16 ,Y;
                -25  , -25,17 ,Y;
                -25  , -25,18 ,Y;
                -25  , -25,19 ,Y;
                -25  , -25,20 ,Y;
                -25  , -25,21 ,Y;
                -25  , -25,22 ,Y];
%% carga puntual 
% P,nudo,direccion X Y G

nudo_carga=[ 5,2,X;
             5,3,X;
             5,4,X;
             5,5,X];
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XYG,1;
            XYG,6;
            XYG,11;
            XYG,15];
           
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
            'EE';%3
            'EE';%4
            'EE';%5
            'EE';%6
            'EE';%7
            'EE';%8
            'EE';%9
            'EE';%10
            'EE';%11
            'EE';%12
            'EE';%13
            'NR';%14
            'EE';%15
            'EE';%16
            'NR';%17
            'EE';%18
            'EE';%19
            'NR';%20
            'EE';%21
            'NR';%22
            'RR';%23
            'RR'];%24
%% coeficioente de balastro para elemetos apoyado sobre lecho elastico        
kw= zeros(24,1);
PD='no';
%% longitud de nudo rigido
lalb=[0,0;
      0,0;
      0,0;
      0,0;
      0,0;
      0,0;
      0,0;
      0,0;
      0,0;
      0,0;
      0,0;
      0,0;
      0,0;
      1,0;
      0,0;
      0,0;
      1,0;
      0,0;
      0,0;
      1,0;
      0,0;
      1,0;
      0,0;
      0,0];
%lalb=[0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2];
%% escalas de Dibujo la estructura y su deformada
esc_def    = 500;          % escalamiento de la deformada
esc_faxial = 0.006;      % escalamiento del diagrama de axiales
esc_V      = 0.01;       % escalamiento del diagrama de cortantes
esc_M      = 0.01;      % escalamiento del diagrama de momentos
fesct=max(max(xyc))*11/16*3;%% escala textos
fesc=max(max(xyc))*0.4/16;%% escala apoyos
[qe_glob,ad,Kdd,Mdd,GLe,d,elementos,nudos,xye,MMexdd,T,Ke,ang,AEe,EIe,Ace,Le,P,iteraciones,ite,puntos_graficas]=porticos(EG,mat_elem,xyc,ninf,sec,sec_ele,elemento_carga,nudo_carga,tipo_apoyo,tipo_conti,nudo_resorte,kw,lalb,PD,esc_def,esc_faxial,esc_V,esc_M,fesct,fesc);