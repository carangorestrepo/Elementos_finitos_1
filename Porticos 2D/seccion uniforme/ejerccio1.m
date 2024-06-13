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

xyc= [0,13;
      4,13;
      0,10;
      4,10;
      7,10;
      0,7;
      4,7;
      7,7;
      11,7;
      0,4;
      4,4;
      7,4;
      11,4;
      0,0;
      4,0;
      7,0;
      11,0];

%% Nudos que conectan elementos
% x y %% [#nuros,2]
ninf=[1,2;
      3,4;
      4,5;
      6,7;
      7,8;
      8,9;
      10,11;
      11,12;
      12,13;
      1,3;
      2,4;
      3,6;
      4,7;
      5,8;
      6,10;
      7,11;
      8,12;
      9,12;
      9,13;
      10,14;
      11,15;
      12,16;
      12,17;
      13,17];
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
sect=3;
% [m^2 area,% m^4 inercia_y,% m^2 area cortante];
%A  I AC  [3 columnas,numero de secciones]
sec = [Aic(sec),Ixic(sec),As2ic(sec); % seccion 1  columnas
       Aiv(sec),Ixiv(sec),As2iv(sec);% seccion 2    vigas
       Aim(sec),Ixim(sec),As2im(sec);% seccion 3   muros
       Aid(sect),Ixid(sect),As2id(sect)];% seccion 4  diagonal
%% tipo de sección para cada elemento
sec_ele = [2;2;2;2;2;2;2;2;2;
           3;
           1;
           3;
           1;1;
           3;
           1;1;
           4;
           1;
           3;
           1;1;
           4;
           1];
           
%% cargas aplicadas distribuida (gdl carga)
% wix   , wfx , direccion   [#cargas,4 columnas]%
%[kN/m ,kN/m,elemeto cargado,direccion X   ] 
%  para cargas proyectadas cos cargas verticales
%  para cargas proyectadas sen cargas horizontales

elemento_carga=[0  , 0,1 ,Y;
                0  , 0,2 ,Y;
                0  , 0,3 ,Y;
                0  , 0,4 ,Y;
                0  , 0,5 ,Y;
                0  , 0,6 ,Y;
                0  , 0,7 ,Y;
                0  , 0,8 ,Y;
                0  , 0,9 ,Y];
%% carga puntual 
% P,nudo,direccion X Y G

nudo_carga=[  607.147519582245,10,X;
              1120.887728459530 ,6,X;
              1176.932114882510 ,3,X;
              672.532637075718 ,1,X];
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XYG,14;
            XYG,15;
            XYG,16;
            XYG,17];
           
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
tipo_conti=['NR';%1
            'NR';%2
            'EE';%3
            'NR';%4
            'EE';%5
            'EE';%6
            'NR';%7
            'EE';%8
            'EE';%9
            'EE';%10
            'EE';%11
            'EE';%12
            'EE';%13
            'EE';%14
            'EE';%15
            'EE';%16
            'EE';%17
            'EE';%18
            'EE';%19
            'EE';%20
            'EE';%21
            'EE';%22
            'EE';%23
            'EE'];%24
%% coeficioente de balastro para elemetos apoyado sobre lecho elastico        
kw= zeros(24,1);
PD='no';
%% longitud de nudo rigido
lalb=[1,0;
      1,0;
      0,0;
      1,0;
      0,0;
      0,0;
      1,0;
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
      0,0;
      0,0;
      0,0;
      0,0;
      0,0];
%lalb=[0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2];
%% escalas de Dibujo la estructura y su deformada
esc_def    = 100;          % escalamiento de la deformada
esc_faxial = 0.001;      % escalamiento del diagrama de axiales
esc_V      = 0.0005;       % escalamiento del diagrama de cortantes
esc_M      = 0.0008;      % escalamiento del diagrama de momentos
fesct=max(max(xyc))*11/16*3;%% escala textos
fesc=max(max(xyc))*0.4/16;%% escala apoyos
[qe_glob,ad]=porticos(EG,mat_elem,xyc,ninf,sec,sec_ele,elemento_carga,nudo_carga,tipo_apoyo,tipo_conti,nudo_resorte,kw,lalb,PD,esc_def,esc_faxial,esc_V,esc_M,fesct,fesc);