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
mat_elem = ones(1,15); %[1,#elemetos]
%% ingreso coordenadas de nudos de estructura % [#coordenadas,2]

xyc= [0,0;0,3.50000000000000;0,6.50000000000000;0,10.3000000000000;6,0;6,3.50000000000000;6,6.50000000000000;6,10.3000000000000;12.5000000000000,0;12.5000000000000,3.50000000000000;18.5000000000000,0;18.5000000000000,3.50000000000000;18.5000000000000,6.50000000000000];


%% Nudos que conectan elementos
% x y %% [#nuros,2]
ninf=[1,2;2,3;3,4;5,6;6,7;7,8;9,10;11,12;12,13;2,6;6,10;10,12;3,7;7,13;4,8];
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
dc = 0.6;% altura viga
bfc = 0.4;% ancho
tw = 9/1000;
tf=  4/1000;
[Aic,Ixic,As2ic] = propiedades_geometrica_perfiles_v1(dc,bfc,tw,tf);

%% columnas
dc = 0.6;% altura viga
bfc = 0.4;% ancho
tw = 9/1000;
tf=  4/1000;
[A2c,Ix2c,As2i2c] = propiedades_geometrica_perfiles_v1(dc,bfc,tw,tf);
%% vigas
dv = 0.7; %altura
bfv = 0.3;% ancho
tw = 9/1000;
tf=  4/1000;
[Aiv,Ixiv,As2iv] = propiedades_geometrica_perfiles_v1(dv,bfv,tw,tf);

seci=2;
sect=6;
% [m^2 area,% m^4 inercia_y,% m^2 area cortante];
%A  I AC  [3 columnas,numero de secciones]
sec = [Aic(seci),Ixic(seci),As2ic(seci); % seccion 1  columnas
       Aiv(seci),Ixiv(seci),As2iv(seci);% seccion 2  vigas
       A2c(seci),Ix2c(seci),As2i2c(seci)];% seccion 3  columnas
   
%% tipo de sección para cada elemento
sec_ele = [1;1;1;1;1;1;1;1;1;2;2;2;2;2;2];
           
%% cargas aplicadas distribuida (gdl carga)
% wix   , wfx , direccion   [#cargas,4 columnas]%
%[kN/m ,kN/m,elemeto cargado,direccion X   ] 
%  para cargas proyectadas cos cargas verticales
%  para cargas proyectadas sen cargas horizontales

elemento_carga=[-50  , -50,10 ,Y;
                -50  , -50,11 ,Y;
                -50  , -50,12 ,Y;
                -50  , -50,13 ,Y;
                -50  , -50,14 ,Y;
                -50  , -50,15 ,Y];
%% carga puntual 
% P,nudo,direccion X Y G

nudo_carga=[ 0,1,X];
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XYG,1;
            XYG,5;
            XYG,9;
            XYG,11];
           
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
            'EE';%14
            'EE'];%24
%% coeficioente de balastro para elemetos apoyado sobre lecho elastico        
kw= zeros(15,1);
PD='no';
%% longitud de nudo rigido
lalb=[0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
%lalb=[0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2];
%% escalas de Dibujo la estructura y su deformada
esc_def    = 500;          % escalamiento de la deformada
esc_faxial = 0.006;      % escalamiento del diagrama de axiales
esc_V      = 0.01;       % escalamiento del diagrama de cortantes
esc_M      = 0.01;      % escalamiento del diagrama de momentos
fesct=max(max(xyc))*11/16*3;%% escala textos
fesc=max(max(xyc))*0.4/16;%% escala apoyos
[qe_glob,ad,Kdd,Mdd,GLe,d,elementos,nudos,xye,MMexdd,T,Ke,ang,AEe,EIe,Ace,Le,P,iteraciones,ite,puntos_graficas]=porticos(EG,mat_elem,xyc,ninf,sec,sec_ele,elemento_carga,nudo_carga,tipo_apoyo,tipo_conti,nudo_resorte,kw,lalb,PD,esc_def,esc_faxial,esc_V,esc_M,fesct,fesc);