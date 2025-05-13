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
Es=24855578;
v=0.25;
Gs=Es/(2*(1+v));
EG=[Es,Gs;%[2,#materiales]
    Es,Gs];
%% tipo de material por elemento 
mat_elem = ones(1,3); %[1,#elemetos]
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
%1Seccion C
%2Seccion rectangular 
%3Circular solida 
%4Seccion I 
%5Seccion circular hueca
%6Tubular hueca 
%7Perfil L 
%8Perfil T
d = 0.5;
bf = 0.3;
tw = 9/1000;
tf=  4/1000;
[Ai,Ixi,As2i] = propiedades_geometrica_perfiles_v1(d,bf,tw,tf);
sec=2;
% [m^2 area,% m^4 inercia_y,% m^2 area cortante];
%A  I AC  [3 columnas,numero de secciones]
sec = [Ai(sec),Ixi(sec),As2i(sec); % seccion 1 As2i(sec);
       0.3*0.3,0.3*0.3^3/12,5/6*0.3*0.3];% seccion 2
%% tipo de sección para cada elemento
sec_ele = ones(5,1);
%% cargas aplicadas distribuida (gdl carga)
% wix   , wfx , direccion   [#cargas,4 columnas]%
%[kN/m ,kN/m,elemeto cargado,direccion X   ] 
%  para cargas proyectadas cos cargas verticales
%  para cargas proyectadas sen cargas horizontales



elemento_carga=[-25  , -30,1 ,Y;
                -28  , -35,2 ,X;
               -35   , -22,3 ,Y];
%% carga puntual 
% P,nudo,direccion X Y G

nudo_carga=[ -30,2,Y;
             -20,3,X;
              35,3,G];
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XYG,1;
            XYG,4];
           
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
            'EE'];%3
%% coeficioente de balastro para elemetos apoyado sobre lecho elastico        
kw= zeros(5,1);
PD='no';
%% longitud de nudo rigido
lalb=zeros(3,2);
%lalb=[0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0,0.4;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2;0.2,0.2];
%% escalas de Dibujo la estructura y su deformada
esc_def    = 500;          % escalamiento de la deformada
esc_faxial = 0.001;      % escalamiento del diagrama de axiales
esc_V      = 0.05;       % escalamiento del diagrama de cortantes
esc_M      = 0.01;      % escalamiento del diagrama de momentos
fesct=max(max(xyc))*11/16*3;%% escala textos
fesc=max(max(xyc))*0.4/16;%% escala apoyos
[qe_glob,ad]=porticos(EG,mat_elem,xyc,ninf,sec,sec_ele,elemento_carga,nudo_carga,tipo_apoyo,tipo_conti,nudo_resorte,kw,lalb,PD,esc_def,esc_faxial,esc_V,esc_M,fesct,fesc);