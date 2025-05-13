function viga(
%clc
%clear
%close all
% se definen algunas constantes que hacen el codigo mas legible
X = 1;
Y = 2;
G = 3;
XY = 12;
XYG=123;

elementos=80;
nudos=elementos+1;
L=4;
x=linspace(0,L,nudos);
ln=L/elementos;
%% Se describen las propiedades de los materiales
%%% E G
Es=24870062;
v=0.25;
Gs=Es/(2*(1+v));
EG=[Es,Gs;%[2,#materiales]
    Es,Gs];
%% tipo de material por elemento 
mat_elem = ones(1,elementos); %[1,#elemetos]
%% ingreso coordenadas de nudos de estructura % [#coordenadas,2]

xyc= [linspace(0,L,nudos)',zeros(nudos,1)];

%% Nudos que conectan elementos
% x y %% [#nuros,2]

ninf=[(1:elementos)',(2:(elementos+1))'];

%%  Propiedades de la geometria secciones 
b1=0.3;     % ancho inicial viga
b2=0.2;     % ancho final viga
nb=1;       % exponente ancho final viga

h1=0.5;     % altura inicial viga
h2=0.8;     % altura final viga
nh=1;       % exponente altura viga

L=4;        % longitud viga
q1=25;      % carga vertical inicial viga
q2=30;      % carga vertical final viga
nq=1;       % exponente carga vertical final viga

b1a=20;     % carga axial inicial viga
b2a=25;     % carga axial final viga
naxi=1;     % exponenete carga axial final viga

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

nbe=ones(elementos,1)*nb;
nhe=ones(elementos,1)*nh;

bx=abs(-(b2-b1)/L^nb*x.^nb-b1);
hx=abs(-(h2-h1)/L^nh*x.^nh-h1);

bh=[bx(1:elementos)',bx(2:(elementos+1))',hx(1:elementos)',hx(2:(elementos+1))',nbe,nhe];
%% cargas aplicadas distribuida (gdl carga)
% wix   , wfx , direccion   [#cargas,4 columnas]%
%[kN/m ,kN/m,elemeto cargado,direccion X   ] 
%  para cargas proyectadas cos cargas verticales
%  para cargas proyectadas sen cargas horizontales

qx =-(q2 -q1 )/L^nq  *x.^nq  -q1;
bax=-(b2a-b1a)/L^naxi*x.^naxi-b1a;

y=ones(elementos,1)*Y;
x=ones(elementos,1)*X;
nqe=ones(elementos,1)*nq;
naxi=ones(elementos,1)*naxi;
elemento_carga=[qx(1:elementos)' ,qx(2:(elementos+1))' ,(1:elementos)',y,nqe,naxi;
                bax(1:elementos)',bax(2:(elementos+1))',(1:elementos)',x,nqe,naxi];
%% carga puntual 
% P,nudo,direccion X Y G

%nudo_carga=[ -30,2,Y;
%             -20,3,X
%             35,3,G];
nudo_carga=[ 0,1,Y]; 
                          
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XYG,1;
            XYG,nudos];
[qe_glob,ad]=porticos(EG,mat_elem,xyc,ninf,sec,elemento_carga,tipo_apoyo,tf,tw,bh,nudo_carga);


