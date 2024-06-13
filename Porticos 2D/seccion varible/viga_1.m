function [u,v,M,V,fax,xx]=viga_1(Es,Gs,sec,q1,q2,b1a,b2a,L,nb,nh,nq,naxi,b1,b2,h1,h2,a,ang,x1,x2,y1,y2,w)
%% escalas de Dibujo la estructura y su deformada
esc_def    = 5000;          % escalamiento de la deformada
esc_faxial = 0.01;      % escalamiento del diagrama de axiales
esc_V      = 0.01;       % escalamiento del diagrama de cortantes
esc_M      = 0.01;      % escalamiento del diagrama de momentos
%clc
%clear
%close all
% se definen algunas constantes que hacen el codigo mas legible
X = 1;
Y = 2;
%G = 3;
%XY = 12;
XYG=123;

elementos=50;
nudos=elementos+1;
%L=4;
x=linspace(0,L,nudos);
%ln=L/elementos;
%% Se describen las propiedades de los materiales
%%% E G
%Es=24870062;
%v=0.25;

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
%b1=0.3;     % ancho inicial viga
%b2=0.2;     % ancho final viga
%nb=1;       % exponente ancho final viga

%h1=0.5;     % altura inicial viga
%h2=0.8;     % altura final viga
%nh=1;       % exponente altura viga

%L=4;        % longitud viga
%q1=25;      % carga vertical inicial viga
%q2=30;      % carga vertical final viga
%nq=1;       % exponente carga vertical final viga

%b1a=20;     % carga axial inicial viga
%b2a=25;     % carga axial final viga
%naxi=1;     % exponenete carga axial final viga

%sec=2;
%seccion C ==1
%Seccion rectangular  ==2
%Circular solida ==3 
%Seccion I ==4
%Seccion circular hueca ==5
%Tubular rectangular hueca ==6
%Perfil T==7
tf=0.09;
tw=0.06;
b1a=w(1);
q1=w(2);
b2a=w(4);
q2=w(5);

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

qx =(q2 -q1 )/L^nq  *x.^nq  +q1;
bax=(b2a-b1a)/L^naxi*x.^naxi+b1a;

y=ones(elementos,1)*Y;
x=ones(elementos,1)*X;
nqe=ones(elementos,1)*nq;
naxie=ones(elementos,1)*naxi;
elemento_carga=[qx(1:elementos)' ,qx(2:(elementos+1))' ,(1:elementos)',y,nqe,naxie;
                bax(1:elementos)',bax(2:(elementos+1))',(1:elementos)',x,nqe,naxie];
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
[u,v,M,V,fax,xx]=porticos_viga(EG,mat_elem,xyc,ninf,sec,elemento_carga,tipo_apoyo,tf,tw,bh,nudo_carga,a);

u=u';
v=v';
M=M';
V=V';
fax=fax';
xx=xx';
xs=xx;
s=xx;

daA = a(1);
da = a(2);
dbA = a(4);
db = a(5);


%v=double(subs(v,x,xs));
%u=double(subs(u,x,xs));
axial=fax;
%V=double(subs(V,x,xs));
%M=double(subs(M,x,xs));
figure(2)
xyv=[0,da;xs',v';L,db];
xyu=[0,daA;xs',u';L,dbA];
xs=[0,xs,L];
% rotacion de la solucion antes de dibujar
T   = [ cos(ang)  -sin(ang)    % matriz de rotacion
        sin(ang)   cos(ang) ];
%% Dibujar de deformada
pos = T*[ xs + esc_def*xyu(:,2)'; esc_def*xyv(:,2)'];
xx = pos(1,:) + x1;
yy = pos(2,:) + y1;
plot([x1 x2], [y1 y2], 'b-', xx, yy, 'r-','LineWidth',2),hold on

%% Dibujar los diagramas de fuerza axial 
figure(3)
pos = T*[ s; esc_faxial*axial ]; % escalamiento del diagrama

ss = pos(1,:) + x1;
aa = pos(2,:) + y1;

plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 aa y2], 'r-','LineWidth',2),hold on;
text(ss(1),   aa(1),   num2str(axial(1,1)));
text(ss(end), aa(end), num2str(axial(1,end)));

%% Dibujar los diagramas de fuerza cortante

pos = T*[ s; esc_V*V ]; % escalamiento del diagrama
ss = pos(1,:) + x1;
vv = pos(2,:) + y1;
figure(4)
plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 vv y2], 'r-','LineWidth',2),hold on;
text(ss(1),   vv(1),   num2str(V(1,1)));
text(ss(end), vv(end), num2str(V(1,end)));
%% Dibujar los diagramas de momento flector


pos = T*[ s; -esc_M*M ]; % escalamiento del diagrama
ss = pos(1,:) + x1;
mm = pos(2,:) + y1;
figure(5)
plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 mm y2], 'r-','LineWidth',2),hold on;
text(ss(1),   mm(1),   num2str(M(1,1)));
text(ss(end), mm(end), num2str(M(1,end)));
[minM,idminM] = min(M); text(ss(idminM), mm(idminM), num2str(minM));
[maxM,idmaxM] = max(M); text(ss(idmaxM), mm(idmaxM), num2str(maxM));

