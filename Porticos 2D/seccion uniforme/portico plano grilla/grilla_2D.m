clear
clc
%% constantes que ayudarán en la lectura del código
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
gama=24;
v=0.25;
Gs=Es/(2*(1+v));
EG=[Es,Gs];%[2,#materiales]];

%% ingreso coordenadas de nudos de estructura % [#coordenadas,2]

x=[6];
y=[5,6,9];
xa= [0,cumsum(x,2)];
ya= [0,cumsum(y,2)];
[Xe,Ye]=meshgrid(xa,ya);%% grilla coordenadas nudos elemento
Bx = reshape(Xe,[],1);
By = reshape(Ye,[],1);
Nx=size(xa,2);
Ny=size(ya,2);
%% nudos elemento
nno  = Nx*Ny; % numero de nodos (numero de filas de xnod)
nnoi=1:nno;
xnod=[Bx,By];

%% nugos grilla elemento
Nudosg = reshape(nnoi,[Ny,Nx]);

nef=(Ny-1)*Nx+(Nx-1)*(Ny-1);
recorridoy=2:1:Ny;
recorridox=2:1:Nx;

LaGc = zeros(nef,2); 
xeg = zeros(nef,2); 
yeg = zeros(nef,2); 
recorridoLaG=1:1:nef;
cg1 = zeros(nef,2); % almacena el centro de gravedad de los EFs
hold on
text(Bx, By, num2str(nnoi'), 'Color', 'k');
e=1;
%% nudos columnas
for ey=1:(Nx)
    for ex=1:(Ny-1)
       %%nudos por elemento
       LaGc(e,:)=Nudosg((recorridoy(ex)-1):recorridoy(ex),ey)'; 
       %coordenadas elemento
       xeg(e,:)=Xe((recorridoy(ex)-1):recorridoy(ex),ey)'; 
       yeg(e,:)=Ye((recorridoy(ex)-1):recorridoy(ex),ey)'; 
       cg1(e,:) = [mean(xeg(e,:)),mean(yeg(e,:))];
       text(cg1(e,X), cg1(e,Y), num2str(recorridoLaG(e)), 'Color', 'b');
       plot(xeg(e,:),yeg(e,:))
       e=e+1;
    end
end
%% nudso vigas
%% nudos columnas
for ey=1:(Nx-1)
    for ex=2:(Ny)
       %%nudos por elemento
       LaGc(e,:)=Nudosg(ex,(recorridox(ey)-1):recorridox(ey)); 
       %coordenadas elemento
       xeg(e,:)=Xe(ex,(recorridox(ey)-1):recorridox(ey)); 
       yeg(e,:)=Ye(ex,(recorridox(ey)-1):recorridox(ey)); 
       cg1(e,:) = [mean(xeg(e,:)),mean(yeg(e,:))];
       text(cg1(e,X), cg1(e,Y), num2str(recorridoLaG(e)), 'Color', 'b');
       plot(xeg(e,:),yeg(e,:))
       e=e+1;
    end

end

axis equal tight
%grid on
title('Malla de elementos finitos');
plot(Xe,Ye, '*r'),hold on%% grafica nodos elemento finito

%% Nudos que conectan elementos
ninf=LaGc;
%% ingreso coordenadas de nudos de estructura % [#coordenadas,2]
xyc=xnod;

columnas=(Ny-1)*Nx;
vigas=(Nx-1)*(Ny-1);
% dimensiones columnas
dc=0.5;
bfc=0.4;
tw=0;
tf=0;
[Aic,Ixic,As2ic] = propiedades_geometrica_perfiles_v1(dc,bfc,tw,tf);
% dimensiones vigas
dv=0.5;
bfv=0.4;

[Aiv,Ixiv,As2iv] = propiedades_geometrica_perfiles_v1(dv,bfv,tw,tf);
sec1=2;

% [m^2 area,% m^4 inercia_y,% m^2 area cortante];
%A  I AC  [3 columnas,numero de secciones]
sec = [Aic(sec1),Ixic(sec1),As2ic(sec1); % seccion 1  columnas
       Aiv(sec1),Ixiv(sec1),As2iv(sec1)];% seccion 2    vigas
%% tipo de secci�n para cada elemento
sec_ele = [ones(columnas,1);ones(vigas,1)*2];

%% cargas aplicadas distribuida (gdl carga)
% wix   , wfx , direccion   [#cargas,4 columnas]%
%[kN/m ,kN/m,elemeto cargado,direccion X   ] 
%  para cargas proyectadas cos cargas verticales
%  para cargas proyectadas sen cargas horizontales

wic=Aic(sec1,1)*gama;
wiv=Aiv(sec1,1)*gama;
wc=ones(columnas,1)*wic;
wv=ones(vigas,1)*wiv;
w=[wc;wv];
nele=(1:nef)';

elemento_carga=[-w  ,-w,nele ,ones(nef,1)*2];

%% carga puntual 
% P,nudo,direccion X Y G

nudo_carga=[ 0,2,Y;
             0,6,Y;
             0,3,Y;
             0,7,Y;
             20,4,X;
             20,8,X];
         
%% tipo de apoyo 
% retriccion,nudo
%%X horizontal, Y vertical, XY horozonal y vertical, XYG horizontal vertica y giro
tipo_apoyo=[XYG,1;
            XYG,5];
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
tipo_conti=char(ones(columnas+vigas,2).*[69,69]);%% coeficioente de balastro para elemetos apoyado sobre lecho elastico        
%% coeficioente de balastro para elemetos apoyado sobre lecho elastico        
kw= zeros(24,1);
PD='no';
%% longitud de nudo rigido
lalb=zeros(columnas+vigas,2);
%% tipo de material por elemento 
mat_elem = ones(1,columnas+vigas); %[1,#elemetos]
%% escalas de Dibujo la estructura y su deformada
esc_def    = 200;          % escalamiento de la deformada
esc_faxial = 0.01;      % escalamiento del diagrama de axiales
esc_V      = 0.02;       % escalamiento del diagrama de cortantes
esc_M      = 0.03;      % escalamiento del diagrama de momentos
fesct=max(max(xyc))*11/16*1;%% escala textos
fesc=max(max(xyc))*0.4/16;%% escala apoyos
[qe_glob,ad,Kdd,Mdd,GLe,d,elementos,nudos,xye,MMexdd,T,Ke,ang,AEe,EIe,Ace,Le,P,iteraciones,ite,puntos_graficas]=porticos(EG,mat_elem,xyc,ninf,sec,sec_ele,elemento_carga,nudo_carga,tipo_apoyo,tipo_conti,nudo_resorte,kw,lalb,PD,esc_def,esc_faxial,esc_V,esc_M,fesct,fesc);     
        
           