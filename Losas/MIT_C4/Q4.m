clear all
close
clc
%% constantes que ayudarán en la lectura del código
X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo

%%
E=4700*sqrt(28)*1000;

nu=0.3;
h=4*0.1;
q  = -(4.6*1.2+1.8*1.6+0.2*24);         % [kN/m^2] carga
escala = 100; % factor de escalamiento de la deformada

Lx=4;
Ly=4;
deltax=0.1;
deltay=0.1;
%%123 empotrada 12%apoyado 0 libre  
EExi=123;
EExf=123;
EEyi=123;
EEyf=123;
kWinkler=500;
n=2;
rho=2.4028;

xqi=0;
xqf=Lx;
yqi=0;
yqf=Ly;

Nx=round(Lx/deltax,0);
dv=round(Nx/n,0);
Nx=dv*n+1;%% numero de nudos X

Ny=round(Ly/deltay,0);
dv=round(Ny/n,0);
Ny=dv*n+1;%% numero de nudos Y

x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);
[Xe,Ye]=meshgrid(x,y);%% grilla coordenadas nudos elemento

plot(Xe,Ye, '*r'),hold on%% grafica nodos elemento finito
Bx = reshape(Xe,[],1);
By = reshape(Ye,[],1);
xnod=[Bx,By];
%% nudos elemento
nnoi =(1:Ny*Nx)';% numero de nodos
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
%text(Bx,By,num2str(nnoi)),hold on
%% nugos grilla elemento
Nudosg = reshape(nnoi,[Ny,Nx]);

recorridox=2:1:Nx;
Nelex=size(recorridox,2);
recorridoy=2:1:Ny;
Neley=size(recorridoy,2);
nef=Nelex*Neley;% numero de EFs (numero de filas de LaG)
LaG=zeros(nef,4);
e=1;
LaGc = cell(nef,1); 
xeg = cell(nef,1); 
yeg = cell(nef,1); 
xe =zeros(nef,4);
ye =zeros(nef,4);
cg = zeros(nef,2); % almacena el centro de gravedad de los EFs
for ey=1:(Neley)
    for ex=1:(Nelex)
       %%nudos por elemento
       LaGc{e}=Nudosg((recorridoy(ey)-1):recorridoy(ey),(recorridox(ex)-1):recorridox(ex)); 
       %coordenadas elemento
       xeg{e}=Xe((recorridoy(ey)-1):recorridoy(ey),(recorridox(ex)-1):recorridox(ex)); 
       yeg{e}=Ye((recorridoy(ey)-1):recorridoy(ey),(recorridox(ex)-1):recorridox(ex)); 
       % se determinan las coordenadas de los nodos el EF e
       LaG(e,:)=[LaGc{e}(1,1),LaGc{e}(1,2),LaGc{e}(2,2),LaGc{e}(2,1)];
       xe(e,:)=[xeg{e}(1,1),xeg{e}(1,2),xeg{e}(2,2),xeg{e}(2,1)];
       ye(e,:)=[yeg{e}(1,1),yeg{e}(1,2),yeg{e}(2,2),yeg{e}(2,1)];
       % se calcula la posición del centro de gravedad del EF e
       cg(e,:) = [mean(xe(e,:)),mean(ye(e,:))];
       % se escribe el numero del EF e
       text(cg(e,X), cg(e,Y), num2str(e), 'Color', 'b');
       plot(xe(e,[1:4,1]),ye(e,[1:4,1]))
       e=e+1;
    end
end
axis equal tight
title('Malla de elementos finitos');
nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)


%% grados de libertad del desplazamiento conocidos y desconocidos
% determino los grados de libertad correspondientes a los bordes
lado_x0 = find(xnod(:,X) == 0);     
lado_y0 = find(xnod(:,Y) == 0);
lado_xLx = find(xnod(:,X) == Lx);     
lado_yLy = find(xnod(:,Y) == Ly);

% Definimos la geometria de la losa
%losa
% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad

%c = [ gdl(lado_x0,ww); gdl(lado_x0,ty);gdl(lado_x0,tx); gdl(lado_y0,ww); gdl(lado_y0,tx);gdl(lado_y0,ty);gdl(lado_xLx,ww);  gdl(lado_xLx,ty);gdl(lado_xLx,tx);gdl(lado_yLy,ww); gdl(lado_yLy,tx); gdl(lado_yLy,ty)];
  

if EExi==123
    cxi = [gdl(lado_x0,ww);
           gdl(lado_x0,ty);
           gdl(lado_x0,tx)]; 
elseif EExi==12
    cxi = [gdl(lado_x0,ww);
           gdl(lado_x0,ty)];
elseif EExi==0  
    cxi=NaN;
end

if EExf==123
    cxf = [gdl(lado_xLx,ww);
           gdl(lado_xLx,ty);
           gdl(lado_xLx,tx)]; 
elseif EExf==12   
    cxf = [gdl(lado_xLx,ww);
           gdl(lado_xLx,ty)];
elseif EExf==0 
    cxf=NaN;
end

if EEyi==123
    cyi = [gdl(lado_y0,ww);
           gdl(lado_y0,ty);
           gdl(lado_y0,tx)]; 
elseif EEyi==12 
    cyi = [gdl(lado_y0,ww);
           gdl(lado_y0,ty)];
elseif EEyi==0    
    cyi=NaN;
end
if EEyf==123
    cyf = [gdl(lado_yLy,ww);
           gdl(lado_yLy,ty);
           gdl(lado_yLy,tx)]; 
elseif EEyf==12 
    cyf = [gdl(lado_yLy,ww);
           gdl(lado_yLy,ty)];
    
elseif EEyf==0  
    cyf=NaN;   
end
c = [cxi;
     cxf;
     cyi;
     cyf];
TF = isnan(c);
f=find(TF==0);
c=c(f,1);
d = setdiff(1:ngdl,c)';





