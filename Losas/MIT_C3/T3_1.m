clear all
close
clc
%% constantes que ayudar�n en la lectura del c�digo
X = 1; Y = 2;
ww= 1; tx= 2; ty= 3; % lectura del codigo
%%

E=210e9;
nu=0.3;
t=0.05;
Lx=2;
Ly=4;
deltax=0.05*2;
deltay=0.025*2;
n=2;
rho=2.3;

Nx=round(Lx/deltax,0);
dv=round(Nx/n,0);
Nx=dv*n+1;%% numero de nudos X

Ny=round(round(Ly/deltay,0)/n,0);
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
nef=Nelex*Neley*2;% numero de EFs (numero de filas de LaG)
LaG1=zeros(nef,3);
LaG2=zeros(nef,3);
LaG=zeros(nef,3);
recorridoLaG=1:2:nef*2;
e=1;
LaGc = cell(nef,1); 
xeg = cell(nef,1); 
yeg = cell(nef,1); 

xe1 =zeros(nef,3);
ye1 =zeros(nef,3);

xe2 =zeros(nef,3);
ye2 =zeros(nef,3);

xe =zeros(nef,3);
ye =zeros(nef,3);

cg1 = zeros(nef,2); % almacena el centro de gravedad de los EFs
cg2 = zeros(nef,2); % almacena el centro de gravedad de los EFs

cg = zeros(nef,2);
for ey=1:(Neley)
    for ex=1:(Nelex)
       %%nudos por elemento
       LaGc{e}=Nudosg((recorridoy(ey)-1):recorridoy(ey),(recorridox(ex)-1):recorridox(ex)); 
       %coordenadas elemento
       xeg{e}=Xe((recorridoy(ey)-1):recorridoy(ey),(recorridox(ex)-1):recorridox(ex)); 
       yeg{e}=Ye((recorridoy(ey)-1):recorridoy(ey),(recorridox(ex)-1):recorridox(ex)); 
       
       % se determinan las coordenadas de los nodos el EF 
       LaG1(e,:)=[LaGc{e}(1,2),LaGc{e}(2,1),LaGc{e}(1,1)];
       LaG2(e,:)=[LaGc{e}(1,2),LaGc{e}(2,2),LaGc{e}(2,1)]; 
       LaG(recorridoLaG(e),:)=LaG1(e,:);
       LaG(recorridoLaG(e)+1,:)=LaG2(e,:);
       
       xe1(e,:)=[xeg{e}(1,2),xeg{e}(2,1),xeg{e}(1,1)];
       xe2(e,:)=[xeg{e}(1,2),xeg{e}(2,2),xeg{e}(2,1)]; 
       
       xe(recorridoLaG(e),:)=xe1(e,:);
       xe(recorridoLaG(e)+1,:)=xe2(e,:);
       
       ye1(e,:)=[yeg{e}(1,2),yeg{e}(2,1),yeg{e}(1,1)];
       ye2(e,:)=[yeg{e}(1,2),yeg{e}(2,2),yeg{e}(2,1)];
       
       ye(recorridoLaG(e),:)=ye1(e,:);
       ye(recorridoLaG(e)+1,:)=ye2(e,:);
       
       ye(recorridoLaG(e),:)=ye1(e,:);
       ye(recorridoLaG(e)+1,:)=ye2(e,:);
       
       %LaG(e,:)=[LaGc{e}(1,1),LaGc{e}(1,2),LaGc{e}(1,3),LaGc{e}(2,3),LaGc{e}(3,3),LaGc{e}(3,2),LaGc{e}(3,1),LaGc{e}(2,1),LaGc{e}(2,2)];
       %xe(e,:)=[xeg{e}(1,1),xeg{e}(1,2),xeg{e}(1,3),xeg{e}(2,3),xeg{e}(3,3),xeg{e}(3,2),xeg{e}(3,1),xeg{e}(2,1),xeg{e}(2,2)];
       %ye(e,:)=[yeg{e}(1,1),yeg{e}(1,2),yeg{e}(1,3),yeg{e}(2,3),yeg{e}(3,3),yeg{e}(3,2),yeg{e}(3,1),yeg{e}(2,1),yeg{e}(2,2)];
       % se calcula la posici�n del centro de gravedad del EF e
       cg1(e,:) = [mean(xe1(e,:)),mean(ye1(e,:))];
       cg2(e,:) = [mean(xe2(e,:)),mean(ye2(e,:))];
       
       cg(recorridoLaG(e),:)=cg1(e,:);
       cg(recorridoLaG(e)+1,:)=cg2(e,:);
       
       % se escribe el numero del EF e
       text(cg1(e,X), cg1(e,Y), num2str(recorridoLaG(e)), 'Color', 'b');
       text(cg2(e,X), cg2(e,Y), num2str(recorridoLaG(e)+1), 'Color', 'b');
       
       plot(xe1(e,[1,2,3,1]),ye1(e,[1,2,3,1]),xe2(e,[1,2,3,1]),ye2(e,[1,2,3,1]))
       e=e+1;
    end
end
axis equal tight
title('Malla de elementos finitos');

mat=ones(nef,1);%% material por elemento finito

%% defino apoyos 
idxNODOr=[Nudosg(1,:)';Nudosg(1,:)'];
dir_desp=[ones(Nx,1);2*ones(Nx,1)];%% defino direcion de desplazamiento
ac=[zeros(Nx,1);zeros(Nx,1)];% desplazamientos conocidos en los apoyos

%% defino cargas
carga=4300/41;
cp=ones(Ny,1)*carga;

dir_cp=ones(Ny,1)*2;
idxNODOp=Nudosg(:,1);


%% cargas sitribuidas 
idxELEM=3121;%% elemento
nodoijk=[41,40];
%% tix	tiy	tjx	tjy	tkx	tky
carga=[0,0,0,0];


nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad
nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)

%ww= 1; tx= 2; ty= 3; % lectura del codigo
%cara izquierda 'RXI''EXI''LXI'
wxi=1;
txi=2;
ap='LXI';
%cara derecha
wxd=1;
txd=2;
%cara inferior
wyi=1;
tyi=3;
%cara superior
wys=1;
tys=3;


%% grados de libertad del desplazamiento conocidos y desconocidos
% determino los grados de libertad correspondientes a los bordes
lado_x0 = find(xnod(:,X) == 0);     lado_y0 = find(xnod(:,Y) == 0);
lado_x2 = find(xnod(:,X) == Lx);     lado_y4 = find(xnod(:,Y) == Ly);


%% empotrado
c = [ gdl(lado_x0,wxi); 
      gdl(lado_x0,ty); 
      gdl(lado_x0,txi);
      gdl(lado_x2,wxd);
      gdl(lado_x2,ty);
      gdl(lado_x2,txd);
      gdl(lado_y0,wyi); 
      gdl(lado_y0,tx);
      gdl(lado_y0,tyi);
      gdl(lado_y4,wys); 
      gdl(lado_y4,tx)
      gdl(lado_y4,tys)];
  
%apoyo derecho apoyado
if strcmp('AXI',ap)==2
c = [ gdl(lado_x0,wxi); 
      gdl(lado_x0,ty); %gdl(lado_x0,txi);
      gdl(lado_x2,wxd);
      gdl(lado_x2,ty);
      gdl(lado_x2,txd);
      gdl(lado_y0,wyi); 
      gdl(lado_y0,tx);
      gdl(lado_y0,tyi);
      gdl(lado_y4,wys); 
      gdl(lado_y4,tx)
      gdl(lado_y4,tys)];
end
%apoyo derecho libre
if strcmp('LXI',ap)==1
c = [ %gdl(lado_x0,wxi);%gdl(lado_x0,ty); %gdl(lado_x0,txi);
      gdl(lado_x2,wxd);
      gdl(lado_x2,ty);
      gdl(lado_x2,txd);
      gdl(lado_y0,wyi); 
      gdl(lado_y0,tx);
      gdl(lado_y0,tyi);
      gdl(lado_y4,wys); 
      gdl(lado_y4,tx)
      gdl(lado_y4,tys)];
end


  