clear all
close all
clc
%% constantes que ayudarán en la lectura del código
X = 1; Y = 2;


%%
E=24870062;
nu=0.25;
t=0.2;
Lx=2;
Ly=2;
deltax=0.05;
deltay=0.025;
n=2;
rho=0;

Nx=round(Lx/deltax,0);
dv=round(Nx/n,0);
Nx=dv*n+3;%% numero de nudos X

Ny=round(round(Ly/deltay,0)/n,0);
dv=round(Ny/n,0);
Ny=dv*n+3;%% numero de nudos Y

x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);
[Xe,Ye]=meshgrid(x,y);%% grilla coordenadas nudos elemento

plot(Xe,Ye, '*r'),hold on%% grafica nodos elemento finito
%Bx = reshape(Xe,[],1);
%By = reshape(Ye,[],1);
%xnod=[Bx,By];
%% nudos elemento
nnoi =(1:Ny*Nx)';% numero de nodos


%% nugos grilla elemento
%Nudosg = reshape(nnoi,[Ny,Nx]);
xegi=NaN(Ny,Nx);
yegi=NaN(Ny,Nx);
Nudosg=zeros(Ny,Nx);
a=0;
for i=1:Nx  
    if (-1)^i==-1
        if i==1
            Nudosg(:,i)=1:Ny;
            xegi(:,i)=Xe(1:Ny,i);
            yegi(:,i)=Ye(1:Ny,i);
        else
            Nudosg(:,i)=(a+1):(Ny+a);
            xegi(:,i)=Xe(1:Ny,i);
            yegi(:,i)=Ye(1:Ny,i);
        end
    else
        a=max(Nudosg(:,i-1));
        pos=1:2:Ny;
        Nudosg(pos,i)=(a+1):1:(a+1+(Ny+1)/2-1);
        xegi(pos,i)=Xe(pos,i);
        yegi(pos,i)=Ye(pos,i);
    end
    a=max(Nudosg(:,i));
end
xnodi = reshape(xegi,[],1);
ynodi = reshape(yegi,[],1);

vxnodii=isnan(xnodi);
vynodii=isnan(ynodi);

f=find(vxnodii==0);

xnod=[xnodi(f,1),ynodi(f,1)];%% coordenasas elementos finitos

By = reshape(Ye,[],1);

recorridox=3:2:Nx;
Nelex=size(recorridox,2);
recorridoy=3:2:Ny;
Neley=size(recorridoy,2);
nef=Nelex*Neley;% numero de EFs (numero de filas de LaG)
LaG=zeros(nef,8);
e=1;
LaGc = cell(nef,1); 
xeg = cell(nef,1); 
yeg = cell(nef,1); 
xe =zeros(nef,8);
ye =zeros(nef,8);
cg = zeros(nef,2); % almacena el centro de gravedad de los EFs
for ey=1:(Neley)
    for ex=1:(Nelex)
       %%nudos por elemento
       LaGc{e}=Nudosg((recorridoy(ey)-2):recorridoy(ey),(recorridox(ex)-2):recorridox(ex)); 
       %coordenadas elemento
       xeg{e}=Xe((recorridoy(ey)-2):recorridoy(ey),(recorridox(ex)-2):recorridox(ex)); 
       yeg{e}=Ye((recorridoy(ey)-2):recorridoy(ey),(recorridox(ex)-2):recorridox(ex)); 
       % se determinan las coordenadas de los nodos el EF e
       LaG(e,:)=[LaGc{e}(1,1),LaGc{e}(1,2),LaGc{e}(1,3),LaGc{e}(2,3),LaGc{e}(3,3),LaGc{e}(3,2),LaGc{e}(3,1),LaGc{e}(2,1)];
       xe(e,:)=[xeg{e}(1,1),xeg{e}(1,2),xeg{e}(1,3),xeg{e}(2,3),xeg{e}(3,3),xeg{e}(3,2),xeg{e}(3,1),xeg{e}(2,1)];
       ye(e,:)=[yeg{e}(1,1),yeg{e}(1,2),yeg{e}(1,3),yeg{e}(2,3),yeg{e}(3,3),yeg{e}(3,2),yeg{e}(3,1),yeg{e}(2,1)];
       % se calcula la posición del centro de gravedad del EF e
       cg(e,:) = [mean(xe(e,:)),mean(ye(e,:))];
       % se escribe el numero del EF e
       text(cg(e,X), cg(e,Y), num2str(e), 'Color', 'b','HorizontalAlignment','right');
       plot(xe(e,[1:8,1]),ye(e,[1:8,1]))
       e=e+1;
    end
end


axis equal tight
title('Malla de elementos finitos');
nno  =Nudosg(Ny,Nx); % numero de nodos (numero de filas de xnod)
sizexnod=size(xnod,1);
idxEF=1:nef;
idxNODO=(1:sizexnod)';% nudos malla elemetos finitos

%text(xnod(:,1),xnod(:,2),num2str(idxNODO)),hold on


mat=ones(nef,1);%% material por elemento finito
%% defino nodos de  apoyos 

idxNODOr=[Nudosg(1,:)';Nudosg(1,:)'];
dir_desp=[ones(Nx,1);2*ones(Nx,1)];%% defino direcion de desplazamiento
ac=[zeros(Nx,1);zeros(Nx,1)];% desplazamientos conocidos en los apoyos

%d = setdiff(1:nno, ap);
%% defino cargas puntuales
carga=4300/43;
cp=ones(Ny,1)*carga;
%cp(Nudosg(:,1)*2-1,1)=carga;
dir_cp=ones(Ny,1)*2;
idxNODOp=Nudosg(:,1);
%% cargas sitribuidas 
idxELEM=421;%% elemento
nodoijk=[43,42,41];
%% tix	tiy	tjx	tjy	tkx	tky
carga=[0,0,0,0,0,0];






