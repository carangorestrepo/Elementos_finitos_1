clear all
close
clc
%% constantes que ayudarán en la lectura del código
X = 1; Y = 2;


%%
Ee=24870062.3;
nue=0.2;
te=0.2;
Lx=2;
Ly=2;
deltax=0.05;
deltay=0.025;
n=2;

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
Bx = reshape(Xe,[],1);
By = reshape(Ye,[],1);
xnod=[Bx,By];
%% nudos elemento
nnoi =(1:Ny*Nx)';% numero de nodos
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
%text(Bx,By,num2str(nnoi)),hold on
%% nugos grilla elemento
Nudosg = reshape(nnoi,[Ny,Nx]);

recorridox=3:2:Nx;
Nelex=size(recorridox,2);
recorridoy=3:2:Ny;
Neley=size(recorridoy,2);
nef=Nelex*Neley;% numero de EFs (numero de filas de LaG)
LaG=zeros(nef,9);
e=1;
LaGc = cell(nef,1); 
xeg = cell(nef,1); 
yeg = cell(nef,1); 
xe =zeros(nef,9);
ye =zeros(nef,9);
cg = zeros(nef,2); % almacena el centro de gravedad de los EFs
for ey=1:(Neley)
    for ex=1:(Nelex)
       %%nudos por elemento
       LaGc{e}=Nudosg((recorridoy(ey)-2):recorridoy(ey),(recorridox(ex)-2):recorridox(ex)); 
       %coordenadas elemento
       xeg{e}=Xe((recorridoy(ey)-2):recorridoy(ey),(recorridox(ex)-2):recorridox(ex)); 
       yeg{e}=Ye((recorridoy(ey)-2):recorridoy(ey),(recorridox(ex)-2):recorridox(ex)); 
       % se determinan las coordenadas de los nodos el EF e
       LaG(e,:)=[LaGc{e}(1,1),LaGc{e}(1,2),LaGc{e}(1,3),LaGc{e}(2,3),LaGc{e}(3,3),LaGc{e}(3,2),LaGc{e}(3,1),LaGc{e}(2,1),LaGc{e}(2,2)];
       xe(e,:)=[xeg{e}(1,1),xeg{e}(1,2),xeg{e}(1,3),xeg{e}(2,3),xeg{e}(3,3),xeg{e}(3,2),xeg{e}(3,1),xeg{e}(2,1),xeg{e}(2,2)];
       ye(e,:)=[yeg{e}(1,1),yeg{e}(1,2),yeg{e}(1,3),yeg{e}(2,3),yeg{e}(3,3),yeg{e}(3,2),yeg{e}(3,1),yeg{e}(2,1),yeg{e}(2,2)];
       % se calcula la posición del centro de gravedad del EF e
       cg(e,:) = [mean(xe(e,:)),mean(ye(e,:))];
       % se escribe el numero del EF e
       text(cg(e,X), cg(e,Y), num2str(e), 'Color', 'b');
       plot(xe(e,[1:8,1]),ye(e,[1:8,1]))
       e=e+1;
    end
end
axis equal tight
title('Malla de elementos finitos');

%% defino apoyos 
ap=Nudosg(1,:);
%% defino cargas
carga=200*23/43;
f=zeros(nno*2,1);

f(Nudosg(:,1)*2-1,1)=carga;


d = setdiff(1:nno, ap);



