clear all
close
clc
%% constantes que ayudarán en la lectura del código
X = 1; Y = 2;

%%
E=24870062;
nu=0.25;
t=0.2;
Lx=2;
Ly=2;

rho=0;

Nx=max(5:4:(10*4+5));
Ny=max(5:4:(10*4+5));

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
text(Bx,By,num2str(nnoi)),hold on
%% nugos grilla elemento
Nudosg = reshape(nnoi,[Ny,Nx]);

recorridox=5:4:Nx;
Nelex=size(recorridox,2);
recorridoy=5:4:Ny;
Neley=size(recorridoy,2);
nef=Nelex*Neley*2;% numero de EFs (numero de filas de LaG)
%LaG=zeros(nef,9);
LaG1=zeros(nef,15);
LaG2=zeros(nef,15);
LaG=zeros(nef,15);
recorridoLaG=1:2:nef*2;
e=1;
LaGc = cell(nef,1); 
xeg = cell(nef,1); 
yeg = cell(nef,1); 


xe1 =zeros(nef,15);
ye1 =zeros(nef,15);

xe2 =zeros(nef,15);
ye2 =zeros(nef,15);


xe =zeros(nef*2,15);
ye =zeros(nef*2,15);

cg1 = zeros(nef,2); % almacena el centro de gravedad de los EFs
cg2 = zeros(nef,2); % almacena el centro de gravedad de los EFs

cg = zeros(nef*2,2);
for ey=1:(Neley)
    for ex=1:(Nelex)
       %%nudos por elemento
       LaGc{e}=Nudosg((recorridoy(ey)-4):recorridoy(ey),(recorridox(ex)-4):recorridox(ex)); 
       %coordenadas elemento
       xeg{e}=Xe((recorridoy(ey)-4):recorridoy(ey),(recorridox(ex)-4):recorridox(ex)); 
       yeg{e}=Ye((recorridoy(ey)-4):recorridoy(ey),(recorridox(ex)-4):recorridox(ex)); 
       % se determinan las coordenadas de los nodos el EF 
       LaG1(e,:)=[LaGc{e}(1,5),LaGc{e}(5,1),LaGc{e}(1,1),...
                  LaGc{e}(2,4),LaGc{e}(3,3),LaGc{e}(4,2)....
                  LaGc{e}(4,1),LaGc{e}(3,1),LaGc{e}(2,1)...
                  LaGc{e}(1,2),LaGc{e}(1,3),LaGc{e}(1,4)...
                  LaGc{e}(2,3),LaGc{e}(3,2),LaGc{e}(2,2)];
              
       LaG2(e,:)=[LaGc{e}(1,5),LaGc{e}(5,5),LaGc{e}(5,1),...
                  LaGc{e}(2,5),LaGc{e}(3,5),LaGc{e}(4,5)....
                  LaGc{e}(5,4),LaGc{e}(5,3),LaGc{e}(5,2)...
                  LaGc{e}(4,2),LaGc{e}(3,3),LaGc{e}(2,4)...
                  LaGc{e}(3,4),LaGc{e}(4,4),LaGc{e}(4,3)];

       LaG(recorridoLaG(e),:)=LaG1(e,:);
       LaG(recorridoLaG(e)+1,:)=LaG2(e,:);
       
       xe1(e,:)=[xeg{e}(1,5),xeg{e}(5,1),xeg{e}(1,1),...
                 xeg{e}(2,4),xeg{e}(3,3),xeg{e}(4,2)....
                 xeg{e}(4,1),xeg{e}(3,1),xeg{e}(2,1)...
                 xeg{e}(1,2),xeg{e}(1,3),xeg{e}(1,4)...
                 xeg{e}(2,3),xeg{e}(3,2),xeg{e}(2,2)];
       
       xe2(e,:)=[xeg{e}(1,5),xeg{e}(5,5),xeg{e}(5,1),...
                 xeg{e}(2,5),xeg{e}(3,5),xeg{e}(4,5)....
                 xeg{e}(5,4),xeg{e}(5,3),xeg{e}(5,2)...
                 xeg{e}(4,2),xeg{e}(3,3),xeg{e}(2,4)...
                 xeg{e}(3,4),xeg{e}(4,4),xeg{e}(4,3)];
       
       xe(recorridoLaG(e),:)=xe1(e,:);
       xe(recorridoLaG(e)+1,:)=xe2(e,:);
       
       ye1(e,:)=[yeg{e}(1,5),yeg{e}(5,1),yeg{e}(1,1),...
                 yeg{e}(2,4),yeg{e}(3,3),yeg{e}(4,2)....
                 yeg{e}(4,1),yeg{e}(3,1),yeg{e}(2,1)...
                 yeg{e}(1,2),yeg{e}(1,3),yeg{e}(1,4)...
                 yeg{e}(2,3),yeg{e}(3,2),yeg{e}(2,2)];
             
       ye2(e,:)=[yeg{e}(1,5),yeg{e}(5,5),yeg{e}(5,1),...
                 yeg{e}(2,5),yeg{e}(3,5),yeg{e}(4,5)....
                 yeg{e}(5,4),yeg{e}(5,3),yeg{e}(5,2)...
                 yeg{e}(4,2),yeg{e}(3,3),yeg{e}(2,4)...
                 yeg{e}(3,4),yeg{e}(4,4),yeg{e}(4,3)];
       
       ye(recorridoLaG(e),:)=ye1(e,:);
       ye(recorridoLaG(e)+1,:)=ye2(e,:);
       
       % se calcula la posición del centro de gravedad del EF e
       cg1(e,:) = [mean(xe1(e,:)),mean(ye1(e,:))];
       cg2(e,:) = [mean(xe2(e,:)),mean(ye2(e,:))];
       
       cg(recorridoLaG(e),:)=cg1(e,:);
       cg(recorridoLaG(e)+1,:)=cg2(e,:);
       
       % se escribe el numero del EF e
       text(cg1(e,X), cg1(e,Y), num2str(recorridoLaG(e)), 'Color', 'b');
       text(cg2(e,X), cg2(e,Y), num2str(recorridoLaG(e)+1), 'Color', 'b');
       
       plot(xe1(e,[1,4,5,6,2,7,8,9,3,10,11,12,1]),ye1(e,[1,4,5,6,2,7,8,9,3,10,11,12,1]),xe2(e,[1,4,5,6,2,7,8,9,3,10,11,12,1]),ye2(e,[1,4,5,6,2,7,8,9,3,10,11,12,1]))
       e=e+1;
    end
end
axis equal tight
title('Malla de elementos finitos');
mat=ones(nef,1);%% material por elemento finito

%% defino apoyos 
idxNODOr=[Nudosg(1,:)';Nudosg(1,:)'];%%nudos donde hay apoyos por direccion 
dir_desp=[ones(Nx,1);2*ones(Nx,1)];%% defino direcion de desplazamiento
ac=[zeros(Nx,1);zeros(Nx,1)];% desplazamientos conocidos en los apoyos

%% defino cargas
carga=4300/45;
cp=ones(Ny,1)*carga;%5 carga a aplicar por nodo

dir_cp=ones(Ny,1)*1;%% direccion carga
idxNODOp=Nudosg(45,:)';%% nudos donde se aplica carga


%% cargas sitribuidas 
idxELEM=221;%% elemento
nodoijk=[45,44,43,42,41];
%% tix	tiy	tjx	tjy	tkx	tky tlx	tly
carga=[0,0,0,0,0,0,0,0,0,0];



