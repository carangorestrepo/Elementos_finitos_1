clear
close
% Coordenadas del polígono cerrado (ejemplo)
xy=[0.5,0;
    6,0;
    6,3;
    6,3;
    5.5,3;
    5.5,6;
    4.4,6;
    4.4,9;
    3.5,9;
    3.5,12;
    3.2,12;
    3.2,15;
    1.72,15;
    1.72,12;
    1.72,9;
    1.233,9;
    1.233,6;
    1,6;1,3;
    0.5,3;
    0.5,0];

x = xy(:,1)';
y = xy(:,2)';
uniquex=unique(x);
uniquey=unique(y);
% Generar una malla rectangular que cubre el polígono
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);
% Definir el tamaño de los elementos de la malla
dx = 0.3; % Tamaño en dirección x
dy = 0.3; % Tamaño en dirección y
xx=vector(dx,uniquex);
yy=vector(dx,uniquey);
% Crear vectores de coordenadas para la malla
[X, Y] = meshgrid(xx, yy);
% Crear una máscara para determinar si cada punto de la malla está dentro del polígono
inPolygon = inpolygon(X, Y, x, y);
xn=sum(inPolygon,2);
yn=sum(inPolygon,1);

uniquexn=unique(xn);
uniqueyn=unique(yn);

nn=sum((uniquexn-1).*((uniqueyn-1)'));%% numero de elemetos finitos
% Filtrar la malla para mantener solo los puntos dentro del polígono
sizexx=size(X);
X1=nan(sizexx);
Y1=nan(sizexx);
[f,c]=find(inPolygon==1);
sizef=size(f,1);
sizeinPolygon=size(inPolygon);
nnoi=1:sizef;
%% nugos grilla elemento
Nudosg=nan(sizeinPolygon);
for i=1:sizef
    Nudosg(f(i),c(i))=nnoi(i);%% nudos elemetos
    X1(f(i),c(i))=X(f(i),c(i));%% coordenadas X elemento
    Y1(f(i),c(i))=Y(f(i),c(i));%% coordenadas X elemento
end
%X_filtered = X(inPolygon);
%Y_filtered = Y(inPolygon);
% Graficar el polígono y la malla
figure;
plot(x, y, '-o');
hold on;
plot(X1, Y1, 'b.');

xlabel('X');
ylabel('Y');
title('Malla generada para un polígono cerrado');
%legend('show');
grid on;
axis equal;


%% nudos elemento
%nnoi =(1:Ny*Nx)';% numero de nodos
nno  = sizef; % numero de nodos (numero de filas de xnod)
%text(Bx,By,num2str(nnoi)),hold on
%% nugos grilla elemento
%Nudosg = reshape(nnoi,[Ny,Nx]);

%recorridox=2:1:Nx;
Nelex=sizexx(1,2);
%recorridoy=2:1:Ny;
Neley=sizexx(1,1);
nef=nn;% numero de EFs (numero de filas de LaG)
LaG=zeros(nef,4);
e=1;
LaGc = cell(nef,1); 
xeg = cell(nef,1); 
yeg = cell(nef,1); 
xe =zeros(nef,4);
ye =zeros(nef,4);
cg = zeros(nef,2); % almacena el centro de gravedad de los EFs
for ey=1:(Neley)
    fx=isnan(X1(ey,:));
    recorridox=find(fx==0);
    recorridox=recorridox(2:end);
    Nelex=size(recorridox,2);
    for ex=1:(Nelex)

       fy=isnan(X1(:,ex));
       recorridoy=find(fy==0)';
       recorridoy=recorridoy(2:end);
       Neley=size(recorridoy,2);
       %%nudos por elemento
       LaGc{e}=Nudosg((recorridoy(ey)-1):recorridoy(ey),(recorridox(ex)-1):recorridox(ex)); 
       %coordenadas elemento
       xeg{e}=X1((recorridoy(ey)-1):recorridoy(ey),(recorridox(ex)-1):recorridox(ex)); 
       yeg{e}=Y1((recorridoy(ey)-1):recorridoy(ey),(recorridox(ex)-1):recorridox(ex)); 
       % se determinan las coordenadas de los nodos el EF e
       LaG(e,:)=[LaGc{e}(1,1),LaGc{e}(1,2),LaGc{e}(2,2),LaGc{e}(2,1)];
       xe(e,:)=[xeg{e}(1,1),xeg{e}(1,2),xeg{e}(2,2),xeg{e}(2,1)];
       ye(e,:)=[yeg{e}(1,1),yeg{e}(1,2),yeg{e}(2,2),yeg{e}(2,1)];
       % se calcula la posiciÃ³n del centro de gravedad del EF e
       cg(e,:) = [mean(xe(e,:)),mean(ye(e,:))];
       % se escribe el numero del EF e
       text(cg(e,1), cg(e,2), num2str(e), 'Color', 'b'), hold on;
       plot(xe(e,[1:4,1]),ye(e,[1:4,1]))
       e=e+1;
       if e==161
           a=1
       end
    end
    an=1
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
idxELEM=1561;%% elemento
nodoijk=[41,40];
%% tix	tiy	tjx	tjy	tkx	tky
carga=[0,0,0,0];



