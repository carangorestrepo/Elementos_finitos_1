clc
clear
wx=20;
wy=20;
x=[1.5,1.5];%% separacion grilla x
y=[1.7,1.7];%% separacion grilla y

Ee=24870062;%4700*sqrt(28)*1000;
Ge=Ee/(2*(1+0.25));% ton/m^2  modulo de elasticidad

puntos_graficas=100;

%% escalas de Dibujo la estructura y su deformada
esc_def    = 10;          % escalamiento de la deformada
esc_faxial = 0.006;      % escalamiento del diagrama de axiales
esc_V      = 0.01;       % escalamiento del diagrama de cortantes
esc_M      = 0.01;      % escalamiento del diagrama de momentos

X=[0,cumsum(x,2)];
Y=[0,cumsum(y,2)];
[xa,ya]=meshgrid(X,Y);

sizex=size(X,2);% número y
sizey=size(Y,2); % número y
xnod=[reshape(xa,[sizex*sizey,1]),reshape(ya,[sizex*sizey,1])];
nno   =size(X,2)*size(Y,2); % número de nodos
nodos = reshape(1:nno,[sizey,sizex]);

C=3*[nodos(:,1);nodos(:,sizex);nodos(1,:)';nodos(sizey,:)'];


nefx=(sizex-1)*sizey;    % numero de EFs X
nefy=(sizey-1)*sizex;    % numero de EFs Y

LaGx=zeros(nefx,2);
recorridox=1:sizey:nefx;

for ex=1:(sizex-1)
    LaGx(recorridox(ex):(recorridox(ex)+sizey-1),:)=nodos(:,ex:ex+1);
end
LaGy=zeros(nefy,2);
recorridoy=1:sizex:nefy;

for ey=1:(sizey-1)
    LaGy(recorridoy(ey):(recorridoy(ey)+sizex-1),:)=nodos(ey:ey+1,:)';
end
LaG=[LaGx;LaGy];
bf=0.1;
d=0.4;
tw=0;
tf=0;

[Ae,Ixe,Iye,As2e,As3e,Je]=propiedades_geometrica_perfiles_v1(d,bf,tw,tf);
nef=nefx+nefy;
seccion=2;%Rectangular
As2=ones(nef,1)*As2e(seccion);
Iz=ones(nef,1)*Ixe(seccion);
Jx=ones(nef,1)*Je(seccion)*0;
E=ones(nef,1)*Ee;
G=ones(nef,1)*Ge;

% gdl: grados de libertad
ngdl = nno*3;          % numero de grados de libertad

xw1=[[0,x];
     [x,0]];
 
yw1=[[0,y];
     [y,0]];
sx=xw1(1,:)/2+xw1(2,:)/2;
sy=yw1(1,:)/2+yw1(2,:)/2;

[Sx,Su]=meshgrid(sx,sy);

Wx=ones(nefx,1)*wx;
Wy=ones(nefy,1)*wy;

qa=[Wx;Wy];
