clc
clear
wx=20;
wy=20;
x=[1.5,2.2,1.7,1.5];%% separacion grilla x
y=[1.7,1.8,1.7,2.2];%% separacion grilla y

Ee=24870062;%4700*sqrt(28)*1000;
Ge=Ee/(2*(1+0.25));% ton/m^2  modulo de elasticidad

puntos_graficas=100;

%% escalas de Dibujo la estructura y su deformada
esc_def    = 1;          % escalamiento de la deformada
esc_faxial = 0.006;      % escalamiento del diagrama de axiales
esc_V      = 0.01;       % escalamiento del diagrama de cortantes
esc_M      = 0.01;      % escalamiento del diagrama de momentos

X=[0,cumsum(x,2)];
Y=[0,cumsum(y,2)];
[xa,ya]=meshgrid(X,Y);

sizex=size(X,2);% n�mero y
sizey=size(Y,2); % n�mero y
xnod=[reshape(xa,[sizex*sizey,1]),reshape(ya,[sizex*sizey,1])];
nno   =size(X,2)*size(Y,2); % n�mero de nodos
nodos = reshape(1:nno,[sizey,sizex]);

C=3*[nodos(:,1);nodos(:,sizex);nodos(1,:)';nodos(sizey,:)'];


nefx=(sizex-1)*sizey;    % numero de EFs en direccion X
nefy=(sizey-1)*sizex;    % numero de EFs en direccion Y

LaGx=zeros(nefx,2); %% nudos que conectan los elemetos Y
recorridox=1:sizey:nefx;

for ex=1:(sizex-1)
    LaGx(recorridox(ex):(recorridox(ex)+sizey-1),:)=nodos(:,ex:ex+1);
end
LaGy=zeros(nefy,2);  %% nudos que conectan los elemetos Y
recorridoy=1:sizex:nefy;

for ey=1:(sizey-1)
    LaGy(recorridoy(ey):(recorridoy(ey)+sizex-1),:)=nodos(ey:ey+1,:)';
end
LaG=[LaGx;LaGy]; %% nudos que conectan los elemetos 
bf=0.1;
d=0.4;
tw=0;
tf=0;



[Ae,Ixe,Iye,As2e,As3e,Je]=propiedades_geometrica_perfiles_v1(d,bf,tw,tf);
nef=nefx+nefy;
seccion=2;%Rectangular

E=ones(nef,1)*Ee;
G=ones(nef,1)*Ge;

EA=Ae(seccion).*E;
Ac2=As2e(seccion).*G;
Ac3=As3e(seccion).*G;
EIz=E.*Ixe(seccion).*ones(nef,1);
EIy=E.*Iye(seccion).*ones(nef,1);
GJ=G.*Je(seccion).*ones(nef,1);

Acy=Ac2;
Acz=Ac3;


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

elecargas=[1:nef]';
xyc=[xnod,zeros(nno,1)];
ninf=LaG;
carga=[qa,qa];
dcarga=ones(1,nef);
carga_nudo=0;
nudos_carga=1;
nudosapoyos=[nodos(:,1);nodos(:,sizex);nodos(1,:)';nodos(sizey,:)']';

%% ingreso tipo de apoyo 
tipoapoyo=ones(size(nudosapoyos,2),1)'*123;

[qe_loc,kmodal,mmodal,gama,GLKM,xye,elementos,nudos,GLe,we,Le,esc_def,esc_faxial,esc_V,esc_M3,esc_M2,matT,T,Ke]= porticos_3D(xyc,ninf,elecargas,carga,dcarga,carga_nudo,nudos_carga,nudosapoyos,tipoapoyo,EA,Ac2,Ac3,EIz,EIy,GJ);

