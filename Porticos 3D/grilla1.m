clc
clear
%% Se describen las propiedades de los materiales
E=4700*sqrt(28)*1000;   % [kPa] modulo de elasticidad del concreto (NSR 10  - art C.8.5) 
G=0.4*E;  % [kPa] modulo de cortante (NSR 10  - art CR8.8.2)
ks=5/6;        % coeficiente de correccion del cortante para seccion rectangular


numero_ejes_z=4;
z=[3,3,3];
zc=[0,cumsum(z)];
numero_de_ejes_y=3;
y=[4,4];
yc=[0,cumsum(y)];
numero_de_ejes_x=4;
x=[5,5,5];
xc=[0,cumsum(x)];
%% Se describen las propiedades de la geometria secciones sectangulares
%% columnas bY
%% columnas hX
c=[0.4,0.4];%% seccion de columnas
v=[0.4,0.4];%% seccion de vigas
sec=[c;v];
%%
[X,Y]=meshgrid(xc,yc);
%[Z]=meshgrid(zc);
xx=reshape(X,numero_de_ejes_x*numero_de_ejes_y,1);%% coordenadas x
yy=reshape(Y,numero_de_ejes_x*numero_de_ejes_y,1);%% coordenadas y
zz=ones(numero_de_ejes_x*numero_de_ejes_y,1);%% coordenadas z

val=1:(numero_ejes_z):(numero_ejes_z*numero_de_ejes_y*numero_de_ejes_x-(numero_ejes_z-1));
xyn=reshape(val,[numero_de_ejes_y,numero_de_ejes_x]);%%nudos en la base del portico
%% separo memoria
elemnetos  = cell(numero_ejes_z,1);
columnas_nivel=(numero_de_ejes_x)*(numero_de_ejes_y);%% # columnas por nivel
columnas_niveles=columnas_nivel*(numero_ejes_z-1);%% # numero total de columnas

recorrido_colu=1:columnas_nivel:columnas_nivel*(numero_ejes_z-1);
recorrido_coor=1:columnas_nivel:columnas_nivel*(numero_ejes_z);
columnas=zeros(columnas_niveles,2);
xyc=zeros(columnas_nivel*(numero_ejes_z),3);
nudos_coodenadas=zeros(columnas_nivel*(numero_ejes_z),1);%%
%%
for i=1:numero_ejes_z    
    elemnetos{i}=xyn+i-1;
    if i<numero_ejes_z
        columnas(recorrido_colu(i):(recorrido_colu(i)+columnas_nivel-1),1:2)=[(val+i-1)',(val+i)'];%% nudos que conectan columnas
    end
    xyc(recorrido_coor(i):(recorrido_coor(i)+columnas_nivel-1),:)=[xx,yy,zz*zc(i)];%% coordenadas nudos
    nudos_coodenadas(recorrido_coor(i):(recorrido_coor(i)+columnas_nivel-1),:)=val+i-1;%% coordenadas nudos
end

vigas_nivel=(numero_de_ejes_x-1)*(numero_de_ejes_y)+(numero_de_ejes_x)*(numero_de_ejes_y-1);
nvigasx=zeros((numero_de_ejes_x-1)*(numero_de_ejes_y),2);%%numero de vigas x por nivel
nvigasy=zeros((numero_de_ejes_x)*(numero_de_ejes_y-1),2);%%numero de vigas y por nivel
recorridox=1:numero_de_ejes_x;
recorridoy=1:numero_de_ejes_y;
recorridoe=1:numero_de_ejes_y:(numero_de_ejes_y*(numero_de_ejes_x-1));
recorridoe2=1:numero_de_ejes_x:(numero_de_ejes_x*(numero_de_ejes_y-1));
recorridoniveles=1:vigas_nivel:(vigas_nivel*(numero_ejes_z-1));
vigas=zeros(vigas_nivel*(numero_ejes_z-1),2);
for i=1:(numero_ejes_z-1)
    for n=1:(numero_de_ejes_x-1)
        nvigasx(recorridoe(n):(recorridoe(n)+numero_de_ejes_y-1),1:2)=elemnetos{2}(:,recorridox(n):(recorridox(n)+1))+i-1;
    end
    for n=1:(numero_de_ejes_y-1)
        nvigasy(recorridoe2(n):(recorridoe2(n)+numero_de_ejes_x-1),1:2)=(elemnetos{2}(recorridoy(n):(recorridoy(n)+1),:))'+i-1;
    end
    nxy=[nvigasx;nvigasy];
    vigas(recorridoniveles(i):(recorridoniveles(i)+vigas_nivel-1),:)=nxy;%% nudos que conectan vigas
end
a=1; 

[B,I] = sort(nudos_coodenadas);

 xyc= xyc(I,:);

%% nudos que conectan elementos
ninf=[columnas;vigas];
%%nudos de apoyos
nudosapoyos=val'; 
%% tipo de restriccion 
tipoapoyo=ones(1,numero_de_ejes_x*numero_de_ejes_y)*123;
%% tipo de sección para cada elemento
tiposec=[ones(1,columnas_niveles,1)*1,ones(1,vigas_nivel*(numero_ejes_z-1))*2];
%% cargas aplicadas distribuida (gdl carga)
carga=[ones(1,columnas_niveles,1)*0,ones(1,vigas_nivel*(numero_ejes_z-1))*20];
elecargas=1:(columnas_niveles+vigas_nivel*(numero_ejes_z-1));
%% direccion de la cargas
dcarga=[ones(1,columnas_niveles,1)*1,ones(1,vigas_nivel*(numero_ejes_z-1))*1];
%% definicion de longitud nudo rigido por elemento
lalb=[[ones(1,columnas_niveles,1)*0,ones(1,vigas_nivel*(numero_ejes_z-1))*0]',[ones(1,columnas_niveles,1)*0,ones(1,vigas_nivel*(numero_ejes_z-1))*0]'];
[qe_loc,kmodal,mmodal,gama,GLKM,xye,elementos,nudos,GLe]= porticos_3D(E,G,ks,xyc,ninf,sec,tiposec,elecargas,carga,dcarga,nudosapoyos,tipoapoyo,lalb);

