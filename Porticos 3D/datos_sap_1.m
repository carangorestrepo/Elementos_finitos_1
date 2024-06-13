clc
clear
close all
datos='sap2000_20_3';
%Connectivity - Frame hoja 1
[NUM1,TXT1,RAW1]=xlsread(datos,'Connectivity - Frame');%%
% Frame Loads - Distributed
[NUM2,TXT2,RAW2]=xlsread(datos,'Frame Loads - Distributed');%%
%Frame Section Assignments 2
[NUM3,TXT3,RAW3]=xlsread(datos,'Frame Section Assignments');%%
%Frame Props 01 - General
[NUM4,TXT4,RAW4]=xlsread(datos,'Frame Props 01 - General');%%
%Joint Coordinates hoja 4
[NUM5,TXT5,RAW5]=xlsread(datos,'Joint Coordinates');%%
%Joint Restraint Assignments hoja 6
[NUM6,TXT6,RAW6]=xlsread(datos,'Joint Restraint Assignments');%%
%MatProp 01 - General
[NUM7,TXT7,RAW7]=xlsread(datos,'MatProp 02 - Basic Mech Props');%%
%Joint Loads - Force
[NUM8,TXT8,RAW8]=xlsread(datos,'Joint Loads - Force');%%
mat=TXT7(4:end,1);
mat_sec=TXT4(4:end,2);
EG_materiales=NUM7(:,3:4);

sizeE=size(EG_materiales,1);
%% nudos que conectan elementos
ninf=str2double(TXT1(4:end,2:3));
nombres_elemetos=(TXT1(4:end,1));
cargas_elemetos=(TXT2(4:end,1));
%% coordenadas de nudos portico
xyc=NUM5(:,1:3);
secciones=TXT4(4:end,1);
sec=NUM4(:,1:2);
%% 
%Channel       1
%Rectangular   2
%Circle        3
%I/Wide Flange 4
%Pipe          5  Seccion circular hueca
%Box/Tube      6  Tubular hueca
%Angle         7  Perfil L 
%Tee           8  Perfil T
tiposeccion=TXT4(4:end,3);
size_tiposeccion=size(TXT4(4:end,3),1);
tiposeccion_1=zeros(size_tiposeccion,1);
tipo_mat=zeros(size_tiposeccion,2);
for i=1:size_tiposeccion
    if strcmp('Channel',tiposeccion(i))==1
        tiposeccion_1(i)=1;
    end
    if strcmp('Rectangular',tiposeccion(i))==1
        tiposeccion_1(i)=2;
    end
    if strcmp('Circle',tiposeccion(i))==1
        tiposeccion_1(i)=3;
    end
    if strcmp('I/Wide Flange',tiposeccion(i))==1
        tiposeccion_1(i)=4;
    end
    if strcmp('Pipe',tiposeccion(i))==1
        tiposeccion_1(i)=5;
    end
    if strcmp('Box/Tube',tiposeccion(i))==1
        tiposeccion_1(i)=6;
    end
    if strcmp('Angle',tiposeccion(i))==1
        tiposeccion_1(i)=7;
    end
    if strcmp('Tee',tiposeccion(i))==1
        tiposeccion_1(i)=8;
    end  
    tf3 = strcmp(mat,mat_sec(i));
    f=find(tf3==1);
    tipo_mat(i,:)=EG_materiales(f,:);
end
%%
secciones_nom=TXT3(4:end,4);
elementos=size(TXT3(4:end,4),1);
tiposec=zeros(elementos,1);
sec1=zeros(elementos,4);
elecargas=(1:elementos)';
nom_seccion=zeros(elementos,1);
carga=zeros(elementos,2);
Ae=zeros(elementos,1);
Ixe=zeros(elementos,1);
Iye=zeros(elementos,1);
Ac2e=zeros(elementos,1);
Ac3e=zeros(elementos,1);
Je=zeros(elementos,1);
E=zeros(elementos,1);
G=zeros(elementos,1);
for i=1:elementos
    tf = strcmp(secciones_nom(i),secciones);
    f=find(tf==1);
    tiposec(i,1)=f;
    nom_seccion(i,1)=tiposeccion_1(f,1); %% tipo de seccion 'Channel', 'Rectangular'
    sec1(i,:)=NUM4(f,1:4);%%dimesiones de secciones
    tf1 = strcmp(nombres_elemetos(i),cargas_elemetos);
    f1=find(tf1==1);
    if isempty(f1)==1
        carga(i,1)=0;
    else
        carga(i,1:2)=NUM2(f1,5:6);
    end
    d=sec1(i,1);
    bf=sec1(i,2);
    tf=sec1(i,3);
    tw=sec1(i,4);
    [A,Ix,Iy,Ac2,Ac3,J]=propiedades_geometrica_perfiles_v1(d,bf,tw,tf);
    Ae(i)=A(nom_seccion(i,1));
    Ixe(i)=Ix(nom_seccion(i,1));
    Iye(i)=Iy(nom_seccion(i,1));
    Ac2e(i)=Ac2(nom_seccion(i,1));
    Ac3e(i)=Ac3(nom_seccion(i,1));
    Je(i)=J(nom_seccion(i,1));
    E(i,1)=tipo_mat(f,1);
    G(i,1)=tipo_mat(f,2);
end
%% ingreso coordenadas de nudos de estructura
nudosapoyos=str2double(TXT6(4:end,1))';
%% ingreso tipo de apoyo 
tipoapoyo=ones(size(nudosapoyos,2),1)'*123;
%% nudos con cargas puntuales
nudos_carga=str2double(TXT8(4:end,1));
carga_nudo=NUM8(:,3);
%% Se describen las propiedades de los materiales
%E=21538106;   % [kPa] modulo de elasticidad del concreto (NSR 10  - art C.8.5) 
%G=8283887;  % [kPa] modulo de cortante (NSR 10  - art CR8.8.2)
%ks=5/6;        % coeficiente de correccion del cortante para seccion rectangular
%sec=[0.3,0.3;0.3,0.5];%%b,h
lalb=zeros(elementos,2);
dcarga=ones(1,elementos);
EA=Ae.*E;
Ac2=Ac2e.*G;
Ac3=Ac3e.*G;
EIz=E.*Ixe;
EIy=E.*Iye;
GJ=G.*Je;

Acy=Ac2;
Acz=Ac3;



[qe_loc,kmodal,mmodal,gama,GLKM,xye,elementos,nudos,GLe,we,Le,esc_def,esc_faxial,esc_V,esc_M3,esc_M2,matT,T,Ke]= porticos_3D(xyc,ninf,sec1,elecargas,elecargas,carga,dcarga,carga_nudo,nudos_carga,nudosapoyos,tipoapoyo,EA,Ac2,Ac3,EIz,EIy,GJ);

