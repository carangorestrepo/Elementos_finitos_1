%% ejerccio 11.3 escamilla 
clc
clear
%% Se describen las propiedades de los materiales
E=22*1000^2;   % [kPa] modulo de elasticidad del concreto (NSR 10  - art C.8.5) 
G=8.8*1000^2;  % [kPa] modulo de cortante (NSR 10  - art CR8.8.2)
ks=5/6;        % coeficiente de correccion del cortante para seccion rectangular
%% ingreso coordenadas de nudos de estructura
xyc=[0,0,3;5,-2,2;0,0,0;-3,3,3];
%% nudos que conectan elementos
ninf=[1,2;4,1;3,1];
%% Se describen las propiedades de la geometria secciones sectangulares
%% columnas bY
%% columnas hX
sec=[0.4,0.30;0.25,0.4;0.4,0.3];%%b,h
%% tipo de sección para cada elemento
tiposec=[1,2,3];
%% cargas aplicadas distribuida (gdl carga)
elecargas=[1,2];
% [kN/m] carga uniforme distribuida vertical
carga=[24,35];
%% direccion de la cargas
dcarga=[1,1];%% 2 en direccón horizontal , 1 vertical 
%% tipo de apoyo de apoyo
nudosapoyos=[3,2,4];%%1 horizontal, 2 vertical, 12 horozonal y vertical, 123 horizontal vertica y giro
%% tipo de restriccion 
tipoapoyo=[123,123,123];
%definicion de longitud nudo rigido por elemento
lalb=[0,0;0,0;0,0];
%lalb=[0.2,0;0,0.15;0,0.4];
[qe_loc,kmodal,mmodal,gama]= porticos_3D(E,G,ks,xyc,ninf,sec,tiposec,elecargas,carga,dcarga,nudosapoyos,tipoapoyo,lalb);