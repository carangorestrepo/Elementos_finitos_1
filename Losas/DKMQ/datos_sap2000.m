clc
clear
close all
%% defino las variables/constantes
X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo

datos='semicircular-1';
[NUM1,TXT1,RAW1]=xlsread(datos,'Joint Coordinates');%%
[NUM2,TXT2,RAW2]=xlsread(datos,'Connectivity - Area');%%

xnod=NUM1(:,1:2);
LaG =str2double(TXT2(4:end,3:6));


%% grados de libertad del desplazamiento conocidos y desconocidos
% determino los grados de libertad correspondientes a los bordes
lado_x0 = find(xnod(:,X) == 0);     
lado_y0 = find(xnod(:,Y) == 0);

% Definimos la geometria de la losa
%losa
% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad
nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)


c = [ gdl(lado_x0,ww); 
      gdl(lado_x0,ty);
      gdl(lado_x0,tx);
      gdl(lado_y0,ww); 
      gdl(lado_y0,tx);
      gdl(lado_y0,ty)];

d = setdiff(1:ngdl,c)';

E  = 10.92;      % [kPa]    modulo de elasticidad = 210GPa
nu = 0.3;         %         coeficiente de Poisson
h  = 0.1;        % [m]     espesor de la losa
q  = -1;         % [kN/m^2] carga
escala = 0.0001; % factor de escalamiento de la deformada

xqi=0;
xqf=5;
yqi=0;
yqf=5;
