       clear all
close
clc
%% constantes que ayudarán en la lectura del código
X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo

%%
E=4700*sqrt(28)*1000;
rho=2.4028;
nu=0.3;
h=4*0.1;
q  = -(4.6*1.2+1.8*1.6+0.2*24);         % [kN/m^2] carga
escala = 100; % factor de escalamiento de la deformada
Lx=4;
Ly=4;

EExi=123;
EExf=123;
EEyi=123;
EEyf=123;
datos='Libro1';
%Connectivity - Area
[NUM1,TXT1,RAW1]=xlsread(datos,'Connectivity - Area');%%
%Joint Coordinates
[NUM2,TXT2,RAW2]=xlsread(datos,'Joint Coordinates');%%


%% nudos que conectan elementos
LaG=str2double(TXT1(4:end,3:6));
%% numero de nudos por elemento
Numero_de_nodos_elem= (NUM1(1:end,1));
%% coordenadas de nodos
xnod=NUM2(:,1:2);


% grados de libertad del desplazamiento conocidos y desconocidos
% determino los grados de libertad correspondientes a los bordes
lado_x0 = find(xnod(:,X) == 0);     
lado_y0 = find(xnod(:,Y) == 0);
lado_xLx = find(xnod(:,X) == Lx);     
lado_yLy = find(xnod(:,Y) == Ly);

% Definimos la geometria de la losa
%losa
% ya tenemos en la memoria las variables
% xnod - posicion (x,y) de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos
nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad

%c = [ gdl(lado_x0,ww); gdl(lado_x0,ty);gdl(lado_x0,tx); gdl(lado_y0,ww); gdl(lado_y0,tx);gdl(lado_y0,ty);gdl(lado_xLx,ww);  gdl(lado_xLx,ty);gdl(lado_xLx,tx);gdl(lado_yLy,ww); gdl(lado_yLy,tx); gdl(lado_yLy,ty)];
  

if EExi==123
    cxi = [gdl(lado_x0,ww);
           gdl(lado_x0,ty);
           gdl(lado_x0,tx)]; 
elseif EExi==12
    cxi = [gdl(lado_x0,ww);
           gdl(lado_x0,ty)];
elseif EExi==0  
    cxi=NaN;
end

if EExf==123
    cxf = [gdl(lado_xLx,ww);
           gdl(lado_xLx,ty);
           gdl(lado_xLx,tx)]; 
elseif EExf==12   
    cxf = [gdl(lado_xLx,ww);
           gdl(lado_xLx,ty)];
elseif EExf==0 
    cxf=NaN;
end

if EEyi==123
    cyi = [gdl(lado_y0,ww);
           gdl(lado_y0,ty);
           gdl(lado_y0,tx)]; 
elseif EEyi==12 
    cyi = [gdl(lado_y0,ww);
           gdl(lado_y0,ty)];
elseif EEyi==0    
    cyi=NaN;
end
if EEyf==123
    cyf = [gdl(lado_yLy,ww);
           gdl(lado_yLy,ty);
           gdl(lado_yLy,tx)]; 
elseif EEyf==12 
    cyf = [gdl(lado_yLy,ww);
           gdl(lado_yLy,ty)];
    
elseif EEyf==0  
    cyf=NaN;   
end
c = [cxi;
     cxf;
     cyi;
     cyf];
TF = isnan(c);
f=find(TF==0);
c=c(f,1);
d = setdiff(1:ngdl,c)';






