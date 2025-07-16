clc
clear
%% constantes que ayudarán en la lectura del código
X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo
nx=5;
ny=6;
sx=1.5;
sy=1.7;
Lx=nx*sx;
Ly=ny*sy;
qx=20;
qy=30;
EExi=123;
EExf=123;
EEyi=123;
EEyf=123;

x=linspace(0,nx*sx,(nx+1));
y=linspace(0,ny*sy,(ny+1));
[Xx,Yy] = meshgrid(x,y);
nudos=1:(nx+1)*(ny+1);
nno=(nx+1)*(ny+1);
nudos_grilla=reshape(nudos,ny+1,nx+1);
plot(Xx,Yy,'b',Xx',Yy','b',Xx,Yy,'*r')

xx=reshape(Xx,(ny+1)*(nx+1),1);%% coordemadas nudos grilla x
yy=reshape(Yy,(ny+1)*(nx+1),1);%% coordemadas nudos grilla y
text(xx,yy,num2str(nudos'));

GLbx=zeros(nx*(ny+1),2);

recorridox=1:(nx+1);
e=1;
for iy=1:(ny+1)
    for ix=1:nx
        GLbx(e,:)=nudos_grilla(iy,recorridox(ix):recorridox(ix+1));
        e=e+1;
    end
end
GLby=zeros(ny*(nx+1),2);
recorridox=1:(ny+1);
e=1;
for ix=1:(nx+1) 
    for  iy=1:ny
        GLby(e,:)=nudos_grilla(recorridox(iy):recorridox(iy+1), ix);
        e=e+1;
    end
end


LaG=[GLbx;GLby];% grados de libertia de elemetoos
xnod=[xx,yy];%% coordenadas de nudos
qe=[qx*ones(nx*(ny+1),2);qy*ones(ny*(nx+1),2)];


nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 3*nno;        % numero de grados de libertad (tres por nodo)
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % nodos vs grados de libertad


%% grados de libertad del desplazamiento conocidos y desconocidos
% determino los grados de libertad correspondientes a los bordes
lado_x0 = nudos_grilla(:,1);     
lado_y0 = nudos_grilla(1,:)';
lado_xLx = nudos_grilla(:,nx+1);     
lado_yLy = nudos_grilla(ny+1,:)';


if EExi==123
    cxi = [lado_x0*3-2;
           lado_x0*3-1;
           lado_x0*3]; 
elseif EExi==3
    cxi = lado_x0*3;
elseif EExi==0  
    cxi=NaN;
end

if EExf==123
    cxf = [lado_xLx*3-2;
          lado_xLx*3-1;
          lado_xLx*3];
elseif EExf==3   
    cxf = lado_xLx*3;

elseif EExf==0 
    cxf=NaN;
end

if EEyi==123
    cyi = [lado_y0*3-2;
          lado_y0*3-1;
          lado_y0*3]; 
elseif EEyi==12 
    cyi = lado_y0*3;
elseif EEyi==0    
    cyi=NaN;
end
if EEyf==123
    cyf = [lado_yLy*3-2;
           lado_yLy*3-1;
           lado_yLy*3]; 
elseif EEyf==12 
    cyf = lado_yLy*3;
    
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
d = setdiff(1:nno*3,c)';