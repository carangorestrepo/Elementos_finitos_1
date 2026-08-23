clear all
close all
clc

%% constantes que ayudaran en la lectura del codigo
X = 1; Y = 2; Z = 3;
ww = 1; tx = 2; ty = 3;

%% ========================================================================
% DATOS DE ENTRADA
% ========================================================================
E = 4700*sqrt(28)*1000;
nu = 0.25;
h = 4*0.1;
q = -(4.6*1.2 + 1.8*1.6 + 0.2*24);   % [kN/m^2]
escala = 100;

Lx = 4;
Ly = 4;
deltax = 0.1;
deltay = 0.1;

%% 123 empotrada, 12 apoyado, 0 libre
EExi = 123;
EExf = 123;
EEyi = 123;
EEyf = 123;

%% Fundacion de Winkler
kWinkler = 500;     % [kN/m^3]

%% Estado membranal inicial para P-Delta / pandeo de placa
% Convencion: Nx0, Ny0 positivos = compresion
Nx0  = 0;           % [kN/m]
Ny0  = 0;           % [kN/m]
Nxy0 = 0;           % [kN/m]

%% Integracion / malla
n = 2;
rho = 2.4028;

%% Zona cargada
xqi = 0;
xqf = Lx;
yqi = 0;
yqf = Ly;

%% ========================================================================
% GENERACION DE LA MALLA
% ========================================================================
Nx = round(Lx/deltax,0);
dv = round(Nx/n,0);
Nx = dv*n + 1;

Ny = round(Ly/deltay,0);
dv = round(Ny/n,0);
Ny = dv*n + 1;

x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
[Xe,Ye] = meshgrid(x,y);

figure
plot(Xe,Ye,'*r'), hold on
Bx = reshape(Xe,[],1);
By = reshape(Ye,[],1);
xnod = [Bx,By];

%% nudos elemento
nnoi = (1:Ny*Nx)';
nno = size(xnod,1);
Nudosg = reshape(nnoi,[Ny,Nx]);

recorridox = 2:Nx;
Nelex = size(recorridox,2);
recorridoy = 2:Ny;
Neley = size(recorridoy,2);
nef = Nelex*Neley;

LaG = zeros(nef,4);
e = 1;
LaGc = cell(nef,1);
xeg = cell(nef,1);
yeg = cell(nef,1);
xe = zeros(nef,4);
ye = zeros(nef,4);
cg = zeros(nef,2);

for ey = 1:Neley
    for ex = 1:Nelex
        LaGc{e} = Nudosg((recorridoy(ey)-1):recorridoy(ey), ...
                         (recorridox(ex)-1):recorridox(ex));
        xeg{e} = Xe((recorridoy(ey)-1):recorridoy(ey), ...
                    (recorridox(ex)-1):recorridox(ex));
        yeg{e} = Ye((recorridoy(ey)-1):recorridoy(ey), ...
                    (recorridox(ex)-1):recorridox(ex));

        LaG(e,:) = [LaGc{e}(1,1), LaGc{e}(1,2), ...
                    LaGc{e}(2,2), LaGc{e}(2,1)];
        xe(e,:) = [xeg{e}(1,1), xeg{e}(1,2), xeg{e}(2,2), xeg{e}(2,1)];
        ye(e,:) = [yeg{e}(1,1), yeg{e}(1,2), yeg{e}(2,2), yeg{e}(2,1)];

        cg(e,:) = [mean(xe(e,:)), mean(ye(e,:))];
        text(cg(e,X),cg(e,Y),num2str(e),'Color','b');
        plot(xe(e,[1:4,1]),ye(e,[1:4,1]))
        e = e+1;
    end
end
axis equal tight
title('Malla de elementos finitos')

%% ========================================================================
% CONDICIONES DE BORDE
% ========================================================================
lado_x0  = find(xnod(:,X) == 0);
lado_y0  = find(xnod(:,Y) == 0);
lado_xLx = find(xnod(:,X) == Lx);
lado_yLy = find(xnod(:,Y) == Ly);

nno  = size(xnod,1);
ngdl = 3*nno;
gdl  = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)'];

if EExi == 123
    cxi = [gdl(lado_x0,ww); gdl(lado_x0,ty); gdl(lado_x0,tx)];
elseif EExi == 12
    cxi = [gdl(lado_x0,ww); gdl(lado_x0,ty)];
else
    cxi = NaN;
end

if EExf == 123
    cxf = [gdl(lado_xLx,ww); gdl(lado_xLx,ty); gdl(lado_xLx,tx)];
elseif EExf == 12
    cxf = [gdl(lado_xLx,ww); gdl(lado_xLx,ty)];
else
    cxf = NaN;
end

if EEyi == 123
    cyi = [gdl(lado_y0,ww); gdl(lado_y0,ty); gdl(lado_y0,tx)];
elseif EEyi == 12
    cyi = [gdl(lado_y0,ww); gdl(lado_y0,ty)];
else
    cyi = NaN;
end

if EEyf == 123
    cyf = [gdl(lado_yLy,ww); gdl(lado_yLy,ty); gdl(lado_yLy,tx)];
elseif EEyf == 12
    cyf = [gdl(lado_yLy,ww); gdl(lado_yLy,ty)];
else
    cyf = NaN;
end

c = [cxi; cxf; cyi; cyf];
c = c(~isnan(c));
c = unique(c);
d = setdiff(1:ngdl,c)';

%% ========================================================================
% RESOLVER LOSA MINDLIN CON FUNCIONES EXACTAS DE TIMOSHENKO EN BORDES
% ========================================================================
resultado = EF_Mindlin_TE(xnod,LaG,gdl,c,d,E,nu,h,q,rho, ...
    kWinkler,Nx0,Ny0,Nxy0,xqi,xqf,yqi,yqf,escala);
