clear all
close
clc

%% === CONSTANTES PARA LECTURA DEL CÓDIGO ===
X = 1; Y = 2; Z = 3;         % Coordenadas espaciales (índices para x, y, z)
ww = 1; tx = 2; ty = 3;     % Grados de libertad por nodo: desplazamiento vertical y rotaciones

%% === PARÁMETROS DE LA LOSA ===
E = 4700 * sqrt(28) * 1000; % Módulo de elasticidad del concreto (Pa)
nu = 0.25;                  % Coeficiente de Poisson
h = 0.1;                    % Espesor de la losa (m)
rho = 2700;             % Densidad (kg/m³)
k = 5/6;                % Factor de corrección por cortante
G = E / (2*(1 + nu));   % Módulo de cortante (MPa)

% Carga distribuida (kN/m^2). Combinación de cargas por área:
q = -(4.6*1.2 + 1.8*1.6 + 0.2*24); 

escala = 100;              % Escala para visualización de deformaciones

%% === DIMENSIONES DE LA LOSA Y MALLA ===
Lx = 1.5;                   % Longitud en x (m)
Ly = 1;                     % Longitud en y (m)
deltax = 0.05;              % Tamaño de elemento en x (m)
deltay = 0.05;              % Tamaño de elemento en y (m)

%% === CONDICIONES DE BORDE ===
% Codificación: 123 = empotrado, 12 = simplemente apoyado, 0 = libre
EExi = 0;   % Borde izquierdo (x = 0)
EExf = 12;  % Borde derecho (x = Lx)
EEyi = 12;  % Borde inferior (y = 0)
EEyf = 12;  % Borde superior (y = Ly)

%% === PARÁMETROS AUXILIARES ===
n = 2;       % Subdivisión para compatibilidad posterior
xqi = 0; xqf = Lx;          % Rango de aplicación de carga (no usado aquí)
yqi = 0; yqf = Ly;

%% === GENERACIÓN DE MALLA ===
% Ajuste del número de nodos para conectividad múltiplo de n
Nx = round(Lx/deltax, 0);
dv = round(Nx/n, 0);
Nx = dv*n + 1;

Ny = round(Ly/deltay, 0);
dv = round(Ny/n, 0);
Ny = dv*n + 1;

x = linspace(0, Lx, Nx);     % Coordenadas nodales en x
y = linspace(0, Ly, Ny);     % Coordenadas nodales en y
[Xe, Ye] = meshgrid(x, y);   % Malla regular de nodos

plot(Xe, Ye, '*r'), hold on  % Visualización de nodos

% Vectorización de coordenadas nodales
Bx = reshape(Xe, [], 1);
By = reshape(Ye, [], 1);
xnod = [Bx, By];             % Matriz de coordenadas nodales

nnoi = (1:Ny*Nx)';           % Lista de nodos
Nudosg = reshape(nnoi, [Ny, Nx]); % Matriz nodal ordenada por filas/columnas

recorridox = 2:Nx;           % Recorrido para generar elementos en x
recorridoy = 2:Ny;           % Recorrido para generar elementos en y
Nelex = length(recorridox); % Número de elementos en x
Neley = length(recorridoy); % Número de elementos en y
nef = Nelex * Neley;        % Número total de elementos finitos

% Inicialización de matrices de conectividad
LaG = zeros(nef, 4);
LaGc = cell(nef, 1);
xeg = cell(nef, 1);
yeg = cell(nef, 1);
xe = zeros(nef, 4);
ye = zeros(nef, 4);
cg = zeros(nef, 2);          % Centroide de cada elemento

e = 1;
for ey = 1:Neley
    for ex = 1:Nelex
        % Nodos del elemento
        LaGc{e} = Nudosg((recorridoy(ey)-1):recorridoy(ey), (recorridox(ex)-1):recorridox(ex));
        
        % Coordenadas del elemento
        xeg{e} = Xe((recorridoy(ey)-1):recorridoy(ey), (recorridox(ex)-1):recorridox(ex));
        yeg{e} = Ye((recorridoy(ey)-1):recorridoy(ey), (recorridox(ex)-1):recorridox(ex));

        % Asignación nodal al elemento en orden antihorario
        LaG(e,:) = [LaGc{e}(1,1), LaGc{e}(1,2), LaGc{e}(2,2), LaGc{e}(2,1)];
        xe(e,:) = [xeg{e}(1,1), xeg{e}(1,2), xeg{e}(2,2), xeg{e}(2,1)];
        ye(e,:) = [yeg{e}(1,1), yeg{e}(1,2), yeg{e}(2,2), yeg{e}(2,1)];

        % Centro de gravedad del EF
        cg(e,:) = [mean(xe(e,:)), mean(ye(e,:))];

        % Visualización y numeración del elemento
        text(cg(e,X), cg(e,Y), num2str(e), 'Color', 'b');
        plot(xe(e,[1:4,1]), ye(e,[1:4,1]))

        e = e + 1;
    end
end

axis equal tight
title('Malla de elementos finitos');

%% === GRADOS DE LIBERTAD POR NODO ===
nno = size(xnod,1);          % Número total de nodos
ngdl = 3 * nno;              % Total de grados de libertad (3 por nodo)
gdl = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)'];

%% === DETECCIÓN DE BORDES ===
lado_x0 = find(xnod(:,X) == 0);
lado_y0 = find(xnod(:,Y) == 0);
lado_xLx = find(xnod(:,X) == Lx);
lado_yLy = find(xnod(:,Y) == Ly);

%% === ASIGNACIÓN DE CONDICIONES DE BORDE ===
if EExi == 123
    cxi = [gdl(lado_x0,ww); gdl(lado_x0,ty); gdl(lado_x0,tx)];
elseif EExi == 12
    cxi = [gdl(lado_x0,ww); gdl(lado_x0,ty)];
elseif EExi == 0
    cxi = [];
end

if EExf == 123
    cxf = [gdl(lado_xLx,ww); gdl(lado_xLx,ty); gdl(lado_xLx,tx)];
elseif EExf == 12
    cxf = [gdl(lado_xLx,ww); gdl(lado_xLx,ty)];
elseif EExf == 0
    cxf = [];
end

if EEyi == 123
    cyi = [gdl(lado_y0,ww); gdl(lado_y0,ty); gdl(lado_y0,tx)];
elseif EEyi == 12
    cyi = [gdl(lado_y0,ww); gdl(lado_y0,ty)];
elseif EEyi == 0
    cyi = [];
end

if EEyf == 123
    cyf = [gdl(lado_yLy,ww); gdl(lado_yLy,ty); gdl(lado_yLy,tx)];
elseif EEyf == 12
    cyf = [gdl(lado_yLy,ww); gdl(lado_yLy,ty)];
elseif EEyf == 0
    cyf = [];
end

% Unión de grados de libertad restringidos y libres
c = [cxi; cxf; cyi; cyf];
c = c(~isnan(c));              % Elimina posibles NaN si existen

d = setdiff(1:ngdl, c)';       % Grados de libertad libres




