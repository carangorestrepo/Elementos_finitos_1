clc; clear; close all;

%% ============================================================
%% 1. DEFINICIÓN GEOMÉTRICA Y PROPIEDADES DEL PROBLEMA
%% ============================================================

Lx = 1.5;      % Largo zapata [m]
Ly = 1.7;        % Ancho zapata [m]
hz = 0.4;      % Espesor zapata [m]

% Dimensiones de la columna
lcx = 0.35;    % Largo columna [m]
lcy = 0.40;    % Ancho columna [m]

% Posición inferior izquierda de la columna dentro de la zapata
cx0 = Lx/2-lcx/2;%1.85;    
cy0 = Ly/2-lcy/2;%1.2;

% Parámetros geotécnicos y estructurales
pd = 2.5;      % Profundidad de desplante [m]
gamas = 18;    % Peso volumétrico suelo [kN/m³]
gamac = 24;    % Peso volumétrico concreto [kN/m³]
Balastro = 6000; % Módulo de balastro [kN/m³]

% Acciones aplicadas en la columna
P  = -120;     % Carga axial [kN]
Mx = 15;       % Momento respecto eje x
My = 20;       % Momento respecto eje y

c = [];

%% Constantes auxiliares para lectura del código
X = 1; Y = 2; Z = 3;       % Índices coordenadas
ww = 1; tx = 2; ty = 3;    % Índices grados de libertad placa Mindlin
escala = 100;              % Factor de escala deformada

%% Propiedades mecánicas del material
E  = 24870062.3;   % Módulo de elasticidad
nu = 0.25;         % Coeficiente de Poisson
h  = 0.4;          % Espesor estructural

rho = 2.4028;      % Densidad equivalente

%% ============================================================
%% 2. CALCULO DE CARGAS DISTRIBUIDAS PERMANENTES
%% ============================================================

% Convencion: cargas hacia abajo negativas
h_suelo = pd-hz;

if h_suelo < 0
    error('pd debe ser mayor o igual que hz.');
end

% Peso de zapata sobre toda la superficie
qzapata = -gamac*hz;

% Peso de suelo solamente fuera de la huella
qsuelo = -gamas*h_suelo;

% Altura real de pedestal sobre la cara superior de la zapata
hped = pd-hz;
qpedestal = -gamac*hped;

% Evitar doble contabilizacion si P ya incluye estos pesos
P_incluye_pedestal = false;
P_incluye_zapata = false;

if P_incluye_pedestal
    qpedestal = 0;
end

if P_incluye_zapata
    qzapata = 0;
end

% Presion uniforme fuera de la huella
qz = qzapata+qsuelo;

% Propiedades exactas de la huella
Ac = lcx*lcy;
xc = cx0+lcx/2;
yc = cy0+lcy/2;
Ix_col = lcx*lcy^3/12;
Iy_col = lcy*lcx^3/12;

% Convencion de momentos
signoMx = 1;
signoMy = 1;

qPM = @(x,y) P/Ac ...
    + signoMy*My/Iy_col.*(x-xc) ...
    + signoMx*Mx/Ix_col.*(y-yc);

%% Definicion geometrica del contorno de la columna
columna = [cx0, cy0;
           cx0+lcx, cy0;
           cx0+lcx, cy0+lcy;
           cx0, cy0+lcy];

%% ============================================================
%% 3. PARÁMETROS DE MALLADO
%% ============================================================

densidad_zapata  = 0.03;   % Tamaño característico elementos en zapata
densidad_columna = 0.03;  % Malla refinada en zona de columna

%% ============================================================
%% 4. GENERACIÓN DE NODOS DEL MODELO
%% ============================================================

% --- Borde exterior zapata ---
borde_zapata = [linspace(0, Lx, round(Lx/densidad_zapata)+1)', ...
                zeros(round(Lx/densidad_zapata)+1, 1);

                Lx*ones(round(Ly/densidad_zapata), 1), ...
                linspace(densidad_zapata, Ly, round(Ly/densidad_zapata))';

                linspace(Lx, 0, round(Lx/densidad_zapata)+1)', ...
                Ly*ones(round(Lx/densidad_zapata)+1, 1);

                zeros(round(Ly/densidad_zapata)-1, 1), ...
                linspace(Ly-densidad_zapata, densidad_zapata, ...
                round(Ly/densidad_zapata)-1)'];

% --- Borde columna ---
borde_columna = [linspace(cx0, cx0+lcx, round(lcx/densidad_columna)+1)', ...
                 cy0*ones(round(lcx/densidad_columna)+1, 1);

                 (cx0+lcx)*ones(round(lcy/densidad_columna), 1), ...
                 linspace(cy0+densidad_columna, cy0+lcy, ...
                 round(lcy/densidad_columna))';

                 linspace(cx0+lcx, cx0, round(lcx/densidad_columna)+1)', ...
                 (cy0+lcy)*ones(round(lcx/densidad_columna)+1, 1);

                 cx0*ones(round(lcy/densidad_columna)-1, 1), ...
                 linspace(cy0+lcy-densidad_columna, ...
                 cy0+densidad_columna, round(lcy/densidad_columna)-1)'];

% --- Nodos interiores zapata (excluyendo columna) ---
[x_zapata, y_zapata] = meshgrid(densidad_zapata:densidad_zapata:Lx-densidad_zapata,...
                                densidad_zapata:densidad_zapata:Ly-densidad_zapata);

puntos_zapata = [x_zapata(:), y_zapata(:)];

% Eliminación de nodos dentro de la columna
en_columna = inpolygon(puntos_zapata(:,1), puntos_zapata(:,2), ...
                       columna(:,1), columna(:,2));

puntos_zapata = puntos_zapata(~en_columna, :);

% --- Nodos interiores columna ---
[x_col, y_col] = meshgrid(cx0+densidad_columna:densidad_columna:cx0+lcx-densidad_columna,...
                          cy0+densidad_columna:densidad_columna:cy0+lcy-densidad_columna);

puntos_columna = [x_col(:), y_col(:)];

% --- Unión total nodos ---
puntos = uniquetol([borde_zapata; borde_columna; ...
                    puntos_zapata; puntos_columna], ...
                    1e-6, 'ByRows', true);

%% ============================================================
%% 5. GENERACIÓN DE ELEMENTOS (TRIANGULACIÓN DELAUNAY)
%% ============================================================

DT   = delaunayTriangulation(puntos);
LaG  = DT.ConnectivityList;   % Conectividad elementos
xnod = DT.Points;             % Coordenadas nodales

%% ============================================================
%% 6. IDENTIFICACIÓN DE ELEMENTOS DE LA COLUMNA
%% ============================================================

% Cálculo centroides
centroides = (xnod(LaG(:,1),:) + xnod(LaG(:,2),:) + xnod(LaG(:,3),:)) / 3;

% Área de cada elemento triangular
area_elem = 1/2 * abs( ...
       xnod(LaG(:,1),1).*(xnod(LaG(:,2),2)-xnod(LaG(:,3),2)) + ...
       xnod(LaG(:,2),1).*(xnod(LaG(:,3),2)-xnod(LaG(:,1),2)) + ...
       xnod(LaG(:,3),1).*(xnod(LaG(:,1),2)-xnod(LaG(:,2),2)));

% Identificación de elementos contenidos en la columna
in_col = inpolygon(centroides(:,1), centroides(:,2), ...
                   columna(:,1), columna(:,2));

%% ============================================================
%% 7. VISUALIZACIÓN DE MALLA
%% ============================================================

figure; hold on; axis equal; grid on;

% Elementos zapata
patch('Faces', LaG(~in_col,:), 'Vertices', xnod,...
      'FaceColor', [0.8 1 0.8], 'EdgeColor', 'k');

% Elementos columna
patch('Faces', LaG(in_col,:), 'Vertices', xnod,...
      'FaceColor', [1 0.8 0.8], 'EdgeColor', 'k');

cencol = LaG(in_col,:);

%% ============================================================
%% 8. PRESION REPRESENTATIVA POR ELEMENTO PARA VISUALIZACION
%% ============================================================

nef = size(LaG,1);
q = zeros(nef,1);

columnas = find(in_col);
zapatass = find(~in_col);
sizecolumnas = length(columnas);

for e = 1:nef
    xce = centroides(e,1);
    yce = centroides(e,2);

    if in_col(e)
        q(e) = qPM(xce,yce)+qzapata+qpedestal;
    else
        q(e) = qzapata+qsuelo;
    end
end

%% Dibujar contornos geometricos
plot([0 Lx Lx 0 0],[0 0 Ly Ly 0],'k-','LineWidth',2);
plot([columna(:,1);columna(1,1)], ...
     [columna(:,2);columna(1,2)],'b--','LineWidth',2);

title(['Malla con densidad zapata = ',num2str(densidad_zapata), ...
       ' m, columna = ',num2str(densidad_columna),' m']);

xlabel('x [m]'); ylabel('y [m]');
legend('Zapata','Columna','Borde zapata','Borde columna', ...
       'Location','northeastoutside');

%% ============================================================
%% 9. ESTADISTICAS DE MALLA
%% ============================================================

fprintf('Densidad zapata: %.2f m',densidad_zapata);
fprintf('Densidad columna: %.2f m',densidad_columna);
fprintf('Numero nodos: %d',size(xnod,1));
fprintf('Numero elementos: %d',size(LaG,1));
fprintf('Area mallada: %.8f m^2',sum(area_elem));
fprintf('Area teorica: %.8f m^2',Lx*Ly);

%% ============================================================
%% 10. VERIFICACION DE P, Mx Y My CON CUADRATURA
%% ============================================================

xw = TriGaussPoints(4);

Pcheck = 0;
Mxcheck = 0;
Mycheck = 0;

for e = columnas.'
    nod = LaG(e,:);
    xe = xnod(nod,1);
    ye = xnod(nod,2);

    detJe = (xe(2)-xe(1))*(ye(3)-ye(1)) ...
          - (xe(3)-xe(1))*(ye(2)-ye(1));

    for pp = 1:size(xw,1)
        xi = xw(pp,1);
        eta = xw(pp,2);
        wg = xw(pp,3);

        Ntri = [1-xi-eta,xi,eta];
        xgp = Ntri*xe;
        ygp = Ntri*ye;

        qgp = qPM(xgp,ygp);
        dA = detJe*wg;

        Pcheck = Pcheck+qgp*dA;
        Mxcheck = Mxcheck+qgp*(ygp-yc)*dA;
        Mycheck = Mycheck+qgp*(xgp-xc)*dA;
    end
end

fprintf('VERIFICACION DE ACCIONES DE COLUMNA');
fprintf('P objetivo/integrado = %.8f / %.8f kN',P,Pcheck);
fprintf('Mx objetivo/integrado = %.8f / %.8f kN*m',Mx,Mxcheck);
fprintf('My objetivo/integrado = %.8f / %.8f kN*m',My,Mycheck);

%% ============================================================
%% 11. EJECUCION DEL ELEMENTO DKMT CORREGIDO
%% ============================================================

run('EF_DKMT_corregido_misma_estructura.m');
