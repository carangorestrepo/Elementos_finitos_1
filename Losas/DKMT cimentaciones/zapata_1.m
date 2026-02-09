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
%% 2. CÁLCULO DE CARGAS DISTRIBUIDAS PERMANENTES
%% ============================================================

qsuelo   = (pd - hz) * gamas;  % Peso del suelo sobre la zapata
qzapata  = hz * gamac;         % Peso propio zapata
qcolumna = -pd * gamac;        % Peso columna transmitido

% Carga uniforme fuera de la zona de columna
qz = -qsuelo - qzapata;

%% Definición geométrica del contorno de la columna
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
area = 1/2 * abs( ...
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
%% 8. DEFINICIÓN DE CARGAS DISTRIBUIDAS
%% ============================================================

nef = size(LaG,1);  % Número total elementos
q   = zeros(nef,1); % Vector presión nodal

columnas  = find(in_col==1);
zapatass  = find(in_col==0);
sizecolumnas=size(columnas,1);
% Carga distribuida zapata
q(zapatass) = qz;

%% Dibujar contornos geométricos
plot([0 Lx Lx 0 0], [0 0 Ly Ly 0], 'k-', 'LineWidth', 2);
plot([columna(:,1); columna(1,1)], ...
     [columna(:,2); columna(1,2)], 'b--', 'LineWidth', 2);

title(['Malla con densidad zapata = ', num2str(densidad_zapata), ...
       ' m, columna = ', num2str(densidad_columna), ' m']);

xlabel('x [m]'); ylabel('y [m]');
legend('Zapata','Columna','Borde zapata','Borde columna',...
       'Location','northeastoutside');

%% ============================================================
%% 9. ESTADÍSTICAS DE MALLA
%% ============================================================

fprintf('Densidad zapata: %.2f m\n', densidad_zapata);
fprintf('Densidad columna: %.2f m\n', densidad_columna);
fprintf('Número nodos: %d\n', size(xnod, 1));
fprintf('Número elementos: %d\n', size(LaG, 1));

%% ============================================================
%% 10. DISTRIBUCIÓN DE CARGA EN LA BASE DE LA COLUMNA
%% ============================================================

areac = area(in_col);          % Áreas elementos columna
centroidesc = (xnod(cencol(:,1),:) + ...
               xnod(cencol(:,2),:) + ...
               xnod(cencol(:,3),:)) / 3;

xi = centroidesc(:,1);
yi = centroidesc(:,2);

% Área total columna
Acol = sum(areac);

% Centroide de la zona columna
xc = sum(areac .* xi) / Acol;
yc = sum(areac .* yi) / Acol;

% Coordenadas relativas al centroide
xr = xi - xc;
yr = yi - yc;

%% Integrales geométricas discretizadas
S0  = sum(areac);
Sxx = sum(areac .* xr.^2);
Syy = sum(areac .* yr.^2);
Sxy = sum(areac .* xr .* yr);

%% Sistema de equilibrio estático
Mmat = [ S0   0    0 ;
         0   Sxx  Sxy;
         0   Sxy  Syy ];

rhs = [P; My; Mx];

coef = Mmat \ rhs;

acc = coef(1);
bcc = coef(2);
ccc = coef(3);

% Distribución presión lineal columna
qcol = acc + bcc*xr + ccc*yr;

% Se adiciona peso propio columna
q(columnas) = qcol + qcolumna;

%% ============================================================
%% 11. VERIFICACIÓN DE EQUILIBRIO
%% ============================================================

Pcheck  = sum(qcol .* areac);
Mxcheck = sum(qcol .* areac .* yr);
Mycheck = sum(qcol .* areac .* xr);

fprintf('Error P  = %.6e\n', (Pcheck-P)/P);
fprintf('Error Mx = %.6e\n', (Mxcheck-Mx)/Mx);
fprintf('Error My = %.6e\n', (Mycheck-My)/My);
