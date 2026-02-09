clc; clear; close all;

%% 1. Definición de geometría
Lx = 2.2;      % Largo zapata [m]
Ly = 3;        % Ancho zapata [m]
hz=0.4;       %espesor zapata
lcx = 0.35;   % Largo columna [m]
lcy = 0.65;   % Ancho columna [m]
cx0 = 1.85;   % Posición x columna [m]
cy0 = 1.2;    % Posición y columna [m]
pd=2.5;       % profundidad de desplante
gamas=18;     % denisidad suelo
gamac=24;     % densidad concreto
Balastro=6000;
P=-150;
Mx=20;
My=15;

c=[];
%% constantes que ayudarÃ¡n en la lectura del cÃ³digo
X = 1; Y = 2; Z = 3; % un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; % lectura del codigo
escala = 100; % factor de escalamiento de la deformada
%%
E=24870062.3;
nu=0.25;
h=0.4;

rho=2.4028;
qsuelo=(pd-hz)*gamas;%% peso propio suelos
qzapata=hz*gamac;%% peso propio zPt
qcolumna=-pd*gamac;
qz=-qsuelo-qzapata;

% Coordenadas de la columna (rectángulo)
columna = [cx0, cy0; cx0+lcx, cy0; cx0+lcx, cy0+lcy; cx0, cy0+lcy];

%% 2. Densidades de malla (editables)
densidad_zapata = 0.1;   % Densidad en la zapata [m]
densidad_columna = 0.05;  % Densidad en la columna [m] (más fina)

%% 3. Generación de puntos
% Puntos en el borde de la zapata (sin duplicados)
borde_zapata = [linspace(0, Lx, round(Lx/densidad_zapata)+1)', zeros(round(Lx/densidad_zapata)+1, 1);
                Lx*ones(round(Ly/densidad_zapata), 1), linspace(densidad_zapata, Ly, round(Ly/densidad_zapata))';
                linspace(Lx, 0, round(Lx/densidad_zapata)+1)', Ly*ones(round(Lx/densidad_zapata)+1, 1);
                zeros(round(Ly/densidad_zapata)-1, 1), linspace(Ly-densidad_zapata, densidad_zapata, round(Ly/densidad_zapata)-1)'];

% Puntos en el borde de la columna (con densidad ajustable)
borde_columna = [linspace(cx0, cx0+lcx, round(lcx/densidad_columna)+1)', cy0*ones(round(lcx/densidad_columna)+1, 1);
                 (cx0+lcx)*ones(round(lcy/densidad_columna), 1), linspace(cy0+densidad_columna, cy0+lcy, round(lcy/densidad_columna))';
                 linspace(cx0+lcx, cx0, round(lcx/densidad_columna)+1)', (cy0+lcy)*ones(round(lcx/densidad_columna)+1, 1);
                 cx0*ones(round(lcy/densidad_columna)-1, 1), linspace(cy0+lcy-densidad_columna, cy0+densidad_columna, round(lcy/densidad_columna)-1)'];

% Puntos interiores en la zapata (excluyendo la columna)
[x_zapata, y_zapata] = meshgrid(densidad_zapata:densidad_zapata:Lx-densidad_zapata, densidad_zapata:densidad_zapata:Ly-densidad_zapata);
puntos_zapata = [x_zapata(:), y_zapata(:)];
en_columna = inpolygon(puntos_zapata(:,1), puntos_zapata(:,2), columna(:,1), columna(:,2));
puntos_zapata = puntos_zapata(~en_columna, :);

% Puntos interiores en la columna (con densidad_columna)
[x_col, y_col] = meshgrid(cx0+densidad_columna:densidad_columna:cx0+lcx-densidad_columna, cy0+densidad_columna:densidad_columna:cy0+lcy-densidad_columna);
puntos_columna = [x_col(:), y_col(:)];

% Unión de todos los puntos (eliminando duplicados)
puntos = uniquetol([borde_zapata; borde_columna; puntos_zapata; puntos_columna], 1e-6, 'ByRows', true);

%% 4. Triangulación de Delaunay
DT = delaunayTriangulation(puntos);
LaG = DT.ConnectivityList;%% nudos que conectan elementos
xnod = DT.Points;%% coordenadas de nodos

%% 5. Identificar elementos dentro de la columna
centroides = (xnod(LaG(:,1),:) + xnod(LaG(:,2),:) + xnod(LaG(:,3),:)) / 3;
xy=[xnod(LaG(:,1),:) , xnod(LaG(:,2),:) , xnod(LaG(:,3),:)];

area=1/2*abs(xnod(LaG(:,1),1).*(xnod(LaG(:,2),2)-xnod(LaG(:,3),2))+...
             xnod(LaG(:,2),1).*(xnod(LaG(:,3),2)-xnod(LaG(:,1),2))+...
             xnod(LaG(:,3),1).*(xnod(LaG(:,1),2)-xnod(LaG(:,2),2)));

in_col = inpolygon(centroides(:,1), centroides(:,2), columna(:,1), columna(:,2));

%% 6. Visualización
figure;
hold on; axis equal; grid on;

% Triángulos fuera de la columna
patch('Faces', LaG(~in_col,:), 'Vertices', xnod, 'FaceColor', [0.8 1 0.8], 'EdgeColor', 'k');

% Triángulos dentro de la columna
patch('Faces', LaG(in_col,:), 'Vertices', xnod, 'FaceColor', [1 0.8 0.8], 'EdgeColor', 'k');
cencol=LaG(in_col,:);
centroidesc = (xnod(cencol(:,1),:) + xnod(cencol(:,2),:) + xnod(cencol(:,3),:)) / 3;

nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)
q=zeros(nef,1);
columnas=find(in_col==1);

sizecolumna=size(columnas,1);
zapatass=find(in_col==0);
sizezapatass=size(zapatass,1);
nef  = size(LaG,1);  % numero de EFs (numero de filas de LaG)

q(zapatass,1)=qz*ones(sizezapatass,1);
% Bordes
plot([0 Lx Lx 0 0], [0 0 Ly Ly 0], 'k-', 'LineWidth', 2); % Zapata
plot([columna(:,1); columna(1,1)], [columna(:,2); columna(1,2)], 'b--', 'LineWidth', 2); % Columna

title(['Malla con densidad en zapata = ', num2str(densidad_zapata), ' m, columna = ', num2str(densidad_columna), ' m']);
xlabel('x [m]'); ylabel('y [m]');
legend('Zapata', 'columna', 'Borde zapata', 'Borde columna','Location'                                                                                                                                                                 ,'northeastoutside');

%% 7. Estadísticas
fprintf('Densidad zapata: %.2f m\n', densidad_zapata);
fprintf('Densidad columna: %.2f m\n', densidad_columna);
fprintf('Número de nodos: %d\n', size(xnod, 1));
fprintf('Número de elementos: %d\n', size(LaG, 1));
%figure
%plot3(centroidescc(:,1),centroidescc(:,2),ones(sizecolumna,1)*qcolumna+qx+qy+qp,'o')

areac = area(in_col);
cent = centroidesc;

xi = cent(:,1);
yi = cent(:,2);

Acol = sum(areac);

xc = sum(areac .* xi) / Acol;
yc = sum(areac .* yi) / Acol;

xr = xi - xc;
yr = yi - yc;

S0  = sum(areac);
Sxx = sum(areac .* xr.^2);
Syy = sum(areac .* yr.^2);
Sxy = sum(areac .* xr .* yr);

Mmat = [ S0   0    0 ;
          0  Sxx  Sxy;
          0  Sxy  Syy ];

rhs = [P; My; Mx];

coef = Mmat \ rhs;

acc = coef(1);
bcc = coef(2);
ccc = coef(3);

qcol = acc + bcc*xr + ccc*yr;

q(columnas) = qcol + qcolumna;
Pcheck  = sum(qcol .* areac);
Mxcheck = sum(qcol .* areac .* yr);
Mycheck = sum(qcol .* areac .* xr);

fprintf('Error P  = %.6e\n', (Pcheck-P)/P);
fprintf('Error Mx = %.6e\n', (Mxcheck-Mx)/Mx);
fprintf('Error My = %.6e\n', (Mycheck-My)/My);