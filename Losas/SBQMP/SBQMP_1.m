%% Limpieza inicial
%clear, clc, close all 

%% Definición de constantes para mejor legibilidad
X = 1; Y = 2; Z = 3;   % Coordenadas
ww = 1; tx = 2; ty = 3; % Grados de libertad por nodo

%% Definición de malla (ejemplo cuadrilátero distorsionado)
%xnod = [0, 0;   % Nodo 1
%        2, 0;   % Nodo 2
%        3, 2;   % Nodo 3
%        1, 3];  % Nodo 4
%LaG = [1 2 3 4]; % Definición del elemento

nno = size(xnod,1);  % Número de nodos
ngdl = 3*nno;        % Grados de libertad totales
gdl = [(1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)']; % Asignación de GDL
nef = size(LaG,1);   % Número de elementos

%% Visualización de la malla
figure;
hold on;
for e = 1:nef
   plot(xnod(LaG(e,[1 2 3 4 1]), X), xnod(LaG(e,[1 2 3 4 1]), Y), 'b-');
   cgx = mean(xnod(LaG(e,:), X));
   cgy = mean(xnod(LaG(e,:), Y));
   text(cgx, cgy, num2str(e), 'Color', 'r');
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
title('Malla de elementos finitos');

%% Matrices constitutivas (según ecuaciones del documento)
% Matriz de rigidez a flexión [D_b] (ecuación 19)
D_b = (E*h^3)/(12*(1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu)/2];

% Matriz de rigidez a cortante [D_s] (ecuación 19)
D_s = k * h * G * eye(2);

% Matriz de rigidez total [D]
D = blkdiag(D_b, D_s);

% Matriz de densidad de masa [T] (ecuación 20)
T = rho * diag([h, h^3/12, h^3/12]);

%% Funciones de forma y derivadas para cuadrilátero bilineal
Nforma = @(xi,eta) [(1-xi)*(1-eta)/4; (1+xi)*(1-eta)/4; 
                   (1+xi)*(1+eta)/4; (1-xi)*(1+eta)/4];

dN_dxi = @(xi,eta) [-(1-eta)/4; (1-eta)/4; (1+eta)/4; -(1+eta)/4];
dN_deta = @(xi,eta) [-(1-xi)/4; -(1+xi)/4; (1+xi)/4; (1-xi)/4];

%% Matrices [P] y [Q] según ecuaciones (9) y (10)
P = @(xi,eta) [...
    1, -xi, -eta, -xi^2/2, -(xi^2*eta)/2 + eta^3/12, -eta^2/2, ...
    -(xi*eta^2)/2 + xi^3/12, -xi*eta/2, xi/2, xi*eta/2, eta/2, xi*eta/2;
    0, 1, 0, xi, xi*eta, 0, xi^2/4 + eta^2/2, eta/2, 1/2, eta/2, 0, -eta/2;
    0, 0, 1, 0, eta^2/4 + xi^2/2, eta, xi*eta, xi/2, 0, -xi/2, 1/2, xi/2];

Q = @(xi,eta) [...
    0,0,0,1,eta,0,xi/2,0,0,0,0,0;
    0,0,0,0,eta/2,1,xi,0,0,0,0,0;
    0,0,0,0,2*xi,0,2*eta,1,0,0,0,0;
    0,0,0,0,0,0,0,0,1,eta,0,0;
    0,0,0,0,0,0,0,0,0,0,1,xi];

%% Integración numérica (Gauss 2x2)
n_gl = 2;
[x_gl, w_gl] = gauss_points(n_gl);

%% Ensamblaje de matrices globales
K = sparse(ngdl, ngdl);
M = sparse(ngdl, ngdl);

for e = 1:nef
    % Coordenadas de los nodos del elemento
    xe = xnod(LaG(e,:), X);
    ye = xnod(LaG(e,:), Y);
    
    % Construcción matriz [C] (ecuación 14)
    C = zeros(12,12);
    % Coordenadas naturales de los 4 nodos de un elemento cuadrilátero (Q4)
    natural_coords = [-1, -1;  % Nodo 1
                       1, -1;  % Nodo 2
                       1,  1;  % Nodo 3
                      -1,  1]; % Nodo 4
    
    %natural_coords=[xe,ye];
    
    for i = 1:4
        xi  = natural_coords(i, 1);
        eta = natural_coords(i, 2);

        % Evalúa la matriz P en (xi, eta)
        P_val = P(xi, eta);  % función que devuelve 3x12

        % Inserta en la matriz C
        C(3*(i-1)+1 : 3*i, :) = P_val;
    end
    
    invC = inv(C');
    
    % Inicialización matrices elementales
    K0 = zeros(12);
    M0 = zeros(12);
    
    for pp = 1:n_gl
        for qq = 1:n_gl
            xi_gl = x_gl(pp);
            eta_gl = x_gl(qq);
            
            % Jacobiano y su determinante
            dx_dxi = dN_dxi(xi_gl, eta_gl)' * xe;
            dy_dxi = dN_dxi(xi_gl, eta_gl)' * ye;
            dx_deta = dN_deta(xi_gl, eta_gl)' * xe;
            dy_deta = dN_deta(xi_gl, eta_gl)' * ye;
            
            J = [dx_dxi, dy_dxi; dx_deta, dy_deta];
            detJ = det(J);
            if detJ <= 0
                error('El determinante del Jacobiano es negativo o cero. Verifica tu malla.');
            end
            % Evaluación de [P] y [Q] en punto de Gauss
            P_val = P(xi_gl, eta_gl);
            Q_val = Q(xi_gl, eta_gl);
            
            % Contribución a matrices elementales (ecuaciones 15-18)
            K0 = K0 + Q_val' * D * Q_val * detJ * w_gl(pp) * w_gl(qq);
            M0 = M0 + P_val' * T * P_val * detJ * w_gl(pp) * w_gl(qq);
        end
    end
    
    % Ensamblaje a matrices globales
    idx = [gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:)];
    K(idx,idx) = K(idx,idx) + invC' * K0 * invC;
    M(idx,idx) = M(idx,idx) + invC' * M0 * invC;
end
fprintf('Condición de K: %.2e\n', condest(K(d,d)));
fprintf('Condición de M: %.2e\n', condest(M(d,d)));

%% Condiciones de contorno (ejemplo: placa empotrada en todos los lados)
% Identificar nodos en los bordes (depende de la malla)
%borde = unique([LaG(:,1); LaG(:,2); LaG(:,3); LaG(:,4)]);
%d = setdiff(1:ngdl, [3*borde-2; 3*borde-1; 3*borde]); % GDL libres

%% Solución del problema de valores propios
opts.tol = 1e-8;
opts.p = 20;
[eigenvectors, eigenvalues] = eigs(K(d,d), M(d,d), 10, 'sm', opts);

frequencies = sqrt(diag(eigenvalues))/(2*pi);
disp('Frecuencias naturales (Hz):');
disp(frequencies(1:6));

%% Visualización de modos de vibración
escala = 0.5; % Factor de escala para visualización
def = zeros(ngdl,1);
def(d,1) = eigenvectors(:,3); % Primer modo

vect_mov = reshape(def,3,nno)'; % [W, ?x, ?y] por nodo

figure;
hold on;
for e = 1:nef
   fill3(xnod(LaG(e,[1 2 3 4 1]),X), ...
         xnod(LaG(e,[1 2 3 4 1]),Y), ...
         escala*vect_mov(LaG(e,[1 2 3 4 1]),ww), ...
         vect_mov(LaG(e,[1 2 3 4 1]),ww));
end
title(sprintf('Primer modo de vibración (escala x%d)',escala));
colormap jet; colorbar;
view(3); axis tight;
xlabel('X'); ylabel('Y'); zlabel('Desplazamiento W');

%% Función auxiliar para puntos de Gauss
function [x, w] = gauss_points(n)
    % Puntos y pesos para cuadratura de Gauss-Legendre
    switch n
        case 1
            x = 0; w = 2;
        case 2
            x = [-1/sqrt(3); 1/sqrt(3)];
            w = [1; 1];
        case 3
            x = [-sqrt(3/5); 0; sqrt(3/5)];
            w = [5/9; 8/9; 5/9];
        otherwise
            error('Solo implementado para n=1,2,3');
    end
end