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

%% Matrices constitutivas
D_b = (E*h^3)/(12*(1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu)/2];
D_s = k * h * G * eye(2);
D = blkdiag(D_b, D_s);
T = rho * diag([h, h^3/12, h^3/12]);

%% Funciones de forma y derivadas para cuadrilátero bilineal
Nforma = @(xi,eta) [(1-xi)*(1-eta)/4; (1+xi)*(1-eta)/4; 
                   (1+xi)*(1+eta)/4; (1-xi)*(1+eta)/4];

dN_dxi = @(xi,eta) [-(1-eta)/4; (1-eta)/4; (1+eta)/4; -(1+eta)/4];
dN_deta = @(xi,eta) [-(1-xi)/4; -(1+xi)/4; (1+xi)/4; (1-xi)/4];

%% Matrices [P] y [Q]
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
    xe = xnod(LaG(e,:), X);
    ye = xnod(LaG(e,:), Y);

    %% Construcción de la matriz C evaluando P en coordenadas naturales
    C = zeros(12,12);
    natural_coords = [-1 -1; 1 -1; 1 1; -1 1]; % Coordenadas naturales de los nodos
    for i = 1:4
        C(3*(i-1)+1:3*i, :) = P(natural_coords(i,1), natural_coords(i,2));
    end
    invC = inv(C');

    K0 = zeros(12);
    M0 = zeros(12);

    for pp = 1:n_gl
        for qq = 1:n_gl
            xi_gl = x_gl(pp);
            eta_gl = x_gl(qq);

            % Funciones de forma y derivadas en coordenadas naturales
            N = Nforma(xi_gl, eta_gl);
            dNxi = dN_dxi(xi_gl, eta_gl);
            dNeta = dN_deta(xi_gl, eta_gl);

            % Evaluación del Jacobiano
            J = [dNxi'; dNeta'] * [xe'; ye']';
            detJ = det(J);

            % Evaluar P y Q en punto de Gauss (coordenadas naturales)
            P_val = P(xi_gl, eta_gl);
            Q_val = Q(xi_gl, eta_gl);

            K0 = K0 + Q_val' * D * Q_val * detJ * w_gl(pp) * w_gl(qq);
            M0 = M0 + P_val' * T * P_val * detJ * w_gl(pp) * w_gl(qq);
        end
    end

    idx = [gdl(LaG(e,1),:) gdl(LaG(e,2),:) gdl(LaG(e,3),:) gdl(LaG(e,4),:)];
    K(idx,idx) = K(idx,idx) + invC' * K0 * invC;
    M(idx,idx) = M(idx,idx) + invC' * M0 * invC;
end

%% Condiciones de contorno: definir 'd' con GDL libres
% borde = ... % definir nodos de borde y GDL restringidos
% d = setdiff(1:ngdl, [3*borde-2; 3*borde-1; 3*borde]);

%% Solución del problema de valores propios
opts.tol = 1e-8;
opts.p = 20;
[eigenvectors, eigenvalues] = eigs(K(d,d), M(d,d), 10, 'sm', opts);
frequencies = sqrt(diag(eigenvalues))/(2*pi);
disp('Frecuencias naturales (Hz):');
disp(frequencies(1:6));

%% Visualización de modos de vibración
escala = 0.5;
def = zeros(ngdl,1);
def(d,1) = eigenvectors(:,2);
vect_mov = reshape(def,3,nno)';

figure;
hold on;
for e = 1:nef
   fill3(xnod(LaG(e,[1 2 3 4 1]),X), ...
         xnod(LaG(e,[1 2 3 4 1]),Y), ...
         escala*vect_mov(LaG(e,[1 2 3 4 1]),ww), ...
         vect_mov(LaG(e,[1 2 3 4 1]),ww));
end
title(sprintf('Primer modo de vibración (escala x%.1f)',escala));
colormap jet; colorbar;
view(3); axis tight;

%% Función auxiliar para puntos de Gauss
function [x, w] = gauss_points(n)
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
