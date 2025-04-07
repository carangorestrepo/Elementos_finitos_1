clear all;
close all;
clc;

% =============================================
% PARÁMETROS DEL MATERIAL Y GEOMETRÍA
% =============================================
E = 1000;               % Módulo de Young (MPa)
nu = 0.3;               % Coeficiente de Poisson
h = 0.1;                % Espesor de la placa (m)
rho = 2700;             % Densidad (kg/m³)
k = 5/6;                % Factor de corrección por cortante
G = E / (2*(1 + nu));   % Módulo de cortante (MPa)

% =============================================
% COORDENADAS NODALES (Elemento cuadrilátero distorsionado)
% =============================================
% Nodos: [x1, y1; x2, y2; x3, y3; x4, y4]
nodes = [0, 0;          % Nodo 1
         2, 0;          % Nodo 2
         3, 2;          % Nodo 3
         1, 3];         % Nodo 4

% =============================================
% FUNCIONES DE FORMA Y DERIVADAS (Cuadrilátero bilineal)
% =============================================
syms xi eta;
N = [(1 - xi)*(1 - eta)/4;   % N1
     (1 + xi)*(1 - eta)/4;   % N2
     (1 + xi)*(1 + eta)/4;   % N3
     (1 - xi)*(1 + eta)/4];  % N4

% Derivadas de las funciones de forma
dN_dxi = diff(N, xi);
dN_deta = diff(N, eta);

% =============================================
% MATRIZ JACOBIANA Y DETERMINANTE
% =============================================
% Coordenadas cartesianas interpoladas
x = sum(N .* nodes(:, 1));
y = sum(N .* nodes(:, 2));

% Derivadas de (x, y) respecto a (xi, eta)
dx_dxi = sum(dN_dxi .* nodes(:, 1));
dy_dxi = sum(dN_dxi .* nodes(:, 2));
dx_deta = sum(dN_deta .* nodes(:, 1));
dy_deta = sum(dN_deta .* nodes(:, 2));

% Matriz Jacobiana
J = [dx_dxi, dy_dxi;
     dx_deta, dy_deta];
detJ = det(J);

% Mostrar Jacobiano en el centro del elemento (xi = 0, eta = 0)
J_center = double(subs(J, [xi, eta], [0, 0]));
detJ_center = double(subs(detJ, [xi, eta], [0, 0]));
disp('Jacobiano en (?=0, ?=0):');
disp(J_center);
disp(['Determinante: ', num2str(detJ_center)]);

% =============================================
% MATRICES DE RIGIDEZ Y MASA
% =============================================
% Matriz de rigidez de flexión [D_b]
D_b = (E*h^3)/(12*(1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu)/2];

% Matriz de rigidez de cortante [D_s]
D_s = k * h * G * eye(2);

% Matriz de rigidez total [D]
D = blkdiag(D_b, D_s);

% Matriz de densidad de masa [T]
T = rho * diag([h, h^3/12, h^3/12]);

% Matrices [P] y [Q] (definidas en el documento)
P =  @(xi,eta)[1, -xi, -eta, -xi^2/2, -(xi^2*eta/2 + eta^3/12), -eta^2/2, -(xi*eta^2/2 + xi^3/12), -xi*eta/2, xi/2, xi*eta/2, eta/2, xi*eta/2;
               0, 1, 0, xi, xi*eta, 0, (xi^2/4 + eta^2/2), eta/2, 1/2, eta/2, 0, -eta/2;
               0, 0, 1, 0, (eta^2/4 + xi^2/2), eta, xi*eta, xi/2, 0, -xi/2, 1/2, xi/2];

Q = @(xi,eta)[0, 0, 0, 1, eta, 0, xi/2, 0, 0, 0, 0, 0;
               0, 0, 0, 0, 0, 0, 0, 1, xi, eta/2, 0, 0;
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2*xi;
               0, 0, 0, 0, 0, 0, 0, 0, 1, eta, 0, 0;
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, xi];

% Matriz [C] (relación entre desplazamientos nodales y parámetros ?)
% Se evalúa [P] en los nodos del elemento
C = zeros(12, 12);
for i = 1:4
    xi_val = nodes(i, 1); % Usar coordenadas naturales (no es correcto, se debe mapear a [-1,1])
    eta_val = nodes(i, 2); % Esto es un placeholder. En la práctica, se usan las coordenadas naturales de los nodos.
    % (Aquí se necesita un mapeo adecuado)
    C(3*(i-1)+1:3*i, :) = P(xi_val, eta_val);
    %C(3*(i-1)+1:3*i, :) = double(subs(P, [xi, eta], [xi_val, eta_val]));
end

% Integración numérica (Gauss 2x2)
n_gl=2;
[gauss_pts, gauss_wts] = gaussQuad(2);
[xi_gl, w_gl] = gausslegendre_quad(n_gl);
K0 = zeros(12, 12);
M0 = zeros(12, 12);

for i = 1:length(gauss_pts)
    xi_gl = gauss_pts(i, 1);
    eta_gl = gauss_pts(i, 2);
    
    % Evaluar Jacobiano y su determinante en el punto de Gauss
    J_val = double(subs(J, {xi, eta},{xi_gl, eta_gl}));
    detJ_val = double(subs(detJ, {xi, eta},{xi_gl, eta_gl}));
    % Evaluar [P] y [Q] en el punto de Gauss
    P_val = P(xi_gl, eta_gl);
    Q_val = Q(xi_gl, eta_gl);
    
    % Contribución a las matrices K0 y M0
    K0 = K0 + Q_val' * D * Q_val * detJ_val * gauss_wts(i);
    M0 = M0 + P_val' * T * P_val * detJ_val * gauss_wts(i);
end

% Matrices de rigidez y masa del elemento
Ke = inv(C') * K0 * inv(C);
Me = inv(C') * M0 * inv(C);

% =============================================
% VIBRACIONES LIBRES (Eigenvalue problem)
% =============================================
% Ensamblar matrices globales (aquí solo para un elemento)
K_global = Ke;
M_global = Me;

% Resolver el problema de autovalores: (K - ?²M)? = 0
[eigenvectors, eigenvalues] = eig(K_global, M_global);
frequencies = sqrt(diag(eigenvalues)) / (2*pi); % Convertir a Hz

disp('Frecuencias naturales (Hz):');
disp(frequencies(1:6)); % Mostrar las primeras 6 frecuencias

% =============================================
% FUNCIÓN AUXILIAR: Cuadratura de Gauss
% =============================================
function [pts, wts] = gaussQuad(n)
    [xi, w] = lgwt(n, -1, 1);
    [xi_grid, eta_grid] = meshgrid(xi, xi);
    pts = [xi_grid(:), eta_grid(:)];
    wts = kron(w, w)';
end

function [x, w] = lgwt(n, a, b)
    [x, w] = gaussLegendre(n);
    x = 0.5*(b - a)*x + 0.5*(b + a);
    w = 0.5*(b - a)*w;
end

function [x, w] = gaussLegendre(n)
    % Puntos y pesos para Gauss-Legendre (implementación simplificada)
    switch n
        case 1
            x = 0;
            w = 2;
        case 2
            x = [-1/sqrt(3); 1/sqrt(3)];
            w = [1; 1];
        case 3
            x = [-sqrt(3/5); 0; sqrt(3/5)];
            w = [5/9; 8/9; 5/9];
        otherwise
            error('Solo se implementa para n=1, 2, 3');
    end
end