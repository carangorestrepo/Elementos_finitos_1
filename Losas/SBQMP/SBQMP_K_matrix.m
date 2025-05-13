function Ke = SBQMP_K_matrix(xy, E, nu, h)

% xy : coordenadas reales (4x2) [x1 y1; x2 y2; x3 y3; x4 y4]
% E  : módulo de Young
% nu : coeficiente de Poisson
% h  : espesor

% --- Propiedades del material ---
k = 5/6;
Db = (E*h^3)/(12*(1 - nu^2)) * [1 nu 0; nu 1 0; 0 0 (1 - nu)/2];
G = E/(2*(1 + nu));
Ds = k * h * G * eye(2);
D = blkdiag(Db, Ds);  % 5x5 matriz constitutiva

% --- Matriz de transformación C (evaluada en nodos reales) ---
C = compute_C(xy);  % 12x12

% --- Puntos de Gauss y pesos (2x2) ---
g = [-1 1]/sqrt(3);
w = [1 1];
Ke_local = zeros(12,12);

% --- Bucle de integración numérica 2x2 ---
for i = 1:2
    for j = 1:2
        xi = g(i); eta = g(j);
        weight = w(i)*w(j);

        % --- Funciones de forma y derivadas ---
        N = 1/4 * [(1 - xi)*(1 - eta);
                   (1 + xi)*(1 - eta);
                   (1 + xi)*(1 + eta);
                   (1 - xi)*(1 + eta)];

        dN_dxi = 1/4 * [
            -(1 - eta), -(1 - xi);
             (1 - eta), -(1 + xi);
             (1 + eta),  (1 + xi);
            -(1 + eta),  (1 - xi)
        ];

        % --- Jacobiano ---
        J = dN_dxi' * xy;
        detJ = det(J);
        [x_gp, y_gp] = shape_to_physical(N, xy);  % (x,y) en punto de Gauss

        % --- Matriz Q evaluada en (x, y) ---
        Q = make_Q(x_gp, y_gp);  % 5x12

        % --- Matriz B ---
        B = Q * inv(C);

        % --- Aportación al Ke ---
        Ke_local = Ke_local + B' * D * B * detJ * weight;
    end
end

% --- Transformación final como en el artículo ---
Ke = C' \ Ke_local / C;

end

%% --- Subfunción: matriz Q(x, y) ---
function Q = make_Q(x, y)
    Q = zeros(5,12);
    Q(1,4)  = 1;
    Q(1,5)  = y;
    Q(1,6)  = x^2;

    Q(2,6)  = 1;
    Q(2,7)  = x;
    Q(2,5)  = y^2;

    Q(3,5)  = 2*x;
    Q(3,7)  = 2*y;
    Q(3,8)  = 1;

    Q(4,9)  = 1;
    Q(4,10) = y;

    Q(5,11) = 1;
    Q(5,12) = x;
end

%% --- Subfunción: coordenadas físicas en punto de Gauss ---
function [x, y] = shape_to_physical(N, xy)
x = sum(N .* xy(:,1));
y = sum(N .* xy(:,2));
end

%% --- Subfunción: matriz de transformación C (ec. 14) ---
function C = compute_C(coords)
C = zeros(12,12);
    for i = 1:4
        xi = coords(i,1);
        yi = coords(i,2);
        i1 = 3*(i-1) + 1;
        C_i = [
            1, -xi, -yi, -xi^2/2, -xi^2*yi/2 + yi^3/12, -yi^2/2, -xi*yi^2/2 + xi^3/12, -xi*yi/2, xi/2, xi*yi/2, yi^2, xi*yi^2;
            0, 1, 0, xi, xi*yi, xi^2/4 + yi^2/2, yi^2, 1/2, yi^2, 0, -yi^2, 0;
            0, 0, 1, yi^2/4 + xi^2/2, yi, xi*yi, xi^2, -xi^2, 0, -xi^2, 1/2, xi^2
        ];
        C(i1:i1+2, :) = C_i;
    end
end
