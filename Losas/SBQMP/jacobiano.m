function [J, detJ] = jacobiano(x_nodes, y_nodes, xi, eta)
    % x_nodes, y_nodes: Coordenadas reales de los 4 nodos del elemento
    % xi, eta: Punto de integración en coordenadas naturales (ej: Gauss)
    
    % Funciones de forma bilineales estándar para mapeo geométrico
    N = @(xi, eta) 0.25 * [(1-xi)*(1-eta);  % N1
                           (1+xi)*(1-eta);  % N2
                           (1+xi)*(1+eta);  % N3
                           (1-xi)*(1+eta)]; % N4
    
    % Derivadas de las funciones de forma respecto a xi y eta
    dN_dxi = 0.25 * [-(1-eta);  % dN1/dxi
                      (1-eta);   % dN2/dxi
                      (1+eta);   % dN3/dxi
                     -(1+eta)];  % dN4/dxi
    
    dN_deta = 0.25 * [-(1-xi);  % dN1/deta
                      -(1+xi);   % dN2/deta
                       (1+xi);   % dN3/deta
                       (1-xi)];  % dN4/deta
    
    % Cálculo del Jacobiano
    dx_dxi = dN_dxi' * x_nodes;  % Suma (dN_i/dxi * x_i)
    dy_dxi = dN_dxi' * y_nodes;
    dx_deta = dN_deta' * x_nodes;
    dy_deta = dN_deta' * y_nodes;
    
    J = [dx_dxi, dy_dxi;
         dx_deta, dy_deta];
    
    detJ = det(J);
end