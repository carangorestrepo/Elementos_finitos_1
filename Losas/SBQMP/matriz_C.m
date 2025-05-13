function C = matriz_C(x, y)
    % Coordenadas de los 4 nodos del elemento (x_i, y_i)
    nodes = [x(1), y(1); x(2), y(2); x(3), y(3); x(4), y(4)];
    C = zeros(12, 12); % 4 nodos × 3 DOF (W, ?x, ?y)
    for i = 1:4
        xi = nodes(i, 1);
        yi = nodes(i, 2);
        
        % Submatriz C_i para el nodo i (Ecuación 14)
        Ci = [1, -xi, -yi, -xi^2/2, -(xi^2*yi/2 + yi^3/12), -yi^2/2, ...
              -(xi*yi^2/2 + xi^3/12), -xi*yi/2, xi/2, xi*yi/2, yi/2, xi*yi/2;
              0, 1, 0, xi, xi*yi, 0, (xi^2/4 + yi^2/2), yi/2, 1/2, yi/2, 0, -yi/2;
              0, 0, 1, 0, (yi^2/4 + xi^2/2), yi, xi*yi, xi/2, 0, -xi/2, 1/2, xi/2];
          
        C(3*(i-1)+1 : 3*i, :) = Ci; % Ensamblar en C global
    end
end