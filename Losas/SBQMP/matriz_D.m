function [D_b,D_s] = matriz_D(E, nu, h, k)
    % k: Factor de corrección por cortante (ej. k = 5/6)
    G = E / (2*(1 + nu));
    D_b = (E*h^3)/(12*(1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    D_s = k*h*G * eye(2);
    %D = blkdiag(D_b, D_s); % Ensamblar matriz bloque-diagonal
end