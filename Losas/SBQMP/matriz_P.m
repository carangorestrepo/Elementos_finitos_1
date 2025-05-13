function P = matriz_P(x, y)
    % x, y: Coordenadas del punto donde se evalúa P (pueden ser xi, eta en naturales)  
    P = [1, -x, -y, -x^2/2, -(x^2*y/2 + y^3/12), -y^2/2, ...
         -(x*y^2/2 + x^3/12), -x*y/2, x/2, x*y/2, y/2, x*y/2;
         0, 1, 0, x, x*y, 0, (x^2/4 + y^2/2), y/2, 1/2, y/2, 0, -y/2;
         0, 0, 1, 0, (y^2/4 + x^2/2), y, x*y, x/2, 0, -x/2, 1/2, x/2];
end