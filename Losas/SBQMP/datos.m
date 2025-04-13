% Coordenadas reales del elemento (puede ser distorsionado)
xy = [0 0;
      2 0;
      2.5 1.5;
      0 2];

E = 2e11;     % Pa
nu = 0.3;
h = 0.1;      % m

Ke = SBQMP_K_matrix(xy, E, nu, h);
disp('Matriz de rigidez Ke:');
disp(Ke);