A = sym(gallery(5));
k=[216760,-306770,105490,-19561,4282.20000000000,-510.880000000000;-306770,668240,-475140,137940,-29375,5385.70000000000;105490,-475140,731370,-493230,159600,-29327;-19561,137940,-493230,749020,-494470,145710;4282.20000000000,-29375,159600,-494470,738110,-511900;-510.880000000000,5385.70000000000,-29327,145710,-511900,889940];
m=[256,0,0,0,0,0;0,256,0,0,0,0;0,0,256,0,0,0;0,0,0,256,0,0;0,0,0,0,256,0;0,0,0,0,0,256];
 

iteraciones=100;
modos=6;

[Phi, Omega, T] = estodola_1(k, m, iteraciones, modos);

tol=1e-10;
 max_iter=100;
 sizek = size(k, 1);

% Inicializar matriz de deflación como matriz identidad
[L, p] = chol(m, 'lower');
Linv = inv(L);
A = Linv * k * Linv';  % Matriz simétrica definida positiva
[lambda, Y] = jacobi_method(A, tol, max_iter);
modes = L' \ Y;
% Frecuencias naturales (rad/s)
omega = sqrt(abs(lambda));  % abs por seguridad numérica
% Periodos (s)
frequencies = omega;
periods = 2 * pi ./ omega;
[omega, In] = sort(omega);     % Ordena frecuencias y obtiene índices
Phi = Phi(:, In);              % Reordena modos según frecuencias
% Normalizar modos con respecto a la matriz de masa (??M? = I)
Phi = Phi ./ repmat(sqrt(diag(Phi' * m * Phi)'), size(Phi, 1), 1);