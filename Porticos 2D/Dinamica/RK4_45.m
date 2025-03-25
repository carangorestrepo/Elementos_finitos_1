% Solución de EDO de segundo orden usando Runge-Kutta-Fehlberg (RKF45)
% Basado en la Figura 15-1, página 287 del libro de Luis Enrique García

clear  % Limpia el espacio de trabajo
clc    % Limpia la ventana de comandos

% Definimos las condiciones iniciales
y = [0;0]; % y(0) = 0, y'(0) = 0
% el centro
h = 0.02;   % Tamaño de paso
x0 = 0;     % Punto inicial
N = 2650;   % Número de iteraciones 

% Prealocamos espacio para almacenar los resultados
%% Ejemplo 14-2 pag 250 LEG
%h = 0.005;   % Tamaño de paso
%x0 = 0;     % Punto inicial
%N = 121;   % Número de iteraciones % Ejemplo 14-2 pag 250 LEG
ys = [y, zeros(2, N)];

% Bucle principal para la integración numérica
for i = 1:N
    % Coeficientes del método RKF45
    a2 = 0.25; a3 = 0.375; a4 = 12/13; a6 = 0.5;
    b21 = 0.25;
    b31 = 3/32; b32 = 9/32;
    b41 = 1932/2197; b42 = -7200/2197; b43 = 7296/2197;
    b51 = 439/216; b52 = -8; b53 = 3680/513; b54 = -845/4104;
    b61 = -8/27; b62 = 2; b63 = -3544/2565; b64 = 1859/4104; b65 = -11/40;
    c1 = 25/216; c3 = 1408/2565; c4 = 2197/4104; c5 = -0.20;
    d1 = 1/360; d3 = -128/4275; d4 = -2197/75240; d5 = 0.02; d6 = 2/55;
    
    % Calculamos los incrementos h ponderados por los coeficientes
    h2 = a2 * h; h3 = a3 * h; h4 = a4 * h; h6 = a6 * h;
    
    % Cálculo de los k's (pendientes intermedias)
    k1 = F(x0, y, i);
    k2 = F(x0 + h2, y + h * b21 * k1, i);
    k3 = F(x0 + h3, y + h * (b31 * k1 + b32 * k2), i);
    k4 = F(x0 + h4, y + h * (b41 * k1 + b42 * k2 + b43 * k3), i);
    k5 = F(x0 + h, y + h * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4), i);
    k6 = F(x0 + h6, y + h * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5), i);
    
    % Cálculo de la nueva estimación de y
    y_nueva = y + h * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5);
    
    % Cálculo del error estimado (aunque no se usa aquí para ajustar h)
    out = d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6;
    
    % Actualizamos el estado para la siguiente iteración
    y = y_nueva;
    ys(:, i + 1) = y; % Guardamos los resultados
    x0 = x0 + h;      % Avanzamos al siguiente punto
end

% Definimos el vector de tiempo
t = 0:h:N * h;

% Graficamos la respuesta del sistema
plot(t, -ys(1,:)); % Se grafica la respuesta del desplazamiento

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definimos la función que describe el sistema dinámico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = F(t, m, i)
    z = m(1, 1); % Desplazamiento
    y = m(2, 1); % Velocidad
    
    % Parámetros del sistema
    mi = 1;  % Masa
    xi = 0.05; % Amortiguamiento % el centro
    %xi = 2/100; % Amortiguamiento % el centro
    T = 0.075;     % Período
    w = 2 * pi / T; % Frecuencia
    k = mi * w^2; % Constante de rigidez
    C = 2 * mi * w * xi; % Coeficiente de amortiguamiento
    
    % Término independiente (fuerza externa)
    a = el_centro;
    %a=ejerccio_14_2LEG;
    ind = a(i); % Fuerza externa dependiente del tiempo
    
    % Cálculo de la EDO
    yp = C / mi * y;  % Término de velocidad
    ys = z * k / mi;   % Término de desplazamiento
    res = [y; ind - yp - ys]; % Devuelve el vector con velocidad y aceleración
end
