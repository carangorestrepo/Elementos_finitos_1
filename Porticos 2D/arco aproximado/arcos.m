%==========================================================================
% LOCKING-FREE CURVED BEAM ELEMENT (Quarter Ring)
% Basado en Lee & Sin / SBQMP curved beam formulation
%
% Implementa:
%  - Matriz de rigidez (Ecs. 22f, 23e, 25, 27, 31)
%  - Vector de cargas distribuidas radiales (Ec. 32)
%  - Transformación energética libre de locking (Ec. 28g)
%==========================================================================

clear; clc;

%% ------------------------------------------------------------------------
% 1. PROPIEDADES GEOMÉTRICAS Y DE MATERIAL
% ------------------------------------------------------------------------
R = 100;                % Radio del arco
E = 1e7;                % Módulo de Young
G = E/2;                % Módulo de corte (suposición del artículo)
k = 5/6;                % Factor de corrección por corte

width = 1.0;            % Ancho de la sección
hR = 1/1000;            % Relación espesor / radio
h  = hR * R;            % Espesor

% Geometría del arco (quarter ring)
An = pi/2;              % Ángulo central
L  = R * An;            % Longitud del arco

A  = width * h;         % Área
I  = (width*h^3)/12;    % Momento de inercia
Ac = k * G * A;         % Rigidez efectiva a corte

% Definición de cargas
q0  = -10;  qL  = -20;    % radial
qt0 =  5;   qtL =  10;    % tangencial
m0  =  2;                 % momento

%% ------------------------------------------------------------------------
% 2. PARÁMETROS DEL MODELO (Ecs. 22f y 23e)
% ------------------------------------------------------------------------
alpha = (E*I)/Ac;
beta  = (I*R^2)/A;

lambda = (4/L^2) * (alpha + I/A + R^2);   % Ec. (22f)
mu     = (4/L^2) * (alpha + R^2);         % Ec. (23e)

%% ------------------------------------------------------------------------
% 3. FUNCIONES CINEMÁTICAS EN LOS EXTREMOS (Ecs. 25a–25c)
% ------------------------------------------------------------------------

% ---------- s = 0 (xi = 0)
xi = 0;

Hth0 = L * [ ...
    xi - 1.5*xi^2 + (2/3)*xi^3, ...
   -0.5*xi^2 + (2/3)*xi^3, ...
    2*xi^2 - (4/3)*xi^3 ];

Hwp0 = R^2 * [ ...
    1 - lambda - 3*xi + 2*xi^2, ...
   -lambda - xi + 2*xi^2, ...
    2*lambda + 4*xi - 4*xi^2 ];

Hu0 = R*L * [ ...
    0.75*mu + (1-mu)*xi - 1.5*xi^2 + (2/3)*xi^3, ...
    0.25*mu - mu*xi - 0.5*xi^2 + (2/3)*xi^3, ...
   -mu + 2*mu*xi + 2*xi^2 - (4/3)*xi^3 ];

% ---------- s = L (xi = 1)
xi = 1;

Hth1 = L * [ ...
    xi - 1.5*xi^2 + (2/3)*xi^3, ...
   -0.5*xi^2 + (2/3)*xi^3, ...
    2*xi^2 - (4/3)*xi^3 ];

Hwp1 = R^2 * [ ...
    1 - lambda - 3*xi + 2*xi^2, ...
   -lambda - xi + 2*xi^2, ...
    2*lambda + 4*xi - 4*xi^2 ];

Hu1 = R*L * [ ...
    0.75*mu + (1-mu)*xi - 1.5*xi^2 + (2/3)*xi^3, ...
    0.25*mu - mu*xi - 0.5*xi^2 + (2/3)*xi^3, ...
   -mu + 2*mu*xi + 2*xi^2 - (4/3)*xi^3 ];

%% ------------------------------------------------------------------------
% 4. MATRIZ Tk (Compatibilidad – Ecs. 25a–25c)
% ------------------------------------------------------------------------
row1 = Hwp1 - Hwp0*cos(L/R) + Hu0*sin(L/R);   % Ec. (25a)
row2 = Hu1  - Hwp0*sin(L/R) - Hu0*cos(L/R);   % Ec. (25b)
row3 = Hth1;                                  % Ec. (25c)

Tk = [row1; row2; row3];

%% ------------------------------------------------------------------------
% 5. MATRIZ Tu (Relación con DOFs físicos)
% ------------------------------------------------------------------------
Tu = zeros(3,6);

% Nodo 1
Tu(1,1) = -cos(L/R);   Tu(1,2) =  sin(L/R);   Tu(1,3) = -R*sin(L/R);
Tu(2,1) = -sin(L/R);   Tu(2,2) = -cos(L/R);   Tu(2,3) = -R*(1-cos(L/R));
Tu(3,3) = -1;

% Nodo 2
Tu(1,4) = 1;
Tu(2,5) = 1;
Tu(3,6) = 1;

%% Transformación energética (Ec. 28g implícita)
T = Tk \ Tu;

%% ------------------------------------------------------------------------
% 6. INTEGRACIÓN NUMÉRICA – RIGIDEZ Y CARGAS (Ecs. 31 y 32)
% ------------------------------------------------------------------------
ngp = 4;
[xg, wg] = gauss_points(ngp);

K_core = zeros(3,3);
R_core = zeros(3,1);



for i = 1:ngp
    xi_g = xg(i);
    s  = (L/2)*(xi_g + 1);
    xi = s/L;
    J  = L/2;

    % Funciones internas (Ec. 27)
    % --- Hk (flexión radial)
    Hk = [1 - 3*xi + 2*xi^2;
         -xi + 2*xi^2;
          4*xi - 4*xi^2];

    % --- Hu (tangencial)
    Hu = R*L * [ ...
        0.75*mu + (1-mu)*xi - 1.5*xi^2 + (2/3)*xi^3;
        0.25*mu - mu*xi - 0.5*xi^2 + (2/3)*xi^3;
       -mu + 2*mu*xi + 2*xi^2 - (4/3)*xi^3 ];

    % --- Hphi (rotación)
    Hphi = L * [ ...
        xi - 1.5*xi^2 + (2/3)*xi^3;
       -0.5*xi^2 + (2/3)*xi^3;
        2*xi^2 - (4/3)*xi^3 ];
    
    Bk = (1/L)*[-3 + 4*xi, -1 + 4*xi, 4 - 8*xi];
    Qk = (1/L^2)*[4, 4, -8];

    % Rigidez interna (Ec. 31)
    K_core = K_core + ...
        (Hk*Hk' + alpha*(Bk'*Bk) + beta*(Qk'*Qk)) * J * wg(i);

    % Cargas
    qr = q0  + (qL  - q0 )*xi;
    qt = qt0 + (qtL - qt0)*xi;

    % Ensamble Ec. (32)
    R_core = R_core ...
           + Hk   * qr * J * wg(i) ...
           + Hu   * qt * J * wg(i) ...
           + Hphi * m0 * J * wg(i);
end

K_core = E*I*K_core;

%% ------------------------------------------------------------------------
% 7. MATRIZ DE RIGIDEZ Y VECTOR DE CARGAS DEL ELEMENTO
% ------------------------------------------------------------------------
Ke = T.' * K_core * T;    % Ec. (30)
Re = T.' * R_core;        % Ec. (32)

%% RESULTADOS
disp('Matriz de rigidez Ke ='); disp(Ke)
disp('Vector de cargas Re =');  disp(Re)

%% ------------------------------------------------------------------------
% FUNCIÓN AUXILIAR: PUNTOS DE GAUSS
% ------------------------------------------------------------------------
function [x, w] = gauss_points(n)
    if n == 4
        x = [-sqrt((3+2*sqrt(6/5))/7);
             -sqrt((3-2*sqrt(6/5))/7);
              sqrt((3-2*sqrt(6/5))/7);
              sqrt((3+2*sqrt(6/5))/7)];
        w = [(18-sqrt(30))/36;
             (18+sqrt(30))/36;
             (18+sqrt(30))/36;
             (18-sqrt(30))/36];
    else
        x = [-1/sqrt(3); 1/sqrt(3)];
        w = [1; 1];
    end
end