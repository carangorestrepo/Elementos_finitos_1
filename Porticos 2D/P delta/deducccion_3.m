clc
clear
syms EI P ksGA x k C1 C2 C3 C4 C5 C6 xi

%% PROPIEDADES DEL MATERIAL Y GEOMETRÍA
E = 24870062.3;      % Módulo de elasticidad [kN/m²]
G = 0.4*E;         % Módulo de cortante [kN/m²]
I = 0.4^4/12;      % Inercia de la sección [m?]
Ae = 0.4^2;        % Área efectiva [m²]
EI_val = E*I;      % Rigidez a flexión [kN·m²]
ksGA_val = Ae*G*5/6; % Rigidez a cortante reducida [kN]

% Carga axial aplicada
P_val = 1000;      % Carga axial [kN]

% Longitud del elemento
L = 4;             % Longitud del elemento [m]
AE_val = Ae*E;     % Rigidez axial [kN]

% Sustitución de valores numéricos
EI = EI_val;
ksGA = ksGA_val;
P = P_val;
AE = AE_val;

% Factor ? (corregido según teoría Timoshenko)
chi = 1/(1 - P/ksGA);  % CORRECCIÓN IMPORTANTE

% Parámetro k (según ecuaciones del documento)
k_val = sqrt(P/(EI*chi));  % = ?

% Solución homogénea para desplazamiento vertical v(x)
% Según documento: v(x) = C1 + C2*x + C3*cos(kx) + C4*sin(kx)
v = C1 + C2*x + C3*cos(k_val*x) + C4*sin(k_val*x);

% Derivadas
dv = diff(v, x);
d2v = diff(v, x, 2);
d3v = diff(v, x, 3);

% Rotación ?(x) - CORRECCIÓN CRÍTICA según Timoshenko
% ? = dv/dx - Q/(ksGA) pero Q = -EI·d³v/dx³ para solución homogénea
theta = dv - (-EI*d3v)/ksGA;

% Simplificamos ?
theta = simplify(theta);

% Momento M (según Timoshenko)
M = -EI * diff(theta, x);

% Cortante Q 
Q = -EI * diff(theta, x, 2);

% Fuerza cortante total V = Q - P*dv/dx
V = Q - P*dv;

% Solución para desplazamiento axial u(x)
% du/dx = P/AE (para carga axial constante)
b = 0; % Sin carga distribuida axial
A = int(b, x) + C5; % Fuerza axial = constante = P
u = int(A/AE, x) + C6;

%% MATRIZ DE RIGIDEZ
K_TE2 = zeros(6);   % Inicialización matriz de rigidez tangente
N_u2 = sym(zeros(1,6));
N_v2 = sym(zeros(1,6)); % Funciones de forma para desplazamiento vertical
N_t2 = sym(zeros(1,6)); % Funciones de forma para rotación

for i = 1:6
    % Resolver constantes para cada condición de contorno unitaria
    [c1,c2,c3,c4,c5,c6] = solve(...
        subs(u,x,0) == (i==1),...    % u1
        subs(v,x,0) == (i==2),...    % v1  
        subs(theta,x,0) == (i==3),... % ?1
        subs(u,x,L) == (i==4),...    % u2
        subs(v,x,L) == (i==5),...    % v2
        subs(theta,x,L) == (i==6),... % ?2
        [C1,C2,C3,C4,C5,C6]);
    
    % Evaluar fuerzas en los extremos (CON SIGNOS CORRECTOS)
    % Convención: fuerzas positivas en dirección positiva de ejes
    K_TE2(:,i) = [
        -subs(A,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0});     % Fx1
        -subs(V,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0});     % Fy1  
        -subs(M,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0});     % M1
         subs(A,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L});     % Fx2
         subs(V,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L});     % Fy2
        -subs(M,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L})      % M2
    ];
    
    % Funciones de forma
    N_t2(i) = subs(theta,{C1,C2,C3,C4,C5,C6},{c1,c2,c3,c4,c5,c6});
    N_v2(i) = subs(v,{C1,C2,C3,C4,C5,C6},{c1,c2,c3,c4,c5,c6});
    N_u2(i) = subs(u,{C1,C2,C3,C4,C5,C6},{c1,c2,c3,c4,c5,c6});
end

% Convertir a numérico
K_TE2 = double(K_TE2);

%% VERIFICACIÓN
fprintf('Matriz de rigidez (6x6):\n');
disp(K_TE2);

fprintf('Determinante: %.6e\n', det(K_TE2));
fprintf('Traza: %.6f\n', trace(K_TE2));

% Verificar simetría
error_simetria = norm(K_TE2 - K_TE2', 'fro')/norm(K_TE2, 'fro');
fprintf('Error de simetría: %.2e\n', error_simetria);

% Derivadas
Bv = diff(N_v2, x);
Bt = diff(N_t2, x);
Bu = diff(N_u2, x);



% Flexión
Kb = double(int(EI * (Bt.' * Bt), x, 0, L));

% Cortante
Bs = (Bv - N_t2);
Ks = double(int(ksGA * (Bs.' * Bs), x, 0, L));

% Axial
Ka =double( int(AE * (Bu.' * Bu), x, 0, L));

% Geométrica
Kg = double(int(P * (Bv.' * Bv), x, 0, L));


K_total1 = Ka + Kb + Ks - Kg;

% Verificación
residuo=K_TE2 - (Ka+Kb+Ks-Kg);

