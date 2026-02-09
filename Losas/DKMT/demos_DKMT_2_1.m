%% ================================================================
%  ELEMENTO TRIANGULAR DKMT – FORMULACIÓN KATILI (2019)
%  Autor base: Usuario
%  Organización y documentación: ChatGPT
%
%  Descripción:
%  Implementación simbólica del operador de deformaciones para
%  placas Mindlin usando el elemento DKMT triangular.
%
%  GDL por nodo:
%     w   -> desplazamiento transversal
%     Bx  -> rotación eje x
%     By  -> rotación eje y
%
%  Total GDL elemento = 9
%
% ================================================================

clc
clear

%% ================================================================
%% 1. DEFINICIÓN DE VARIABLES SIMBÓLICAS
%% ================================================================

syms xi_gl eta_gl
syms x1 x2 x3 y1 y2 y3
syms nu h

% Parámetros geométricos lados elemento
syms L4 L5 L6
syms C4 C5 C6 S4 S5 S6
syms phi4 phi5 phi6

% Variables auxiliares
syms x21 x32 x13 y21 y32 y13

% Vectores parámetros lado
Ck = [C4;C5;C6];
Sk = [S4;S5;S6];
Lk = [L4;L5;L6];
phi_k = [phi4;phi5;phi6];

%% ================================================================
%% 2. FUNCIONES DE FORMA DEL ELEMENTO
%% ================================================================
% ---- Funciones principales ----
% (Tabla 1 Katili)

Nforma = [ 1-xi_gl-eta_gl
           xi_gl
           eta_gl ];

% Derivadas naturales
dN_dxi = [ -1
            1
            0 ];

dN_deta = [ -1
             0
             1 ];

% ---- Funciones secundarias ----
Pforma = [ 4*(1-xi_gl-eta_gl)*xi_gl
           4*xi_gl*eta_gl
           4*(1-xi_gl-eta_gl)*eta_gl ];

dP_dxi  = [ 4 - 8*xi_gl - 4*eta_gl
            4*eta_gl
           -4*eta_gl ];

dP_deta = [ -4*xi_gl
             4*xi_gl
             4 - 4*xi_gl - 8*eta_gl ];

%% ================================================================
%% 3. MATRIZ JACOBIANA
%% ================================================================
% (Transformación coordenadas naturales ? físicas)

xe = [x1;x2;x3];
ye = [y1;y2;y3];

% (Ecuación 6)
dx_dxi  = sum(dN_dxi  .* xe);
dy_dxi  = sum(dN_dxi  .* ye);

dx_deta = sum(dN_deta .* xe);
dy_deta = sum(dN_deta .* ye);

Je = [ dx_dxi   dy_dxi
       dx_deta  dy_deta ];

% Determinante Jacobiano
det_Je = det(Je);

% Inversa
inv_Je = inv(Je);

j11 = inv_Je(1,1);
j12 = inv_Je(1,2);
j21 = inv_Je(2,1);
j22 = inv_Je(2,2);

%% ================================================================
%% 4. MATRIZ DE INTERPOLACIÓN DE DESPLAZAMIENTOS
%% ================================================================
% (Ecuación 10)

N = sym(zeros(3,9));

recorrido = 1:3:12;

for i=1:3
    N(:,recorrido(i):recorrido(i)+2) = ...
        diag([Nforma(i) Nforma(i) Nforma(i)]);
end

%% ================================================================
%% 5. MATRIZ Bb_beta  (Ecuación 12)
%% ================================================================

Bb_beta = sym(zeros(3,9));

for i=1:3

    % Derivadas físicas
    dNi_dx = j11*dN_dxi(i) + j12*dN_deta(i);
    dNi_dy = j21*dN_dxi(i) + j22*dN_deta(i);

    Bb_beta(:,recorrido(i):recorrido(i)+2) = ...
        [ 0   dNi_dx   0
          0   0        dNi_dy
          0   dNi_dy   dNi_dx ];
end

%% ================================================================
%% 6. MATRIZ Bb_dbeta (Ecuación 13)
%% ================================================================

Bb_dbeta = sym(zeros(3,3));

for k=1:3

    dPk_dx = j11*dP_dxi(k) + j12*dP_deta(k);
    dPk_dy = j21*dP_dxi(k) + j22*dP_deta(k);

    Bb_dbeta(:,k) = ...
        [ dPk_dx*Ck(k)
          dPk_dy*Sk(k)
          dPk_dy*Ck(k) + dPk_dx*Sk(k) ];
end

%% ================================================================
%% 7. MATRIZ A_dbeta (Ecuación 39)
%% ================================================================

A_dbeta = diag((2/3) * Lk .* (1 + phi_k));

%% ================================================================
%% 8. RELACIÓN DEFORMACIONES DE CORTE DKMT
%% ================================================================
% Derivación de parámetro dBsk
% Basado en:
%    - Compatibilidad de cortante
%    - Integración sobre lado del elemento
%
% (Ecuaciones 22a, 32, 35b, 40b)

syms s wi wj
syms Bxi Bxj Byi Byj
syms Bni Bnj Bsi Bsj dBsk phik Ck Sk Lk

% Interpolación rotaciones tangenciales
Bsi = Ck*Bxi + Sk*Byi;
Bsj = Ck*Bxj + Sk*Byj;

Bs = (1-s/Lk)*Bsi + (s/Lk)*Bsj ...
     + 4*s/Lk*(1-s/Lk)*dBsk;

ws = (1-s/Lk)*wi + (s/Lk)*wj;

% Cortante DKMT
gamazk = -(2/3)*phik*dBsk;        % (22a)
gamask = Bs + diff(ws,s);         % (32)

% Condición compatibilidad
ec = int(gamask-gamazk,s,0,Lk) == 0;

dBsk_r = solve(ec,dBsk);          % (35b)

%% ================================================================
%% 9. MATRIZ An (Ecuaciones 41–44)
%% ================================================================

syms w1 w2 w3
syms Bx1 Bx2 Bx3 By1 By2 By3

% Definición cosenos directores lados
C4 = x21/L4;   S4 = y21/L4;
C5 = x32/L5;   S5 = y32/L5;
C6 = x13/L6;   S6 = y13/L6;

% Sustituciones
dbs4 = subs(dBsk_r,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk,phik},...
            {w2,w1,Bx2,Bx1,By2,By1,L4,C4,S4,phi4});

dbs5 = subs(dBsk_r,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk,phik},...
            {w3,w2,Bx3,Bx2,By3,By2,L5,C5,S5,phi5});

dbs6 = subs(dBsk_r,{wj,wi,Bxj,Bxi,Byj,Byi,Lk,Ck,Sk,phik},...
            {w1,w3,Bx1,Bx3,By1,By3,L6,C6,S6,phi6});

Un = [ w1 Bx1 By1 w2 Bx2 By2 w3 Bx3 By3 ].';

An = simplify(equationsToMatrix([dbs4;dbs5;dbs6],Un));

Adb = diag((2/3)*[L4*(1+phi4) L5*(1+phi5) L6*(1+phi6)]);

Aw = simplify(Adb*An);

An = A_dbeta \ Aw;

%% ================================================================
%% 10. MATRIZ DEFORMACIÓN FINAL Bb
%% ================================================================

Bb = simplify(Bb_beta + Bb_dbeta*An);

%% ================================================================
%% 11. INTERPOLACIÓN CORTANTE GLOBAL
%% ================================================================

syms gxz1 gyz1 gxz2 gyz2 gxz3 gyz3
syms gsz4 gsz5 gsz6

N1 = 1-xi_gl-eta_gl;
N2 = xi_gl;
N3 = eta_gl;

A1 = C4*S6 - C6*S4;
A2 = C5*S4 - C4*S5;
A3 = C6*S5 - C5*S6;

% Interpolación cortante nodal
gxz1 = ( S6*gsz4 - S4*gsz6)/A1;
gyz1 = (-C6*gsz4 + C4*gsz6)/A1;

gxz2 = ( S4*gsz5 - S5*gsz4)/A2;
gyz2 = (-C4*gsz5 + C5*gsz4)/A2;

gxz3 = ( S5*gsz6 - S6*gsz5)/A3;
gyz3 = (-C5*gsz6 + C6*gsz5)/A3;

gxzbar = N1*gxz1 + N2*gxz2 + N3*gxz3;
gyzbar = N1*gyz1 + N2*gyz2 + N3*gyz3;

M6 = equationsToMatrix([gxzbar;gyzbar],...
                       [gxz1;gyz1;gxz2;gyz2;gxz3;gyz3]);

M5 = equationsToMatrix([gxz1;gyz1;gxz2;gyz2;gxz3;gyz3],...
                       [gsz4;gsz5;gsz6]);

M8 = diag(-(2/3)*[phi4 phi5 phi6]);

Ngamma = M6*M5;

Bs_dbeta = Ngamma*M8;

Bs = simplify(Bs_dbeta*An);

