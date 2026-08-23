%% ========================================================================
% MATRICES SIMBÓLICAS Kb Y Ks DEL ELEMENTO MITC3+
%
% Basado en:
%   Katili, Maknun, Batoz y Katili
%   A comparative formulation of T3gamma_s, DST, DKMT and MITC3+
%
% GDL nodales:
%
%   un = [w1; Bx1; By1;
%         w2; Bx2; By2;
%         w3; Bx3; By3]
%
% Variables internas:
%
%   dbeta = [dBx7; dBy7]
%
% Matriz ampliada:
%
%   K_ampliada = [K11  K12
%                  K21  K22]
%
% Matriz condensada:
%
%   K = K11 - K12*(K22\K21)
%
% Separación condensada:
%
%   Kb = integral(Bb_cond.'*Hb*Bb_cond dA)
%   Ks = integral(Bs_cond.'*Hs*Bs_cond dA)
%
%   K = Kb + Ks
% ========================================================================

clc
clear
format compact


%% ========================================================================
% 1. VARIABLES SIMBÓLICAS
% ========================================================================

syms xi eta real

syms x1 x2 x3 real
syms y1 y2 y3 real

syms E nu h kappa real
syms d real

assume(E > 0)
assume(h > 0)
assume(kappa > 0)

assumeAlso(nu > -1)
assumeAlso(nu < 1/2)


%% ========================================================================
% 2. GRADOS DE LIBERTAD
% ========================================================================

syms w1 Bx1 By1 real
syms w2 Bx2 By2 real
syms w3 Bx3 By3 real

syms dBx7 dBy7 real

un = [w1
      Bx1
      By1
      w2
      Bx2
      By2
      w3
      Bx3
      By3];

dbeta = [dBx7
         dBy7];

q_ampliado = [un
              dbeta];


%% ========================================================================
% 3. FUNCIONES DE FORMA LINEALES
% ========================================================================

N1 = 1-xi-eta;
N2 = xi;
N3 = eta;

Nforma = [N1
          N2
          N3];

dN_dxi = diff(Nforma,xi);

dN_deta = diff(Nforma,eta);


%% ========================================================================
% 4. FUNCIÓN BURBUJA DEL MITC3+
%
%   P7 = 27*N1*N2*N3
%
% Propiedades:
%
%   P7 = 0 sobre los tres lados.
%   P7 = 1 en el centroide xi = eta = 1/3.
% ========================================================================

P7 = expand(27*N1*N2*N3);

dP7_dxi  = simplify(diff(P7,xi));
dP7_deta = simplify(diff(P7,eta));

disp('Función burbuja P7:')
pretty(P7)

disp('Derivada P7,xi:')
pretty(dP7_dxi)

disp('Derivada P7,eta:')
pretty(dP7_deta)


%% ========================================================================
% 5. INTERPOLACIÓN DE LA GEOMETRÍA
% ========================================================================

x = simplify(N1*x1 + N2*x2 + N3*x3);

y = simplify(N1*y1 + N2*y2 + N3*y3);


%% ========================================================================
% 6. MATRIZ JACOBIANA
%
% Convención:
%
%       [dx/dxi    dy/dxi ]
%   J = [                  ]
%       [dx/deta   dy/deta]
%
% Entonces:
%
%   [f,x]     [f,xi ]
%   [   ] = j [     ]
%   [f,y]     [f,eta]
%
% donde:
%
%   j = inv(J)
% ========================================================================

dx_dxi  = diff(x,xi);
dy_dxi  = diff(y,xi);

dx_deta = diff(x,eta);
dy_deta = diff(y,eta);

J = [dx_dxi,  dy_dxi
     dx_deta, dy_deta];

detJ = simplify(det(J));

j = simplify(inv(J));

j11 = j(1,1);
j12 = j(1,2);
j21 = j(2,1);
j22 = j(2,2);

disp('Matriz jacobiana:')
pretty(J)

disp('Determinante jacobiano:')
pretty(detJ)


%% ========================================================================
% 7. DERIVADAS FÍSICAS DE Ni
% ========================================================================

dN_dx = simplify( ...
         j11*dN_dxi + j12*dN_deta);

dN_dy = simplify( ...
         j21*dN_dxi + j22*dN_deta);


%% ========================================================================
% 8. OPERADOR NODAL DE FLEXIÓN Bb_beta
%
%   chi_nodal = Bb_beta*un
%
% con:
%
%   chi_x  = Bx,x
%   chi_y  = By,y
%   chi_xy = Bx,y + By,x
% ========================================================================

Bb_beta = sym(zeros(3,9));

recorrido = 1:3:9;

for i = 1:3

    columnas = recorrido(i):recorrido(i)+2;

    Bb_beta(:,columnas) = ...
        [0, dN_dx(i), 0
         0, 0,         dN_dy(i)
         0, dN_dy(i), dN_dx(i)];

end

Bb_beta = simplify(Bb_beta);


%% ========================================================================
% 9. OPERADOR INTERNO DE FLEXIÓN Bb_dbeta
%
% Rotaciones internas:
%
%   Bx_burbuja = P7*dBx7
%   By_burbuja = P7*dBy7
%
% Por tanto:
%
%                   [P7,x    0  ]
%   Bb_dbeta =      [0      P7,y]
%                   [P7,y   P7,x]
% ========================================================================

dP7_dx = simplify( ...
          j11*dP7_dxi + j12*dP7_deta);

dP7_dy = simplify( ...
          j21*dP7_dxi + j22*dP7_deta);

Bb_dbeta = [dP7_dx, 0
            0,       dP7_dy
            dP7_dy, dP7_dx];

Bb_dbeta = simplify(Bb_dbeta);


%% ========================================================================
% 10. COMPROBACIÓN DEL CAMPO DE CURVATURAS
%
%   chi = Bb_beta*un + Bb_dbeta*dbeta
% ========================================================================

chi = simplify(Bb_beta*un + Bb_dbeta*dbeta);


%% ========================================================================
% 11. MATRICES CONSTITUTIVAS
% ========================================================================

G = E/(2*(1+nu));

Db = E*h^3/(12*(1-nu^2));

Hb = Db*[1,  nu,             0
         nu, 1,              0
         0,  0, (1-nu)/2];

Ds = kappa*G*h;

Hs = Ds*eye(2);

Hb = simplify(Hb);
Hs = simplify(Hs);


%% ========================================================================
% 12. DIFERENCIAS DE COORDENADAS
%
% Convención del documento:
%
%   xji = xj-xi
%   yji = yj-yi
% ========================================================================

x21 = x2-x1;
x32 = x3-x2;
x13 = x1-x3;

y21 = y2-y1;
y32 = y3-y2;
y13 = y1-y3;


%% ========================================================================
% 13. OPERADOR NATURAL DE CORTANTE DEL MITC3
%
%   Bs_nat_MITC3 =
%
%       [Bs1  Bs2  Bs3]
%
% Cada bloque tiene dimensiones 2x3.
%
% Las filas corresponden a:
%
%   gamma_xi
%   gamma_eta
% ========================================================================

Bs1_MITC3 = [ ...
    -1,  (x21+eta*x32)/2,  (y21+eta*y32)/2
    -1, -(x13+xi*x32)/2,  -(y13+xi*y32)/2 ];

Bs2_MITC3 = [ ...
     1,  (x21+eta*x13)/2,  (y21+eta*y13)/2
     0,       -xi*x13/2,        -xi*y13/2 ];

Bs3_MITC3 = [ ...
     0,        eta*x21/2,         eta*y21/2
     1, -(x13+xi*x21)/2,  -(y13+xi*y21)/2 ];

Bs_nat_MITC3 = [Bs1_MITC3, ...
                Bs2_MITC3, ...
                Bs3_MITC3];

Bs_nat_MITC3 = simplify(Bs_nat_MITC3);


%% ========================================================================
% 14. CORRECCIÓN PLUS DEL OPERADOR NODAL
%
% Según el documento:
%
%   d_hat = 1/6 - d
%
% con:
%
%   d aproximadamente igual a 10^(-4)
%
% La matriz adicional conserva:
%
%   gamma_xi  lineal respecto a eta
%   gamma_eta lineal respecto a xi
% ========================================================================

d_hat = sym(1)/6-d;

Bs1_plus = -d_hat*[ ...
    0, 3*eta*(1-eta)*x32, 3*eta*(1-eta)*y32
    0, xi*(1-3*xi)*x32,   xi*(1-3*xi)*y32 ];

Bs2_plus = -d_hat*[ ...
    0, 3*eta*(1-eta)*x13, 3*eta*(1-eta)*y13
    0, xi*(1-3*xi)*x13,   xi*(1-3*xi)*y13 ];

Bs3_plus = -d_hat*[ ...
    0, 3*eta*(1-eta)*x21, 3*eta*(1-eta)*y21
    0, xi*(1-3*xi)*x21,   xi*(1-3*xi)*y21 ];

Bs_nat_plus = [Bs1_plus, ...
               Bs2_plus, ...
               Bs3_plus];

Bs_nat_plus = simplify(Bs_nat_plus);


%% ========================================================================
% 15. OPERADOR NATURAL NODAL DEL MITC3+
%
%   Bs_nat_beta =
%
%       Bs_nat_MITC3 + Bs_nat_plus
% ========================================================================

Bs_nat_beta = simplify( ...
              Bs_nat_MITC3 + Bs_nat_plus);


%% ========================================================================
% 16. OPERADOR NATURAL DE CORTANTE ASOCIADO A LAS VARIABLES INTERNAS
%
% Ecuación del MITC3+:
%
%                     1 [ x21    y21  ]
%   Bs_nat_dbeta = -----[             ]
%                     2 [-x13   -y13  ]
% ========================================================================

Bs_nat_dbeta = sym(1)/2*[ ...
     x21,  y21
    -x13, -y13 ];

Bs_nat_dbeta = simplify(Bs_nat_dbeta);


%% ========================================================================
% 17. TRANSFORMACIÓN A COMPONENTES CARTESIANAS
%
%   [gamma_x]       [gamma_xi ]
%   [       ] = j * [         ]
%   [gamma_y]       [gamma_eta]
%
% Por tanto:
%
%   Bs_beta  = j*Bs_nat_beta
%   Bs_dbeta = j*Bs_nat_dbeta
% ========================================================================

Bs_beta = simplify(j*Bs_nat_beta);

Bs_dbeta = simplify(j*Bs_nat_dbeta);


%% ========================================================================
% 18. CAMPO COMPLETO DE CORTANTE
%
%   gamma = Bs_beta*un + Bs_dbeta*dbeta
% ========================================================================

gamma = simplify( ...
        Bs_beta*un + Bs_dbeta*dbeta);


%% ========================================================================
% 19. INTEGRANDOS DE FLEXIÓN: MATRIZ AMPLIADA
%
%   Kb_ampliada =
%
%       [Kb11  Kb12
%        Kb21  Kb22]
% ========================================================================

integrando_Kb11 = simplify( ...
                   Bb_beta.'*Hb*Bb_beta*detJ);

integrando_Kb12 = simplify( ...
                   Bb_beta.'*Hb*Bb_dbeta*detJ);

integrando_Kb21 = simplify( ...
                   Bb_dbeta.'*Hb*Bb_beta*detJ);

integrando_Kb22 = simplify( ...
                   Bb_dbeta.'*Hb*Bb_dbeta*detJ);


%% ========================================================================
% 20. INTEGRANDOS DE CORTANTE: MATRIZ AMPLIADA
%
%   Ks_ampliada =
%
%       [Ks11  Ks12
%        Ks21  Ks22]
% ========================================================================

integrando_Ks11 = simplify( ...
                   Bs_beta.'*Hs*Bs_beta*detJ);

integrando_Ks12 = simplify( ...
                   Bs_beta.'*Hs*Bs_dbeta*detJ);

integrando_Ks21 = simplify( ...
                   Bs_dbeta.'*Hs*Bs_beta*detJ);

integrando_Ks22 = simplify( ...
                   Bs_dbeta.'*Hs*Bs_dbeta*detJ);


%% ========================================================================
% 21. INTEGRACIÓN SIMBÓLICA EXACTA
%
% Dominio natural:
%
%   0 <= xi <= 1
%   0 <= eta <= 1-xi
% ========================================================================

disp('Integrando bloques de flexión...')

Kb11 = int(int(integrando_Kb11,eta,0,1-xi),xi,0,1);

Kb12 = int(int(integrando_Kb12,eta,0,1-xi),xi,0,1);

Kb21 = int(int(integrando_Kb21,eta,0,1-xi),xi,0,1);

Kb22 = int(int(integrando_Kb22,eta,0,1-xi),xi,0,1);


disp('Integrando bloques de cortante...')

Ks11 = int(int(integrando_Ks11,eta,0,1-xi),xi,0,1);

Ks12 = int(int(integrando_Ks12,eta,0,1-xi),xi,0,1);

Ks21 = int(int(integrando_Ks21,eta,0,1-xi),xi,0,1);

Ks22 = int(int(integrando_Ks22,eta,0,1-xi),xi,0,1);


%% ========================================================================
% 22. SIMPLIFICACIÓN DE LOS BLOQUES
% ========================================================================

Kb11 = simplify(Kb11,'Steps',20);
Kb12 = simplify(Kb12,'Steps',20);
Kb21 = simplify(Kb21,'Steps',20);
Kb22 = simplify(Kb22,'Steps',20);

Ks11 = simplify(Ks11,'Steps',20);
Ks12 = simplify(Ks12,'Steps',20);
Ks21 = simplify(Ks21,'Steps',20);
Ks22 = simplify(Ks22,'Steps',20);


%% ========================================================================
% 23. COMPROBACIÓN DE ORTOGONALIDAD DE FLEXIÓN
%
% Debido a que P7 es cero sobre el contorno:
%
%   integral_A P7,x dA = 0
%   integral_A P7,y dA = 0
%
% y se debe obtener:
%
%   Kb12 = 0
%   Kb21 = 0
% ========================================================================

error_ortogonalidad_12 = simplify(Kb12);

error_ortogonalidad_21 = simplify(Kb21);


%% ========================================================================
% 24. MATRICES AMPLIADAS DE FLEXIÓN Y CORTANTE
% ========================================================================

Kb_ampliada = [Kb11, Kb12
               Kb21, Kb22];

Ks_ampliada = [Ks11, Ks12
               Ks21, Ks22];

K_ampliada = simplify(Kb_ampliada + Ks_ampliada);


%% ========================================================================
% 25. BLOQUES TOTALES PARA CONDENSACIÓN
% ========================================================================

K11 = simplify(Kb11+Ks11);

K12 = simplify(Kb12+Ks12);

K21 = simplify(Kb21+Ks21);

K22 = simplify(Kb22+Ks22);


%% ========================================================================
% 26. TRANSFORMACIÓN DE CONDENSACIÓN
%
% Equilibrio de las variables internas:
%
%   K21*un + K22*dbeta = 0
%
% Por tanto:
%
%   dbeta = Tcond*un
%
% donde:
%
%   Tcond = -(K22\K21)
% ========================================================================

Tcond = simplify(-(K22\K21),'Steps',20);


%% ========================================================================
% 27. OPERADORES CONDENSADOS
%
% Curvaturas:
%
%   chi = Bb_beta*un + Bb_dbeta*dbeta
%
% Como:
%
%   dbeta = Tcond*un
%
% entonces:
%
%   Bb_cond = Bb_beta + Bb_dbeta*Tcond
%
%
% Cortantes:
%
%   gamma = Bs_beta*un + Bs_dbeta*dbeta
%
% entonces:
%
%   Bs_cond = Bs_beta + Bs_dbeta*Tcond
% ========================================================================

Bb_cond = simplify( ...
          Bb_beta + Bb_dbeta*Tcond);

Bs_cond = simplify( ...
          Bs_beta + Bs_dbeta*Tcond);


%% ========================================================================
% 28. MATRIZ CONDENSADA DE FLEXIÓN
%
% Esta es la contribución de flexión después de aplicar a las rotaciones
% internas la misma transformación de condensación obtenida de la energía
% total.
% ========================================================================

integrando_Kb_cond = simplify( ...
                      Bb_cond.'*Hb*Bb_cond*detJ);

Kb = int( ...
     int(integrando_Kb_cond,eta,0,1-xi), ...
     xi,0,1);

Kb = simplify(Kb,'Steps',20);


%% ========================================================================
% 29. MATRIZ CONDENSADA DE CORTANTE
% ========================================================================

integrando_Ks_cond = simplify( ...
                      Bs_cond.'*Hs*Bs_cond*detJ);

Ks = int( ...
     int(integrando_Ks_cond,eta,0,1-xi), ...
     xi,0,1);

Ks = simplify(Ks,'Steps',20);


%% ========================================================================
% 30. MATRIZ TOTAL CONDENSADA
% ========================================================================

Ke = simplify(Kb+Ks,'Steps',20);


%% ========================================================================
% 31. COMPROBACIÓN MEDIANTE COMPLEMENTO DE SCHUR
%
% La matriz total condensada también debe ser:
%
%   Ke_Schur = K11 - K12*(K22\K21)
% ========================================================================

Ke_Schur = simplify( ...
           K11-K12*(K22\K21), ...
           'Steps',20);

error_condensacion = simplify(Ke-Ke_Schur);


%% ========================================================================
% 32. COMPROBACIONES DE SIMETRÍA
% ========================================================================

error_simetria_Kb = simplify(Kb-Kb.');

error_simetria_Ks = simplify(Ks-Ks.');

error_simetria_Ke = simplify(Ke-Ke.');


%% ========================================================================
% 33. RESULTADOS
% ========================================================================

disp(' ')
disp('===============================================================')
disp('Kb11: flexión nodal 9x9')
disp('===============================================================')
disp(Kb11)

disp(' ')
disp('===============================================================')
disp('Kb12: acoplamiento flexión nodal-interno 9x2')
disp('Debe ser simbólicamente cero')
disp('===============================================================')
disp(Kb12)

disp(' ')
disp('===============================================================')
disp('Kb22: flexión interna 2x2')
disp('===============================================================')
disp(Kb22)

disp(' ')
disp('===============================================================')
disp('Ks11: cortante nodal 9x9')
disp('===============================================================')
disp(Ks11)

disp(' ')
disp('===============================================================')
disp('Ks12: acoplamiento cortante nodal-interno 9x2')
disp('===============================================================')
disp(Ks12)

disp(' ')
disp('===============================================================')
disp('Ks22: cortante interno 2x2')
disp('===============================================================')
disp(Ks22)

disp(' ')
disp('===============================================================')
disp('Matriz condensada de flexión Kb, 9x9')
disp('===============================================================')
disp(Kb)

disp(' ')
disp('===============================================================')
disp('Matriz condensada de cortante Ks, 9x9')
disp('===============================================================')
disp(Ks)

disp(' ')
disp('===============================================================')
disp('Matriz total MITC3+ Ke = Kb + Ks')
disp('===============================================================')
disp(Ke)

disp(' ')
disp('Error Ke - Ke_Schur:')
disp(error_condensacion)

disp('Error de simetría de Kb:')
disp(error_simetria_Kb)

disp('Error de simetría de Ks:')
disp(error_simetria_Ks)

disp('Error de simetría de Ke:')
disp(error_simetria_Ke)

disp('Comprobación Kb12 = 0:')
disp(error_ortogonalidad_12)