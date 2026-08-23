clc
clear
format long g

%% ========================================================================
% VIGA DE TIMOSHENKO + TORSION NO UNIFORME + ALABEO DE VLASOV
%
% Se obtiene:
%
%   1. Matriz de rigidez exacta del elemento.
%
%   2. Funciones de forma exactas:
%
%          Nu(x)    -> desplazamiento axial
%          Nv(x)    -> desplazamiento transversal
%          Nt(x)    -> giro por flexion
%          Nphi(x)  -> giro torsional
%          Npsi(x)  -> amplitud de alabeo psi = d(phi)/dx
%
%   3. Fuerzas nodales equivalentes para:
%
%          qx(x) -> carga axial distribuida
%          qy(x) -> carga transversal distribuida
%          mx(x) -> momento torsor distribuido
%
%   4. Formulacion matricial:
%
%          d = B*C + dp
%          f = S*C + fp
%
%      de donde:
%
%          K = S/B
%
%      y:
%
%          f = K*d + fp - K*dp
%
%      Si se adopta:
%
%          fint = K*d - feq
%
%      entonces:
%
%          feq = K*dp - fp
%
% ------------------------------------------------------------------------
%
% Grados de libertad locales:
%
%      d = [u1;
%           v1;
%           t1;
%           phi1;
%           psi1;
%           u2;
%           v2;
%           t2;
%           phi2;
%           psi2]
%
% donde:
%
%      phi = giro torsional
%
%      psi = d(phi)/dx
%
%            coordenada generalizada asociada al alabeo
%
% ------------------------------------------------------------------------
%
% Estado interno:
%
%      s = [N;
%           V;
%           M;
%           T;
%           Bw]
%
% donde:
%
%      N  = fuerza axial
%      V  = fuerza cortante
%      M  = momento flector
%      T  = momento torsor total
%      Bw = bimomento
%
% ------------------------------------------------------------------------
%
% ECUACIONES DIFERENCIALES ADOPTADAS
%
% FLEXION TIMOSHENKO
%
%      N' = qx
%      V' = qy
%      M' = V
%
%      u' = N/EA
%      t' = M/EI
%      v' = t - V/Ks
%
% TORSION + ALABEO DE VLASOV
%
%      T'   = mx
%
%      phi' = psi
%
%      psi' = Bw/EIw
%
%      Bw'  = GJ*psi - T
%
% Estas dos ultimas expresiones conducen a:
%
%      Bw = EIw*phi''
%
%      T  = GJ*phi' - EIw*phi'''
%
% y por tanto, para mx = 0:
%
%      EIw*phi'''' - GJ*phi'' = 0
%
% donde:
%
%      Ks  = kappa*G*A
%      GJ  = rigidez torsional de Saint-Venant
%      EIw = E*Iw = rigidez de alabeo
%
% ========================================================================


%% ========================================================================
% DEFINICION DE VARIABLES SIMBOLICAS
% ========================================================================
syms x
syms C1 C2 C3 C4 C5 C6
syms C7 C8 C9 C10
%% ========================================================================
% PROPIEDADES DEL MATERIAL Y LA SECCION
% ========================================================================
g = 9.8066502;             % [m/s^2] aceleracion de la gravedad
E = 24870062.3;            % [kPa] modulo de elasticidad
G = 0.4*E;                 % [kPa] modulo de cortante
L = 4;                     % [m] longitud del elemento
Ae = 0.4^2;                % [m^2] area
I = 0.4^4/12;              % [m^4] inercia flexional
%% ------------------------------------------------------------------------
% PROPIEDADES DE TORSION Y ALABEO
%
% IMPORTANTE:
%
% J e Iw deben provenir de las propiedades reales de la seccion.
%
% Los valores siguientes son solamente valores de ejemplo para que
% el codigo pueda ejecutarse.
% -------------------------------------------------------------------------
J  = 0.001;                % [m^4] constante torsional Saint-Venant
Iw = 1.0e-5;               % [m^6] constante de alabeo
%% ========================================================================
% RIGIDECES
% ========================================================================
EA = E*Ae;                 % [kN]
EI = E*I;                  % [kN*m^2]
kappa = 5/6;
Ks = kappa*G*Ae;           % [kN]
GJ = G*J;                  % [kN*m^2]
EIw = E*Iw;                % [kN*m^4]
%% ========================================================================
% PARAMETRO CARACTERISTICO DE TORSION CON ALABEO
% ========================================================================
lambda_w = sqrt(GJ/EIw);
%% ========================================================================
% PARTE I
%
% SOLUCION HOMOGENEA
%
% Se utiliza para obtener:
%
%   - matriz de rigidez
%   - funciones de forma
%
% Para ello:
%
%       qx = 0
%       qy = 0
%       mx = 0
%
% ========================================================================
qx = sym(0);
qy = sym(0);
mx = sym(0);
%% ========================================================================
% ECUACIONES DIFERENCIALES HOMOGENEAS
% ========================================================================
%% ------------------------------------------------------------------------
% FUERZA CORTANTE
%
%       V' = qy
%
% con qy = 0:
%
%       V = C1
% -------------------------------------------------------------------------
V = int(qy,x) + C1;
%% ------------------------------------------------------------------------
% MOMENTO FLECTOR
%
%       M' = V
% -------------------------------------------------------------------------
M = int(V,x) + C2;
%% ------------------------------------------------------------------------
% GIRO POR FLEXION
%
%       t' = M/EI
% -------------------------------------------------------------------------
t = int(M/EI,x) + C3;
%% ------------------------------------------------------------------------
% DESPLAZAMIENTO TRANSVERSAL
%
% Timoshenko:
%
%       v' = t - V/Ks
% -------------------------------------------------------------------------
v = int(t - V/Ks,x) + C4;
%% ------------------------------------------------------------------------
% FUERZA AXIAL
%
%       N' = qx
%
% con qx = 0:
%
%       N = C5
% -------------------------------------------------------------------------
N = int(qx,x) + C5;
%% ------------------------------------------------------------------------
% DESPLAZAMIENTO AXIAL
%
%       u' = N/EA
% -------------------------------------------------------------------------
u = int(N/EA,x) + C6;
%% ========================================================================
% TORSION CON ALABEO
% ========================================================================
%
% Para mx = 0:
%
%       EIw*phi'''' - GJ*phi'' = 0
%
% cuya solucion general exacta es:
%
%       phi =
%           C7
%         + C8*x
%         + C9*cosh(lambda_w*x)
%         + C10*sinh(lambda_w*x)
%
% ========================================================================
phi = C7 ...
    + C8*x ...
    + C9*cosh(lambda_w*x) ...
    + C10*sinh(lambda_w*x);

%% ------------------------------------------------------------------------
% AMPLITUD DE ALABEO
%
%       psi = phi'
% -------------------------------------------------------------------------
psi = diff(phi,x);
%% ------------------------------------------------------------------------
% BIMOMENTO
%
%       Bw = EIw*phi''
% -------------------------------------------------------------------------
Bw = EIw*diff(phi,x,2);
%% ------------------------------------------------------------------------
% MOMENTO TORSOR TOTAL
%
%       T = GJ*phi' - EIw*phi'''
% -------------------------------------------------------------------------
T = GJ*diff(phi,x) - EIw*diff(phi,x,3);
%% ========================================================================
% SIMPLIFICACION SIMBOLICA
% ========================================================================
N   = simplify(N);
V   = simplify(V);
M   = simplify(M);
u   = simplify(u);
v   = simplify(v);
t   = simplify(t);
phi = simplify(phi);
psi = simplify(psi);
T   = simplify(T);
Bw  = simplify(Bw);
%% ========================================================================
% COMPROBACION DE LAS ECUACIONES DE TORSION-ALABEO
% ========================================================================

error_phi_1 = simplify(diff(phi,x) - psi);

error_phi_2 = simplify(diff(psi,x) - Bw/EIw);

error_phi_3 = simplify(diff(Bw,x) - (GJ*psi - T));

error_phi_4 = simplify(diff(T,x) - mx);


%% ========================================================================
% MATRIZ DE RIGIDEZ MEDIANTE DESPLAZAMIENTOS UNITARIOS
% ========================================================================
FNE = sym(zeros(10,1));
K_TE2 = sym(zeros(10,10));
%% ========================================================================
% FUNCIONES DE FORMA
% ========================================================================
Nu   = sym(zeros(1,10));
Nv   = sym(zeros(1,10));
Nt   = sym(zeros(1,10));
Nphi = sym(zeros(1,10));
Npsi = sym(zeros(1,10));
%% ========================================================================
% CONSTANTES DE INTEGRACION
% ========================================================================
Constantes = [C1 C2 C3 C4 C5 C6 C7 C8 C9 C10];
%% ------------------------------------------------------------------------
% Cada iteracion impone uno de los diez desplazamientos nodales igual a 1
% y los otros nueve iguales a cero.
%
% La respuesta correspondiente genera una columna de la matriz Ke.
% -------------------------------------------------------------------------
for i = 1:10
    [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10] = solve( ...
        subs(u,x,0)   == (i==1), ...
        subs(v,x,0)   == (i==2), ...
        subs(t,x,0)   == (i==3), ...
        subs(phi,x,0) == (i==4), ...
        subs(psi,x,0) == (i==5), ...
        subs(u,x,L)   == (i==6), ...
        subs(v,x,L)   == (i==7), ...
        subs(t,x,L)   == (i==8), ...
        subs(phi,x,L) == (i==9), ...
        subs(psi,x,L) == (i==10), ...
        [C1,C2,C3,C4,C5,C6,C7,C8,C9,C10]);
    ci = {c1,c2,c3,c4,c5,c6,c7,c8,c9,c10};
    %% ====================================================================
    % FUNCIONES DE FORMA
    %
    % Cada columna corresponde a la respuesta generada por un GDL unitario.
    %
    % Por tanto:
    %
    %       u(x)   = Nu(x)*d
    %
    %       v(x)   = Nv(x)*d
    %
    %       t(x)   = Nt(x)*d
    %
    %       phi(x) = Nphi(x)*d
    %
    %       psi(x) = Npsi(x)*d
    %
    % ====================================================================
    Nu(i) = subs(u,Constantes,ci);
    Nv(i) = subs(v,Constantes,ci);
    Nt(i) = subs(t,Constantes,ci);
    Nphi(i) = subs(phi,Constantes,ci);
    Npsi(i) = subs(psi,Constantes,ci);
    %% ====================================================================
    % FUERZAS NODALES
    %
    % Convencion utilizada:
    %
    %       f =
    %
    %       [-N(0);
    %         V(0);
    %        -M(0);
    %        -T(0);
    %        -Bw(0);
    %
    %         N(L);
    %        -V(L);
    %         M(L);
    %         T(L);
    %         Bw(L)]
    %
    % ====================================================================
    K_TE2(:,i) = [ ...
        -subs(N, ...
        [Constantes x], ...
        [ci {0}]);
         subs(V, ...
        [Constantes x], ...
        [ci {0}]);
        -subs(M, ...
        [Constantes x], ...
        [ci {0}]);
        -subs(T, ...
        [Constantes x], ...
        [ci {0}]);
        -subs(Bw, ...
        [Constantes x], ...
        [ci {0}]);
         subs(N, ...
        [Constantes x], ...
        [ci {L}]);
        -subs(V, ...
        [Constantes x], ...
        [ci {L}]);
         subs(M, ...
        [Constantes x], ...
        [ci {L}]);
         subs(T, ...
        [Constantes x], ...
        [ci {L}]);
         subs(Bw, ...
        [Constantes x], ...
        [ci {L}])
        ];
end
%% ========================================================================
% SIMPLIFICACION DE FUNCIONES DE FORMA
% ========================================================================
Nu = simplify(Nu);
Nv = simplify(Nv);
Nt = simplify(Nt);
Nphi = simplify(Nphi);
Npsi = simplify(Npsi);
%% ========================================================================
% COMPROBACION CINEMATICA DEL ALABEO
%
% Debe cumplirse:
%
%       Npsi = d(Nphi)/dx
% ========================================================================
error_Npsi = simplify(Npsi - diff(Nphi,x));


%% ========================================================================
% SIMPLIFICACION DE MATRIZ DE RIGIDEZ
% ========================================================================
K_TE2 = simplify(K_TE2);
K_TE2 = double(K_TE2);
%% ========================================================================
% SIMETRIZACION NUMERICA
% ========================================================================

K_TE2 = 0.5*(K_TE2 + K_TE2.');


%% ========================================================================
% MOSTRAR MATRIZ DE RIGIDEZ
% ========================================================================

disp(' ')
disp('============================================================')
disp(' MATRIZ DE RIGIDEZ TIMOSHENKO + TORSION + ALABEO')
disp('============================================================')

disp(K_TE2)


%% ========================================================================
% PARTE II
%
% CARGAS DISTRIBUIDAS TRAPEZOIDALES
% ========================================================================


%% ------------------------------------------------------------------------
% CARGA AXIAL
% -------------------------------------------------------------------------
ba = 25;                   % [kN/m] inicial
bb = 35;                   % [kN/m] final
%% ------------------------------------------------------------------------
% CARGA TRANSVERSAL
% -------------------------------------------------------------------------
wa = 25;                   % [kN/m] inicial
wb = 30;                   % [kN/m] final
%% ------------------------------------------------------------------------
% MOMENTO TORSOR DISTRIBUIDO
% -------------------------------------------------------------------------
ma = 5;                    % [kN*m/m] inicial
mb = 8;                    % [kN*m/m] final
%% ========================================================================
% CARGAS TRAPEZOIDALES
% ========================================================================
b = ba + (bb-ba)*x/L;
w = wa + (wb-wa)*x/L;
mt = ma + (mb-ma)*x/L;
%% ========================================================================
% SOLUCION PARTICULAR TRANSVERSAL
% ========================================================================
%% Vq
%
%       Vq' = w
%
Vq = int(w,x);
%% Mq
%
%       Mq' = Vq
%
Mq = int(Vq,x);
%% tq
%
%       tq' = Mq/EI
%
tq = int(Mq/EI,x);
%% vq
%
%       vq' = tq - Vq/Ks
%
vq = int(tq - Vq/Ks,x);
%% ========================================================================
% SOLUCION PARTICULAR AXIAL
% ========================================================================


%% Nq
%
%       Nq' = b
%
Nq = int(b,x);
%% uq
%
%       uq' = Nq/EA
%
uq = int(Nq/EA,x);

%% ========================================================================
% SOLUCION PARTICULAR TORSION-ALABEO
%
% Se resuelve:
%
%       Tq'   = mt
%
%       Bwq'  = GJ*psiq - Tq
%
%       psiq' = Bwq/EIw
%
%       phiq' = psiq
%
% Para construir una solucion particular se toman constantes de
% integracion iguales a cero.
%
% ========================================================================

syms Tqf(x) Bwqf(x) psiqf(x) phiqf(x)


%% ------------------------------------------------------------------------
% Para evitar que dsolve introduzca constantes arbitrarias que luego
% interfieran con C7...C10, se obtiene una solucion particular imponiendo
% condiciones iniciales nulas.
% -------------------------------------------------------------------------

Tq = int(mt,x,0,x);


%% ------------------------------------------------------------------------
% El subsistema puede escribirse:
%
%       EIw*phiq'''' - GJ*phiq'' = -mt
%
% porque:
%
%       T = GJ*phi' - EIw*phi'''
%
% y:
%
%       T' = mt
%
% luego:
%
%       GJ*phi'' - EIw*phi'''' = mt
%
%       EIw*phi'''' - GJ*phi'' = -mt
%
% -------------------------------------------------------------------------
ode_phi_q = EIw*diff(phiqf,x,4) ...
          - GJ*diff(phiqf,x,2) ...
          == -mt;
%% ------------------------------------------------------------------------
% Condiciones iniciales nulas para obtener solamente una solucion
% particular:
%
%       phiq(0)   = 0
%       phiq'(0)  = 0
%       phiq''(0) = 0
%       phiq'''(0)= 0
%
% -------------------------------------------------------------------------
cond_phi_q = [ ...
    phiqf(0) == 0, ...
    subs(diff(phiqf,x),x,0) == 0, ...
    subs(diff(phiqf,x,2),x,0) == 0, ...
    subs(diff(phiqf,x,3),x,0) == 0];
phiq_sol = dsolve(ode_phi_q,cond_phi_q);
%% ------------------------------------------------------------------------
% Convertir a expresion simbolica
% -------------------------------------------------------------------------
phiq = simplify(phiq_sol);
%% ------------------------------------------------------------------------
% Variables particulares asociadas
% -------------------------------------------------------------------------
psiq = simplify(diff(phiq,x));
Bwq = simplify(EIw*diff(phiq,x,2));
Tq = simplify(GJ*diff(phiq,x) ...
            - EIw*diff(phiq,x,3));
%% ========================================================================
% SIMPLIFICACION
% ========================================================================
Vq = simplify(Vq);
Mq = simplify(Mq);
tq = simplify(tq);
vq = simplify(vq);
Nq = simplify(Nq);
uq = simplify(uq);
phiq = simplify(phiq);
psiq = simplify(psiq);
Bwq = simplify(Bwq);
Tq = simplify(Tq);
%% ========================================================================
% COMPROBACION DE SOLUCION PARTICULAR TORSIONAL
% ========================================================================
error_Tq = simplify(diff(Tq,x) - mt);
%% ========================================================================
% FUERZAS NODALES OBTENIDAS DIRECTAMENTE IMPONIENDO d = 0
%
% Se impone:
%
%       u_total   = 0
%       v_total   = 0
%       t_total   = 0
%       phi_total = 0
%       psi_total = 0
%
% en ambos extremos.
%
% La solucion total es:
%
%       u_total   = u   + uq
%       v_total   = v   + vq
%       t_total   = t   + tq
%       phi_total = phi + phiq
%       psi_total = psi + psiq
%
% ========================================================================
[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10] = solve( ...
    subs(u+uq,x,0)       == 0, ...
    subs(v+vq,x,0)       == 0, ...
    subs(t+tq,x,0)       == 0, ...
    subs(phi+phiq,x,0)   == 0, ...
    subs(psi+psiq,x,0)   == 0, ...
    subs(u+uq,x,L)       == 0, ...
    subs(v+vq,x,L)       == 0, ...
    subs(t+tq,x,L)       == 0, ...
    subs(phi+phiq,x,L)   == 0, ...
    subs(psi+psiq,x,L)   == 0, ...
    [C1,C2,C3,C4,C5,C6,C7,C8,C9,C10]);
ci = {c1,c2,c3,c4,c5,c6,c7,c8,c9,c10};
%% ========================================================================
% FUERZAS DE EMPOTRAMIENTO PERFECTO
%
% Convencion:
%
%       [-N1;
%         V1;
%        -M1;
%        -T1;
%        -Bw1;
%
%         N2;
%        -V2;
%         M2;
%         T2;
%         Bw2]
%
% ========================================================================
FNE(:) = [ ...
    -subs(Nq+N, ...
    [Constantes x], ...
    [ci {0}]);
     subs(Vq+V, ...
    [Constantes x], ...
    [ci {0}]);
    -subs(Mq+M, ...
    [Constantes x], ...
    [ci {0}]);
    -subs(Tq+T, ...
    [Constantes x], ...
    [ci {0}]);
    -subs(Bwq+Bw, ...
    [Constantes x], ...
    [ci {0}]);
     subs(Nq+N, ...
    [Constantes x], ...
    [ci {L}]);
    -subs(Vq+V, ...
    [Constantes x], ...
    [ci {L}]);
     subs(Mq+M, ...
    [Constantes x], ...
    [ci {L}]);
     subs(Tq+T, ...
    [Constantes x], ...
    [ci {L}]);
     subs(Bwq+Bw, ...
    [Constantes x], ...
    [ci {L}])
    ];

FNE = double(FNE);

disp(' ')
disp('============================================================')
disp(' FUERZAS DE EMPOTRAMIENTO PERFECTO')
disp('============================================================')
disp(FNE)

%% ========================================================================
% PARTE III
%
% FORMULACION MATRICIAL GENERAL
%
% Se escribe:
%
%       d = B*C + dp
%
%       f = S*C + fp
%
% donde:
%
%       C  = constantes de integracion
%
%       dp = desplazamientos particulares
%
%       fp = fuerzas particulares
%
% ========================================================================
Un = [C1;
      C2;
      C3;
      C4;
      C5;
      C6;
      C7;
      C8;
      C9;
      C10];
%% ========================================================================
% MATRIZ B
%
% Relacion:
%
%       d = B*C
%
% para la solucion homogenea.
% ========================================================================
Bd = equationsToMatrix( ...
    [u;
     v;
     t;
     phi;
     psi], ...
     Un);
%% Nodo inicial
B0 = simplify(subs(Bd,x,0));
%% Nodo final
BL = simplify(subs(Bd,x,L));
%% Matriz completa
B = [B0;
     BL];
B = double(B);
%% ========================================================================
% MATRIZ S
%
% Estado interno:
%
%       s = [N;
%            V;
%            M;
%            T;
%            Bw]
%
% ========================================================================
Sf = equationsToMatrix( ...
    [N;
     V;
     M;
     T;
     Bw], ...
     Un);
%% Nodo inicial

S0 = simplify(subs(Sf,x,0));


%% Nodo final

SL = simplify(subs(Sf,x,L));


%% ========================================================================
% MATRICES DE SIGNOS
%
% Nodo inicial:
%
%       [-N;
%         V;
%        -M;
%        -T;
%        -Bw]
%
% Nodo final:
%
%       [ N;
%        -V;
%         M;
%         T;
%         Bw]
%
% ========================================================================
Ri = [-1  0  0  0  0;
       0  1  0  0  0;
       0  0 -1  0  0;
       0  0  0 -1  0;
       0  0  0  0 -1];


Rj = [ 1  0  0  0  0;
       0 -1  0  0  0;
       0  0  1  0  0;
       0  0  0  1  0;
       0  0  0  0  1];


%% Matriz S completa

S = [Ri*S0;
     Rj*SL];

S = double(S);


%% ========================================================================
% MATRIZ DE RIGIDEZ
%
% De:
%
%       d = B*C
%
%       C = inv(B)*d
%
% y:
%
%       f = S*C
%
% resulta:
%
%       f = S*inv(B)*d
%
% por tanto:
%
%       K = S/B
%
% ========================================================================

K = S/B;


%% Elimina pequeños errores numericos de asimetria

K = 0.5*(K+K.');


%% ========================================================================
% COMPARACION DE AMBAS MATRICES
% ========================================================================

error_K = max(max(abs(K-K_TE2)));

fprintf('\nError maximo entre K y K_TE2 = %.12e\n',error_K);


%% ========================================================================
% VECTOR PARTICULAR DE DESPLAZAMIENTOS
%
% La solucion total puede escribirse:
%
%       d = B*C + dp
%
% con:
%
%       dp =
%
%       [uq(0);
%        vq(0);
%        tq(0);
%        phiq(0);
%        psiq(0);
%
%        uq(L);
%        vq(L);
%        tq(L);
%        phiq(L);
%        psiq(L)]
%
% ========================================================================
dp0 = [subs(uq,x,0);
       subs(vq,x,0);
       subs(tq,x,0);
       subs(phiq,x,0);
       subs(psiq,x,0)];
dpL = [subs(uq,x,L);
       subs(vq,x,L);
       subs(tq,x,L);
       subs(phiq,x,L);
       subs(psiq,x,L)];
dp = [dp0;
      dpL];
dp = double(dp);
disp(' ')
disp('============================================================')
disp(' VECTOR PARTICULAR DE DESPLAZAMIENTOS dp')
disp('============================================================')

disp(dp)

%% ========================================================================
% VECTOR PARTICULAR DE FUERZAS
%
% El estado particular es:
%
%       sp = [Nq;
%             Vq;
%             Mq;
%             Tq;
%             Bwq]
%
% ========================================================================
sp0 = [subs(Nq,x,0);
       subs(Vq,x,0);
       subs(Mq,x,0);
       subs(Tq,x,0);
       subs(Bwq,x,0)];
spL = [subs(Nq,x,L);
       subs(Vq,x,L);
       subs(Mq,x,L);
       subs(Tq,x,L);
       subs(Bwq,x,L)];

%% Aplicar la convencion de signos

fp = [Ri*sp0;
      Rj*spL];
fp = double(fp);
%% ========================================================================
% FUERZAS DE EMPOTRAMIENTO PERFECTO
%
% Tenemos:
%
%       d = B*C + dp
%
% por tanto:
%
%       C = B\(d-dp)
%
% y:
%
%       f = S*C + fp
%
% Sustituyendo:
%
%       f = S/B*d - S/B*dp + fp
%
%       f = K*d + fp - K*dp
%
% Si:
%
%       d = 0
%
% entonces:
%
%       f0 = fp - K*dp
%
% ========================================================================
f0 = fp - K*dp;
%% ========================================================================
% FUERZAS NODALES EQUIVALENTES
%
% Si se adopta:
%
%       fint = K*d - feq
%
% comparando con:
%
%       f = K*d + fp - K*dp
%
% se obtiene:
%
%       -feq = fp - K*dp
%
% por tanto:
%
%       feq = K*dp - fp
%
% y:
%
%       feq = -f0
%
% ========================================================================

feq = K*dp - fp;
disp(' ')
disp('============================================================')
disp(' FUERZAS NODALES EQUIVALENTES')
disp('============================================================')

disp(feq)
%% ========================================================================
% COMPROBACION
%
% La solucion directa FNE corresponde a:
%
%       FNE = f0
%
% ========================================================================

error_f0 = max(abs(FNE-f0));

fprintf('\nError maximo FNE directa vs f0 = %.12e\n',error_f0);


%% ========================================================================
% COMPROBACION
%
% Como:
%
%       feq = -f0
%
% debe cumplirse:
%
%       feq + FNE = 0
%
% ========================================================================

error_feq = max(abs(feq+FNE));

fprintf('\nError maximo feq + FNE = %.12e\n',error_feq);


%% ========================================================================
% PARTE IV
%
% MATRIZ CINEMATICA BDEF
%
% A partir de las funciones de forma exactas:
%
%       u   = Nu*d
%
%       v   = Nv*d
%
%       t   = Nt*d
%
%       phi = Nphi*d
%
%       psi = Npsi*d
%
% Las deformaciones generalizadas son:
%
%       epsilon = du/dx
%
%       chi     = dt/dx
%
%       gamma   = t - dv/dx
%
%       gamma_t = dphi/dx = psi
%
%       chi_w   = dpsi/dx = phi''
%
% Por tanto:
%
%       e = Bdef*d
%
% ========================================================================
Baxial = diff(Nu,x);
Bflexion = diff(Nt,x);
Bcortante = Nt - diff(Nv,x);
Btorsion = diff(Nphi,x);
Balabeo = diff(Npsi,x);
Bdef = [Baxial;
        Bflexion;
        Bcortante;
        Btorsion;
        Balabeo];
Bdef = simplify(Bdef);
%% ========================================================================
% MATRIZ CONSTITUTIVA
%
%       D =
%
%       [EA    0    0    0     0
%
%         0   EI    0    0     0
%
%         0    0   Ks    0     0
%
%         0    0    0   GJ     0
%
%         0    0    0    0    EIw]
%
% ========================================================================

D = [EA  0   0   0    0;
      0  EI   0   0    0;
      0   0  Ks   0    0;
      0   0   0  GJ    0;
      0   0   0   0   EIw];
%% ========================================================================
% COMPROBACION DE LA MATRIZ DE RIGIDEZ MEDIANTE
%
%       K = integral(Bdef' * D * Bdef dx)
%
% ========================================================================

K_B = int(Bdef.'*D*Bdef,x,0,L);

K_B = double(K_B);

K_B = 0.5*(K_B+K_B.');


disp(' ')
disp('============================================================')
disp(' MATRIZ K OBTENIDA MEDIANTE INTEGRAL B''*D*B')
disp('============================================================')

disp(K_B)


error_KB = max(max(abs(K_B-K)));

fprintf('\nError maximo K_B vs K = %.12e\n',error_KB);


%% ========================================================================
% PARTE V
%
% FUERZAS NODALES MEDIANTE INTEGRACION DE FUNCIONES DE FORMA
%
% Para las cargas distribuidas:
%
%       feq_N = integral(Nu'   * b  dx)
%
%       feq_V = integral(Nv'   * qy dx)
%
%       feq_T = integral(Nphi' * mt dx)
%
% En el caso del momento torsor distribuido mx, el trabajo virtual es:
%
%       delta W = integral(mx * delta(phi) dx)
%
% por tanto:
%
%       feq_T = integral(Nphi' * mx dx)
%
% ========================================================================


%% ------------------------------------------------------------------------
% Fuerza axial
% -------------------------------------------------------------------------

feq_N = int(Nu.'*b,x,0,L);


%% ------------------------------------------------------------------------
% Fuerza transversal
%
% Se conservan las dos alternativas de signo del codigo original para
% comprobar la convencion.
% -------------------------------------------------------------------------

feq_V_pos = int(Nv.'*w,x,0,L);

feq_V_neg = int(Nv.'*(-w),x,0,L);


%% ------------------------------------------------------------------------
% Momento torsor distribuido
% -------------------------------------------------------------------------

feq_T_pos = int(Nphi.'*mt,x,0,L);

feq_T_neg = int(Nphi.'*(-mt),x,0,L);


%% Convertir a numerico
feq_N = double(feq_N);
feq_V_pos = double(feq_V_pos);
feq_V_neg = double(feq_V_neg);
feq_T_pos = double(feq_T_pos);
feq_T_neg = double(feq_T_neg);
%% ========================================================================
% COMBINACIONES DE SIGNO
% ========================================================================
feq_int_1 = feq_N + feq_V_pos + feq_T_pos;
feq_int_2 = feq_N + feq_V_pos + feq_T_neg;
feq_int_3 = feq_N + feq_V_neg + feq_T_pos;
feq_int_4 = feq_N + feq_V_neg + feq_T_neg;
%% ========================================================================
% COMPARACION AUTOMATICA CON feq OBTENIDO POR ECUACIONES DIFERENCIALES
% ========================================================================
errores_integracion = [ ...
    max(abs(feq_int_1-feq));
    max(abs(feq_int_2-feq));
    max(abs(feq_int_3-feq));
    max(abs(feq_int_4-feq))];
[error_minimo,caso_signos] = min(errores_integracion);
disp(' ')
disp('============================================================')
disp(' COMPARACION FNE POR INTEGRACION DE FUNCIONES DE FORMA')
disp('============================================================')
fprintf('Caso de signos que reproduce feq = %d\n',caso_signos);
fprintf('Error maximo = %.12e\n',error_minimo);
%% ========================================================================
% PARTE VI
%
% IMPRESION DE BLOQUES DE LA MATRIZ
%
% Esto permite visualizar por separado:
%
%       axial
%       flexion Timoshenko
%       torsion-alabeo
%
% ========================================================================
gdof_axial = [1 6];
gdof_flexion = [2 3 7 8];
gdof_warping = [4 5 9 10];
K_axial = K(gdof_axial,gdof_axial);
K_flexion = K(gdof_flexion,gdof_flexion);
K_warping = K(gdof_warping,gdof_warping);
disp(' ')
disp('============================================================')
disp(' BLOQUE AXIAL')
disp('============================================================')
disp(K_axial)

disp(' ')
disp('============================================================')
disp(' BLOQUE FLEXION TIMOSHENKO')
disp('============================================================')

disp(K_flexion)


disp(' ')
disp('============================================================')
disp(' BLOQUE TORSION + ALABEO')
disp('============================================================')

disp(K_warping)


%% ========================================================================
% PARTE VII
%
% VERIFICACION DE PROPIEDADES DE LAS FUNCIONES DE FORMA
% ========================================================================


%% ------------------------------------------------------------------------
% Valores nodales de Nphi
% -------------------------------------------------------------------------

Nphi_0 = simplify(subs(Nphi,x,0));

Nphi_L = simplify(subs(Nphi,x,L));


%% ------------------------------------------------------------------------
% Valores nodales de Npsi
% -------------------------------------------------------------------------

Npsi_0 = simplify(subs(Npsi,x,0));

Npsi_L = simplify(subs(Npsi,x,L));


disp(' ')
disp('============================================================')
disp(' Nphi EN LOS NODOS')
disp('============================================================')

disp(double(Nphi_0))
disp(double(Nphi_L))


disp(' ')
disp('============================================================')
disp(' Npsi EN LOS NODOS')
disp('============================================================')

disp(double(Npsi_0))
disp(double(Npsi_L))


%% ========================================================================
% PARTE VIII
%
% RESULTADOS PRINCIPALES
% ========================================================================

disp(' ')
disp('============================================================')
disp(' MATRIZ DE RIGIDEZ FINAL K')
disp('============================================================')

disp(K)


disp(' ')
disp('============================================================')
disp(' VECTOR DE FUERZAS NODALES EQUIVALENTES feq')
disp('============================================================')

disp(feq)


disp(' ')
disp('============================================================')
disp(' VECTOR DE FUERZAS DE EMPOTRAMIENTO f0')
disp('============================================================')

disp(f0)


%% ========================================================================
% FIN
% ========================================================================