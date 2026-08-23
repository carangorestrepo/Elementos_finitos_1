clc
clear
format long g
%% ========================================================================
% VIGA DE TIMOSHENKO
% Se obtiene:
%   1. Matriz de rigidez exacta del elemento.
%   2. Funciones de forma exactas:
%          Nu(x) -> desplazamiento axial
%          Nv(x) -> desplazamiento transversal
%          Nt(x) -> giro
%   3. Fuerzas nodales equivalentes para cargas trapezoidales:
%          qx(x) -> carga axial
%          qy(x) -> carga transversal
%   4. Formulación matricial:
%          d = B*C + dp
%          f = S*C + fp
%      de donde:
%          K = S/B
%      y:
%          f = K*d + fp - K*dp
%      Si se adopta:
%          fint = K*d - feq
%      entonces:
%          feq = K*dp - fp
% ------------------------------------------------------------------------
% Grados de libertad locales:
%      d = [u1;
%           v1;
%           t1;
%           u2;
%           v2;
%           t2]
% ------------------------------------------------------------------------
% Estado interno:
%      s = [N;
%           V;
%           M]
% ------------------------------------------------------------------------
% Ecuaciones diferenciales adoptadas:
%      N' = qx
%      V' = qy
%      M' = V
%      u' = N/EA
%      t' = M/EI
%      v' = t - V/Ks
% donde:
%      Ks = kappa*G*A
%
% ========================================================================
%% ========================================================================
% DEFINICIÓN DE VARIABLES SIMBÓLICAS
% ========================================================================
syms x
syms C1 C2 C3 C4 C5 C6
%% ========================================================================
% PROPIEDADES DEL MATERIAL Y LA SECCIÓN
% ========================================================================
g = 9.8066502;             % [m/s^2] aceleración de la gravedad
E = 24870062.3;            % [kPa] módulo de elasticidad del concreto
G = 0.4*E;                 % [kPa] módulo de cortante
L = 4;                     % [m] longitud del elemento
Ae = 0.4^2;                % [m^2] área de la sección
I = 0.4^4/12;              % [m^4] momento de inercia
%% Rigideces
EA = E*Ae;                 % [kN] rigidez axial
EI = E*I;                  % [kN*m^2] rigidez a flexión
kappa = 5/6;               % factor de corrección por cortante
Ks = kappa*G*Ae;           % [kN] rigidez efectiva a cortante kGA
%% ========================================================================
% PARTE I
% SOLUCIÓN HOMOGÉNEA
% Se utiliza para obtener:
%   - matriz de rigidez
%   - funciones de forma
% Para ello:
%       qx = 0
%       qy = 0
% ========================================================================
qx = 0;
qy = 0;
%% ========================================================================
% ECUACIONES DIFERENCIALES HOMOGÉNEAS
% ========================================================================
%% Fuerza cortante
%       V' = qy
% con qy = 0:
%       V = C1
V = int(qy,x) + C1;
%% Momento flector
%       M' = V
M = int(V,x) + C2;
%% Giro
%       t' = M/EI
t = int(M/EI,x) + C3;
%% Desplazamiento transversal
% Timoshenko:
%       v' = t - V/Ks
v = int(t - V/Ks,x) + C4;
%% Fuerza axial
%       N' = qx
% con qx = 0:
%       N = C5
N = int(qx,x) + C5;
%% Desplazamiento axial
%       u' = N/EA
u = int(N/EA,x) + C6;
%% Simplificación simbólica
N = simplify(N);
V = simplify(V);
M = simplify(M);
u = simplify(u);
v = simplify(v);
t = simplify(t);
%% ========================================================================
% MATRIZ DE RIGIDEZ MEDIANTE DESPLAZAMIENTOS UNITARIOS
% ========================================================================
FNE = sym(zeros(6,1));
K_TE2 = sym(zeros(6,6));
Nu = sym(zeros(1,6));
Nv = sym(zeros(1,6));
Nt = sym(zeros(1,6));
%% ------------------------------------------------------------------------
% Cada iteración impone uno de los seis desplazamientos nodales igual a 1
% y los otros cinco iguales a cero.
% La respuesta correspondiente genera una columna de la matriz Ke.
% -----------------------------------------------------------------------
for i = 1:6
    [c1,c2,c3,c4,c5,c6] = solve( ...
        subs(u,x,0) == (i==1), ...
        subs(v,x,0) == (i==2), ...
        subs(t,x,0) == (i==3), ...
        subs(u,x,L) == (i==4), ...
        subs(v,x,L) == (i==5), ...
        subs(t,x,L) == (i==6), ...
        [C1,C2,C3,C4,C5,C6]);
    %% ====================================================================
    % FUNCIONES DE FORMA
    % Cada columna corresponde a la respuesta generada por un GDL unitario.
    %    % Por tanto:
    %       u(x) = Nu(x)*d
    %       v(x) = Nv(x)*d
    %       t(x) = Nt(x)*d
    % ====================================================================
    Nt(i) = subs(t, {C1,C2,C3,C4,C5,C6},{c1,c2,c3,c4,c5,c6});
    Nv(i) = subs(v, {C1,C2,C3,C4,C5,C6},{c1,c2,c3,c4,c5,c6});
    Nu(i) = subs(u, {C1,C2,C3,C4,C5,C6},{c1,c2,c3,c4,c5,c6});
    %% ====================================================================
    % FUERZAS NODALES
    % Convención utilizada:
    %       f =
    %       [-N(0);
    %         V(0);
    %        -M(0);
    %         N(L);
    %        -V(L);
    %         M(L)]
    % ====================================================================
    K_TE2(:,i) = [ ...
        -subs(N, ...
        {C1,C2,C3,C4,C5,C6,x}, ...
        {c1,c2,c3,c4,c5,c6,0});
         subs(V, ...
        {C1,C2,C3,C4,C5,C6,x}, ...
        {c1,c2,c3,c4,c5,c6,0});
        -subs(M, ...
        {C1,C2,C3,C4,C5,C6,x}, ...
        {c1,c2,c3,c4,c5,c6,0});
         subs(N, ...
        {C1,C2,C3,C4,C5,C6,x}, ...
        {c1,c2,c3,c4,c5,c6,L});
        -subs(V, ...
        {C1,C2,C3,C4,C5,C6,x}, ...
        {c1,c2,c3,c4,c5,c6,L});
         subs(M, ...
        {C1,C2,C3,C4,C5,C6,x}, ...
        {c1,c2,c3,c4,c5,c6,L})
        ];
end
%% Simplificación de funciones de forma
Nu = simplify(Nu);
Nv = simplify(Nv);
Nt = simplify(Nt);
%% Simplificación de matriz de rigidez
K_TE2 = double(simplify(K_TE2));
%% Simetrización numérica
K_TE2 = 0.5*(K_TE2 + K_TE2.');
%% ========================================================================
% PARTE II
% CARGAS DISTRIBUIDAS TRAPEZOIDALES
% ========================================================================
ba = 25;                   % [kN/m] carga axial inicial
bb = 35;                   % [kN/m] carga axial final
wa = 25;                   % [kN/m] carga transversal inicial
wb = 30;                   % [kN/m] carga transversal final
%% Carga axial trapezoidal
b = ba + (bb-ba)*x/L;
%% Carga transversal trapezoidal
w = wa + (wb-wa)*x/L;
%% ========================================================================
% SOLUCIÓN PARTICULAR TRANSVERSAL
% ========================================================================
%% Vq
%       Vq' = w
Vq = int(w,x);
%% Mq
%       Mq' = Vq
Mq = int(Vq,x);
%% tq
%       tq' = Mq/EI
tq = int(Mq/EI,x);
%% vq
%       vq' = tq - Vq/Ks
vq = int(tq - Vq/Ks,x);
%% ========================================================================
% SOLUCIÓN PARTICULAR AXIAL
% ========================================================================
%% Nq
%       Nq' = b
Nq = int(b,x);
%% uq
%       uq' = Nq/EA
uq = int(Nq/EA,x);
%% Simplificación
Vq = simplify(Vq);
Mq = simplify(Mq);
tq = simplify(tq);
vq = simplify(vq);
Nq = simplify(Nq);
uq = simplify(uq);
%% ========================================================================
% FUERZAS NODALES OBTENIDAS DIRECTAMENTE IMPONIENDO d = 0
% Se impone:
%       u_total = 0
%       v_total = 0
%       t_total = 0
% en ambos extremos.
% La solución total es:
%       u_total = u + uq
%       v_total = v + vq
%       t_total = t + tq
% ========================================================================
[c1,c2,c3,c4,c5,c6] = solve( ...
    subs(u+uq,x,0) == 0, ...
    subs(v+vq,x,0) == 0, ...
    subs(t+tq,x,0) == 0, ...
    subs(u+uq,x,L) == 0, ...
    subs(v+vq,x,L) == 0, ...
    subs(t+tq,x,L) == 0, ...
    [C1,C2,C3,C4,C5,C6]);
%% ========================================================================
% FUERZAS DE EMPOTRAMIENTO PERFECTO
% Se conserva la misma convención:
%       [-N1;
%         V1;
%        -M1;
%         N2;
%        -V2;
%         M2]
% ========================================================================
FNE(:) = [ ...
    -subs(Nq+N, ...
    {C1,C2,C3,C4,C5,C6,x}, ...
    {c1,c2,c3,c4,c5,c6,0});
     subs(Vq+V, ...
    {C1,C2,C3,C4,C5,C6,x}, ...
    {c1,c2,c3,c4,c5,c6,0});
    -subs(Mq+M, ...
    {C1,C2,C3,C4,C5,C6,x}, ...
    {c1,c2,c3,c4,c5,c6,0});
     subs(Nq+N, ...
    {C1,C2,C3,C4,C5,C6,x}, ...
    {c1,c2,c3,c4,c5,c6,L});
    -subs(Vq+V, ...
    {C1,C2,C3,C4,C5,C6,x}, ...
    {c1,c2,c3,c4,c5,c6,L});
     subs(Mq+M, ...
    {C1,C2,C3,C4,C5,C6,x}, ...
    {c1,c2,c3,c4,c5,c6,L})
    ];
FNE = double(FNE);

%% ========================================================================
% PARTE III
% FORMULACIÓN MATRICIAL GENERAL
% Se escribe:
%       d = B*C + dp
%       f = S*C + fp
% donde:
%       C  = constantes de integración
%       dp = desplazamientos particulares
%       fp = fuerzas particulares
% ========================================================================
Un = [C1;
      C2;
      C3;
      C4;
      C5;
      C6];
%% ========================================================================
% MATRIZ B
%
% Relación:
%
%       d = B*C
%
% para la solución homogénea.
% ========================================================================
Bd = equationsToMatrix([u;
                        v;
                        t],Un);
%% Nodo inicial
B0 = double(subs(Bd,x,0));
%% Nodo final
BL = double(subs(Bd,x,L));
B  = [B0;
      BL];
%% ========================================================================
% MATRIZ S
% Estado interno:
%       s = [N;
%            V;
%            M]
% ========================================================================
Sf = equationsToMatrix([N;
                        V;
                        M],Un);
%% Nodo inicial
S0 = double(subs(Sf,x,0));
%% Nodo final
SL = double(subs(Sf,x,L));
%% ========================================================================
% MATRICES DE SIGNOS
% Nodo inicial:
%       [-N;
%         V;
%        -M]
% Nodo final:
%       [ N;
%        -V;
%         M]
% ========================================================================
Ri = [-1  0  0;
       0  1  0;
       0  0 -1];
Rj = [ 1  0  0;
       0 -1  0;
       0  0  1];
%% Matriz S completa

S = [Ri*S0;
     Rj*SL];
%% ========================================================================
% MATRIZ DE RIGIDEZ
% De:
%       d = B*C
%       C = inv(B)*d
% y:
%       f = S*C
% resulta:
%       f = S*inv(B)*d
% por tanto:
%       K = S/B
% ========================================================================
K = S/B;
%% Elimina pequeños errores numéricos de asimetría
K = 0.5*(K+K.');
%% ========================================================================
% COMPARACIÓN DE AMBAS MATRICES
% ========================================================================
error_K = max(max(abs(K-K_TE2)));
fprintf('\nError máximo entre K y K_TE2 = %.12e\n',error_K);
%% ========================================================================
% VECTOR PARTICULAR DE DESPLAZAMIENTOS
% La solución total puede escribirse:
%       d = B*C + dp
% con:
%       dp =
%       [uq(0);
%        vq(0);
%        tq(0);
%        uq(L);
%        vq(L);
%        tq(L)]
% ========================================================================
dp0 = [subs(uq,x,0);
       subs(vq,x,0);
       subs(tq,x,0)];
dpL = [subs(uq,x,L);
       subs(vq,x,L);
       subs(tq,x,L)];
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
% El estado particular es:
%       sp = [Nq;
%             Vq;
%             Mq]
% ========================================================================
sp0 = [subs(Nq,x,0);
       subs(Vq,x,0);
       subs(Mq,x,0)];
spL = [subs(Nq,x,L);
       subs(Vq,x,L);
       subs(Mq,x,L)];
%% Aplicar la convención de signos
fp = [Ri*sp0;
      Rj*spL];
fp = double(fp);
%% ========================================================================
% FUERZAS DE EMPOTRAMIENTO PERFECTO
% Tenemos:
%       d = B*C + dp
% por tanto:
%       C = B\(d-dp)
% y:
%       f = S*C + fp
% Sustituyendo:
%       f = S/B*d - S/B*dp + fp
%       f = K*d + fp - K*dp
% Si d = 0:
%       f0 = fp - K*dp
% ========================================================================
f0 = fp - K*dp;
%% ========================================================================
% FUERZAS NODALES EQUIVALENTES
% Si se adopta la ecuación del elemento:
%
%       fint = K*d - feq
% comparando con:
%       f = K*d + fp - K*dp
% se obtiene:
%       -feq = fp - K*dp
% por tanto:
%       feq = K*dp - fp
% y:
%       feq = -f0
% ========================================================================
feq = K*dp - fp;
%% ========================================================================
% COMPROBACIÓN:
% La solución directa FNE corresponde a las fuerzas de empotramiento:
%       FNE = f0
% ========================================================================
error_f0 = max(abs(FNE-f0));
fprintf('\nError máximo FNE directa vs f0 = %.12e\n',error_f0);
%% ========================================================================
% COMPROBACIÓN:
% Como:
%       feq = -f0
% debe cumplirse:
%       feq + FNE = 0
% =======================================================================
%% ========================================================================
% PARTE IV
% MATRIZ CINEMÁTICA BDEF
%
% A partir de las funciones de forma exactas:
%       u = Nu*d
%       v = Nv*d
%       t = Nt*d
% Las deformaciones son:
%       epsilon = du/dx
%       chi     = dt/dx
%       gamma   = t - dv/dx
% Por tanto:
%       e = Bdef*d
% ========================================================================
Baxial = diff(Nu,x);
Bflexion = diff(Nt,x);
Bcortante = Nt - diff(Nv,x);
Bdef = [Baxial;
        Bflexion;
        Bcortante];
Bdef = simplify(Bdef);
%% ========================================================================
% MATRIZ CONSTITUTIVA
%       D =
%       [EA   0    0
%         0   EI   0
%         0    0   Ks]
% ========================================================================
D = [EA 0  0;
     0  EI 0;
     0  0  Ks];
%% ========================================================================
% COMPROBACIÓN DE LA MATRIZ DE RIGIDEZ MEDIANTE:
%       K = integral(Bdef' * D * Bdef dx)
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
fprintf('\nError máximo K_B vs K = %.12e\n',error_KB);
%% ========================================================================
% PARTE V
% FUERZAS NODALES MEDIANTE INTEGRACIÓN DE FUNCIONES DE FORMA
%
% Para las cargas distribuidas:
%       feq_N = integral(Nu' * b dx)
%       feq_V = integral(Nv' * qv dx)
% El signo de qv debe ser consistente con la convención adoptada para v.
% En esta formulación:
%       V' = w
% mientras que las fuerzas nodales verticales adoptadas son:
%       [ V(0);
%        -V(L)]
% Por ello se verifica explícitamente el signo comparando con feq.
%
% ========================================================================
feq_N = int(Nu.'*b,x,0,L);
feq_V_pos = int(Nv.'*w,x,0,L);
feq_V_neg = int(Nv.'*(-w),x,0,L);
feq_N = double(feq_N);
feq_V_pos = double(feq_V_pos);
feq_V_neg = double(feq_V_neg);
%% Vector con carga transversal positiva
feq_int_pos = feq_N + feq_V_pos;
%% Vector con carga transversal negativa
feq_int_neg = feq_N + feq_V_neg;







