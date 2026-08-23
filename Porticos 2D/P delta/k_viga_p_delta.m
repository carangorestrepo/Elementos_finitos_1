clc
clear
format long g

%% ========================================================================
% VIGA DE TIMOSHENKO CON EFECTO P-DELTA
%
% Este programa obtiene:
%
%   1. Solucion exacta homogenea del elemento.
%
%   2. Matriz de rigidez exacta mediante desplazamientos nodales
%      unitarios.
%
%   3. Funciones de forma exactas:
%
%          Nu(x) -> desplazamiento axial
%          Nv(x) -> desplazamiento transversal
%          Nt(x) -> giro
%
%   4. Fuerzas de empotramiento perfecto.
%
%   5. Fuerzas nodales equivalentes.
%
%   6. Formulacion matricial:
%
%          d = B*C + dp
%
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
%   7. Verificacion energetica:
%
%          K = integral(Bdef'*D*Bdef dx)
%
%            - integral(P*(Nv')'*(Nv') dx)
%
%   8. Verificacion de las fuerzas nodales mediante:
%
%          feq = integral(N'*q dx)
%
% ------------------------------------------------------------------------
% GRADOS DE LIBERTAD LOCALES
%
%       d = [u1;
%            v1;
%            t1;
%            u2;
%            v2;
%            t2]
%
% ------------------------------------------------------------------------
% ESTADO INTERNO
%
%       s = [N;
%            V;
%            M]
%
% donde:
%
%       N = fuerza axial incremental
%
%       V = fuerza transversal total
%
%       M = momento flector
%
% ------------------------------------------------------------------------
% CONVENCION DE FUERZAS NODALES
%
%       f = [-N(0);
%             V(0);
%            -M(0);
%             N(L);
%            -V(L);
%             M(L)]
%
% ------------------------------------------------------------------------
% ECUACIONES DE EQUILIBRIO
%
% Axial:
%
%       N' = -b(x)
%
% Transversal con P-Delta:
%
%       V' = -w(x)
%
%       M' = V - P*v'
%
% ------------------------------------------------------------------------
% ECUACIONES CONSTITUTIVAS / CINEMATICAS
%
%       u' = N/EA
%
%       t' = M/EI
%
% Para Timoshenko:
%
%       Q = Ks*(t-v')
%
% donde:
%
%       Q = cortante constitutivo
%
% Con efecto P-Delta:
%
%       Q = M'
%
%       V = Q + P*v'
%
% por tanto:
%
%       Q = V - P*v'
%
% Sustituyendo en:
%
%       Q = Ks*(t-v')
%
% resulta:
%
%       V-P*v' = Ks*(t-v')
%
% y:
%
%       t = (1-P/Ks)*v' + V/Ks
%
% Definiendo:
%
%       phiP = 1-P/Ks
%
% queda:
%
%       t = phiP*v' + V/Ks
%
% ------------------------------------------------------------------------
% ECUACION DIFERENCIAL TRANSVERSAL
%
% Eliminando M,V,t:
%
%       EI*phiP*v'''' + P*v''
%
%                             EI
%                       = -w + --*w''
%                             Ks
%
% Para carga trapezoidal:
%
%       w'' = 0
%
% por lo que:
%
%       EI*phiP*v'''' + P*v'' = -w
%
% Para la solucion homogenea:
%
%       EI*phiP*v'''' + P*v'' = 0
%
% ========================================================================
%% ========================================================================
% VARIABLES SIMBOLICAS
% ========================================================================
syms x
syms C1 C2 C3 C4 C5 C6
%% ========================================================================
% PROPIEDADES DEL MATERIAL Y DE LA SECCION
% ========================================================================
g = 9.8066502;             % [m/s^2] aceleracion de la gravedad
E = 24870062.3;            % [kPa] modulo de elasticidad
G = 0.4*E;                 % [kPa] modulo de cortante
L = 4;                     % [m] longitud del elemento
Ae = 0.4^2;                % [m^2] area
I = 0.4^4/12;              % [m^4] momento de inercia
%% ========================================================================
% FUERZA AXIAL INICIAL
%
% P es la fuerza axial inicial que genera el efecto P-Delta.
%
% Convencion:
%
%       P > 0  -> compresion
%
% En este desarrollo P se considera CONSTANTE a lo largo del elemento,
% correspondiente por ejemplo a una carga axial nodal.
%
% Las cargas axiales distribuidas b(x) que aparecen posteriormente se
% consideran cargas externas incrementales y NO modifican el estado
% inicial P utilizado para construir la rigidez geometrica.
% ========================================================================
P = 1000;
%% ========================================================================
% RIGIDECES
% ========================================================================
EA = E*Ae;                 % [kN]
EI = E*I;                  % [kN*m^2]
kappa = 5/6;               % coeficiente de correccion por cortante
Ks = kappa*G*Ae;           % [kN] rigidez efectiva a cortante
%% ========================================================================
% PARAMETRO P-DELTA
% De:
%       V-P*v' = Ks*(t-v')
% se obtiene:
%       t = (1-P/Ks)*v' + V/Ks
% Definiendo:
%       phiP = 1-P/Ks
% La ecuacion homogenea es:
%       EI*phiP*v'''' + P*v'' = 0
% o:
%       v'''' + la^2*v'' = 0
% donde:
%                    P
%       la^2 = ----------------
%                EI*phiP
% ========================================================================
phiP = 1 - P/Ks;
la = sqrt(P/(EI*phiP));

%% ========================================================================
% COMPROBACION DEL PARAMETRO phiP
% ========================================================================
if phiP <= 0
    warning('Se tiene P >= Ks. Revisar los datos y la formulacion.')
end
%% ========================================================================
% FORMA DE LA SOLUCION HOMOGENEA
% La ecuacion caracteristica:
%       r^4 + la^2*r^2 = 0
% puede escribirse:
%       r^2*(r^2+la^2) = 0
% Sus raices son:
%       r1 = 0
%       r2 = 0
%       r3 = +i*la
%       r4 = -i*la
% por tanto:
%       v =
%       C1 + C2*x
%       + C3*cos(la*x)
%       + C4*sin(la*x)
% ========================================================================
%% ========================================================================
% PARTE I
% SOLUCION HOMOGENEA
% ========================================================================
% Para obtener la matriz de rigidez:
%       b(x) = 0
%       w(x) = 0
% ========================================================================
%% ========================================================================
% DESPLAZAMIENTO TRANSVERSAL HOMOGENEO
% Solucion de:
%       EI*phiP*v'''' + P*v'' = 0
% ========================================================================
v = C1 + C2*x ...
  + C3*cos(la*x) ...
  + C4*sin(la*x);
%% ========================================================================
% MOMENTO HOMOGENEO
% Como:
%       t = phiP*v' + V/Ks
% y para la solucion homogenea:
%       V' = 0
% entonces:
%       t' = phiP*v''
% Como:
%       M = EI*t'
% resulta:
%       M = EI*phiP*v''
% ========================================================================
M = EI*phiP*diff(v,x,2);
%% ========================================================================
% CORTANTE HOMOGENEO
% Equilibrio con P-Delta:
%       M' = V-P*v'
% despejando:
%       V = M' + P*v'
% ========================================================================
V = diff(M,x) + P*diff(v,x);
%% ========================================================================
% CORTANTE CONSTITUTIVO HOMOGENEO
% Se puede definir:
%       Q = M'
% y por equilibrio:
%       Q = V-P*v'
% ========================================================================
Q = diff(M,x);
Q = simplify(Q);
%% ========================================================================
% GIRO HOMOGENEO
% De:
%       V-P*v' = Ks*(t-v')
% resulta
%       t = (1-P/Ks)*v' + V/Ks
% es decir:
%       t = phiP*v' + V/Ks
% Equivalentemente:
%       t = v' + Q/Ks
% ========================================================================
t = phiP*diff(v,x) + V/Ks;
%% ========================================================================
% PARTE AXIAL HOMOGENEA
% Sin carga axial distribuida:
%       N' = 0
% entonces:
%       N = C5
% ========================================================================
N = C5;
%% Desplazamiento axial
%       u' = N/EA
u = int(N/EA,x) + C6;
%% ========================================================================
% SIMPLIFICAR SOLUCION HOMOGENEA
% ========================================================================
N = simplify(N);
V = simplify(V);
M = simplify(M);
u = simplify(u);
v = simplify(v);
t = simplify(t);
%% ========================================================================
% COMPROBACION DE LA SOLUCION HOMOGENEA
% Debe cumplirse:
%       V' = 0
%       M' - V + P*v' = 0
%       t' - M/EI = 0
%       V-P*v' - Ks*(t-v') = 0
% ========================================================================
res_h1 = simplify(diff(V,x));
res_h2 = simplify( ...
         diff(M,x)-V+P*diff(v,x));
res_h3 = simplify( ...
         diff(t,x)-M/EI);
res_h4 = simplify( ...
         V-P*diff(v,x)-Ks*(t-diff(v,x)));
%% ========================================================================
% MATRIZ DE RIGIDEZ MEDIANTE DESPLAZAMIENTOS UNITARIOS
% ========================================================================
K_TE2 = sym(zeros(6,6));
Nu = sym(zeros(1,6));
Nv = sym(zeros(1,6));
Nt = sym(zeros(1,6));
%% ========================================================================
% CADA COLUMNA SE OBTIENE IMPONIENDO UN GDL UNITARIO
% ========================================================================

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
    % FUNCIONES DE FORMA EXACTAS
    %       u(x) = Nu*d
    %       v(x) = Nv*d
    %       t(x) = Nt*d
    % ====================================================================
    Nu(i) = subs(u, ...
        {C1,C2,C3,C4,C5,C6}, ...
        {c1,c2,c3,c4,c5,c6});
    Nv(i) = subs(v, ...
        {C1,C2,C3,C4,C5,C6}, ...
        {c1,c2,c3,c4,c5,c6});
    Nt(i) = subs(t, ...
        {C1,C2,C3,C4,C5,C6}, ...
        {c1,c2,c3,c4,c5,c6});
    %% ====================================================================
    % FUERZAS NODALES
    % Convencion:
    %       [-N1;
    %         V1;
    %        -M1;
    %         N2;
    %        -V2;
    %         M2]
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
%% ========================================================================
% SIMPLIFICAR FUNCIONES DE FORMA
% ========================================================================
Nu = simplify(Nu);
Nv = simplify(Nv);
Nt = simplify(Nt);
%% ========================================================================
% MATRIZ DE RIGIDEZ NUMERICA
% ========================================================================
K_TE2 = double(K_TE2);

%% Eliminar pequenos errores de redondeo
K_TE2 = 0.5*(K_TE2 + K_TE2.');
%% ========================================================================
% PARTE II
% CARGAS DISTRIBUIDAS
% ========================================================================
ba = 25;                   % [kN/m] carga axial inicial
bb = 35;                   % [kN/m] carga axial final
wa = 25;                   % [kN/m] carga transversal inicial
wb = 30;                   % [kN/m] carga transversal final
%% ========================================================================
% CARGA AXIAL TRAPEZOIDAL
% Positiva en la direccion +u
% ========================================================================
b = ba + (bb-ba)*x/L;
%% ========================================================================
% CARGA TRANSVERSAL TRAPEZOIDAL
% Positiva en la direccion +v
% =======================================================================
w = wa + (wb-wa)*x/L;
%% ========================================================================
% SOLUCION PARTICULAR TRANSVERSAL
% EDO:
%       EI*phiP*v'''' + P*v''
%
%                             EI
%                       = -w + --*w''
%                             Ks
% Para carga trapezoidal:
%
%       w = wa + (wb-wa)*x/L
% por tanto:
%       w'' = 
% y:
%       EI*phiP*v'''' + P*v'' = -w
%
% Para una solucion particular cubica:
%       vq'''' = 0
% entonces:
%       P*vq'' = -w
% Una solucion particular conveniente es:
%       vq =
%       -wa*x^2/(2*P)
%       -(wb-wa)*x^3/(6*P*L)
% ========================================================================
vq = -wa*x^2/(2*P) - (wb-wa)/L*x^3/(6*P);
%% ========================================================================
% MOMENTO PARTICULAR
% Como:
%       t = phiP*v' + V/Ks
% entonces:
%       t' = phiP*v'' + V'/Ks
% y:
%       V' = -w
% resulta:
%       t' = phiP*v'' - w/Ks
% Como:
%       M = EI*t'
% se obtiene:
%       Mq = EI*(phiP*vq'' - w/Ks)
% ========================================================================
Mq = EI*(phiP*diff(vq,x,2) - w/Ks);
Mq = simplify(Mq);
%% ========================================================================
% CORTANTE PARTICULAR
% Equilibrio:
%       Mq' = Vq-P*vq'
% despejando:
%       Vq = Mq' + P*vq'
% ========================================================================
Vq = diff(Mq,x) + P*diff(vq,x);
Vq = simplify(Vq);
%% ========================================================================
% CORTANTE CONSTITUTIVO PARTICULAR
%       Qq = Mq'
% y:
%       Qq = Vq-P*vq'
% ========================================================================
Qq = diff(Mq,x);
Qq = simplify(Qq);
%% ========================================================================
% GIRO PARTICULAR
% De:
%       Vq-P*vq' = Ks*(tq-vq')
% resulta:
%       tq = phiP*vq' + Vq/Ks
% Equivalentemente:
%       tq = vq' + Qq/Ks
% ========================================================================
tq = phiP*diff(vq,x) + Vq/Ks;
tq = simplify(tq);
%% Simplificar
vq = simplify(vq);
%% ========================================================================
% COMPROBACION TRANSVERSAL
%
% Debe cumplirse:
%
%       Vq' + w = 0
%
%       Mq' - Vq + P*vq' = 0
%
%       tq' - Mq/EI = 0
%
%       Vq-P*vq' - Ks*(tq-vq') = 0
%
% ========================================================================

res_trans_1 = simplify( ...
              diff(Vq,x)+w);

res_trans_2 = simplify( ...
              diff(Mq,x)-Vq+P*diff(vq,x));

res_trans_3 = simplify( ...
              diff(tq,x)-Mq/EI);

res_trans_4 = simplify( ...
              Vq-P*diff(vq,x) ...
              -Ks*(tq-diff(vq,x)));
disp([res_trans_1;
      res_trans_2;
      res_trans_3;
      res_trans_4]);
%% ========================================================================
% COMPROBACION DE LAS DOS EXPRESIONES DEL GIRO
%       tq = phiP*vq' + Vq/Ks
%       tq = vq' + Qq/Ks
% ========================================================================
res_tq = simplify( ...
         tq-(diff(vq,x)+Qq/Ks));
disp(' ')
disp('Comprobacion alternativa del giro tq')
disp(res_tq);
%% ========================================================================
% COMPROBACION DE LA EDO DE CUARTO ORDEN
%       EI*phiP*vq'''' + P*vq''
%                             EI
%                       = -w + --*w''
%                             Ks
% ========================================================================
res_EDO = simplify( ...
          EI*phiP*diff(vq,x,4) ...
          +P*diff(vq,x,2) ...
          +w ...
          -EI/Ks*diff(w,x,2));
disp(' ')
disp('Comprobacion EDO transversal')
disp(res_EDO);
%% ========================================================================
% SOLUCION PARTICULAR AXIAL
% Para una carga b(x) positiva en +u:
%       N' + b = 0
% por tanto:
%       Nq' = -b
% =======================================================================
Nq = -int(b,x);
%% Desplazamiento axial particular
%       uq' = Nq/EA
uq = int(Nq/EA,x);
%% Simplificar
Nq = simplify(Nq);
uq = simplify(uq);
%% ========================================================================
% COMPROBACION AXIAL
%       Nq' + b = 0
% =======================================================================
res_axial = simplify(diff(Nq,x)+b);
%% ========================================================================
% PARTE III
% FUERZAS DE EMPOTRAMIENTO MEDIANTE SOLUCION DIRECTA
% Solucion total:
%       uT = u + uq
%       vT = v + vq
%       tT = t + tq
% Para empotramiento perfecto:
%       d = 0
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
% FUERZAS DE EMPOTRAMIENTO PERFECTO DIRECTAS
%       f0 =
%       [-N1;
%         V1;
%        -M1;
%         N2;
%        -V2;
%         M2]
% ========================================================================
F0_directa = [ ...

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
F0_directa = double(F0_directa);
%% ========================================================================
% PARTE IV
% FORMULACION MATRICIAL GENERAL
%       d = B*C + dp
%
%       f = S*C + fp
% ========================================================================
Un = [C1;
      C2;
      C3;
      C4;
      C5;
      C6];
%% ========================================================================
% MATRIZ B
% Relacion entre las constantes C y los desplazamientos nodales.
%       d_h = B*C
% ========================================================================

Bd = equationsToMatrix( ...
     [u;
      v;
      t], ...
      Un);
%% En x = 0
B0 = double(subs(Bd,x,0));
%% En x = L
BL = double(subs(Bd,x,L));
%% Matriz completa

B = [B0;
     BL];

%% ========================================================================
% MATRIZ S
% Estado homogeneo:
%       s = [N;
%            V;
%            M]
% ========================================================================

Sf = equationsToMatrix( ...
     [N;
      V;
      M], ...
      Un);
%% Estado en x = 0
S0 = double(subs(Sf,x,0));
%% Estado en x = L
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
%% Matriz S
S = [Ri*S0;
     Rj*SL];
%% ========================================================================
% MATRIZ DE RIGIDEZ
%       d = B*C
%       f = S*C
% entonces:
%       C = B\d
%       f = S/B*d
% por tanto:
%       K = S/B
% ========================================================================
K = S/B;
%% Eliminar errores pequenos de redondeo
K = 0.5*(K+K.');
%% ========================================================================
% VECTOR PARTICULAR DE DESPLAZAMIENTOS
%       d = B*C + dp
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
%% ========================================================================
% VECTOR PARTICULAR DE FUERZAS
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
%% Aplicar convencion nodal
fp = [Ri*sp0;
      Rj*spL];
fp = double(fp);
%% ========================================================================
% FUERZAS DE EMPOTRAMIENTO PERFECTO
% Partimos de:
%       d = B*C + dp
% entonces:
%       C = B\(d-dp)
% y:
%       f = S*C + fp
% sustituyendo:
%       f = K*d - K*dp + fp
% Para:
%       d = 0
% resulta:
%       f0 = fp - K*dp
% ========================================================================
f0 = fp - K*dp;
%% ========================================================================
% FUERZAS NODALES EQUIVALENTES
% Si:
%       fint = K*d - feq
% y:
%       f = K*d + fp - K*dp
% entonces:
%       feq = K*dp - fp
%       feq = -f0
% ========================================================================
feq = K*dp - fp;

%% ========================================================================
% PARTE V
% MATRIZ CINEMATICA Bdef
%
% A partir de:
%
%       u = Nu*d
%
%       v = Nv*d
%
%       t = Nt*d
%
% tenemos:
%
% Deformacion axial:
%
%       epsilon = u'
%
% Curvatura:
%
%       chi = t'
%
% Deformacion por cortante:
%
%       gamma = t - v'
%
% ========================================================================
Baxial = diff(Nu,x);
Bflexion = diff(Nt,x);
Bcortante = Nt - diff(Nv,x);
Bdef = [Baxial;
        Bflexion;
        Bcortante];
Bdef = simplify(Bdef);
%% ========================================================================
% MATRIZ CONSTITUTIVA DE LA VIGA
% ========================================================================
D = [EA 0  0;
     0  EI 0;
     0  0  Ks];
%% ========================================================================
% PARTE VI
% MATRIZ DE RIGIDEZ POR ENERGIA
%
% Para Timoshenko + P-Delta:
%
%       Pi =
%
%       1/2 integral[
%           EA*(u')^2
%         + EI*(t')^2
%         + Ks*(t-v')^2
%         - P*(v')^2
%       ] dx
% Por tanto:
%
%       Ktotal = Kviga - KP
%
% donde:
%
%       Kviga =
%
%          integral(Bdef'*D*Bdef dx)
%
% y:
%
%       KP =
%
%          integral(P*(Nv')'*(Nv') dx)
%
% ========================================================================
%% Rigidez de la viga
K_beam = int(Bdef.'*D*Bdef,x,0,L);
%% ========================================================================
% RIGIDEZ GEOMETRICA P-DELTA
%
% La energia geometrica depende de la pendiente v':
%
%       Ug = -1/2*integral(P*(v')^2 dx)
%
% Como:
%
%       v = Nv*d
%
% entonces:
%
%       v' = Nv'*d
%
% por tanto:
%
%       KP = integral(P*(Nv')^T*(Nv') dx)
%
% IMPORTANTE:
%
% NO:
%
%       P*Nv^T*Nv
%
% sino:
%
%       P*(Nv')^T*(Nv')
%
% ========================================================================

Bgeo = diff(Nv,x);
K_P = int(P*(Bgeo.'*Bgeo),x,0,L);
%% Rigidez total
K_B = K_beam - K_P;
%% Convertir a numerico
K_beam = double(K_beam);
K_P = double(K_P);
K_B = double(K_B);
%% Simetrizacion numerica
K_beam = 0.5*(K_beam+K_beam.');
K_P = 0.5*(K_P+K_P.');
K_B = 0.5*(K_B+K_B.');

%% ========================================================================
% PARTE VII
% FUERZAS NODALES CONSISTENTES MEDIANTE FUNCIONES DE FORMA

% Carga externa:
%       b(x) -> axial
%
%       w(x) -> transversal
%
% IMPORTANTE:
%
% La fuerza axial inicial P que produce el efecto P-Delta
% NO se introduce nuevamente en el vector de cargas externas.
%
% Su efecto ya esta incorporado en:
%
%       - las funciones de forma exactas
%
%       - la matriz de rigidez del elemento
%
%       - la rigidez geometrica
%
% Las cargas b(x) y w(x) son las cargas distribuidas externas
% utilizadas para calcular las fuerzas nodales equivalentes.
%
% ========================================================================
%% ========================================================================
% CARGA AXIAL CONSISTENTE
%
%       feq_axial = integral(Nu'*b dx)
%
% ========================================================================
feq_axial = int(Nu.'*b,x,0,L);


%% ========================================================================
% CARGA TRANSVERSAL CONSISTENTE
%
%       feq_trans = integral(Nv'*w dx)
%
% ========================================================================

feq_trans = int(Nv.'*w,x,0,L);


%% Convertir

feq_axial = double(feq_axial);

feq_trans = double(feq_trans);


%% Vector total

feq_int = feq_axial + feq_trans;

