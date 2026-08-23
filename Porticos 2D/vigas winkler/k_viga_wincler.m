clc
clear
format long g
%% ========================================================================
% VIGA DE TIMOSHENKO SOBRE FUNDACION ELASTICA DE WINKLER
% Este programa obtiene:
%
%   1. Solucion exacta homogénea del elemento.
%   2. Matriz de rigidez exacta mediante desplazamientos nodales
%      unitarios.
%   3. Funciones de forma exactas:
%          Nu(x) -> desplazamiento axial
%          Nv(x) -> desplazamiento transversal
%          Nt(x) -> giro
%   4. Fuerzas de empotramiento perfecto.
%
%   5. Fuerzas nodales equivalentes.
%
%   6. Formulacion matricial:
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
%   7. Verificacion energética:
%          K = integral(Bdef'*D*Bdef dx)
%            + integral(k*Nv'*Nv dx)
%   8. Verificacion de las fuerzas nodales mediante:
%          feq = integral(N'*q dx)
% ------------------------------------------------------------------------
% GRADOS DE LIBERTAD LOCALES
%       d = [u1;
%            v1;
%            t1;
%            u2;
%            v2;
%            t2]
% ------------------------------------------------------------------------
% ESTADO INTERNO
%       s = [N;
%            V;
%            M]
% ------------------------------------------------------------------------
% CONVENCION DE FUERZAS NODALES
%
%       f = [-N(0);
%             V(0);
%            -M(0);
%             N(L);
%            -V(L);
%             M(L)]
% ------------------------------------------------------------------------
% ECUACIONES DE EQUILIBRIO
% Axial:
%       N' = -b(x)
% Transversal con Winkler:
%       V' = -k*v + w(x)
%       M' = V
% ------------------------------------------------------------------------
% ECUACIONES CONSTITUTIVAS / CINEMATICAS
%       u' = N/EA
%       t' = M/EI
%       v' = t - V/Ks
% donde:
%       Ks = kappa*G*A
% ------------------------------------------------------------------------
% ECUACION DIFERENCIAL TRANSVERSAL
% Eliminando M,V,t:
%       v'''' - (k/Ks)*v'' + (k/EI)*v
%                 w          w''
%             = ----- - -----------
%                 EI         Ks
% Para carga lineal:
%       w'' = 0
% por lo que:
%       v'''' - (k/Ks)*v'' + (k/EI)*v = w/EI
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
% FUNDACION DE WINKLER
% ========================================================================
k = 500;                   % [kN/m^2]
                           %
                           % Reaccion por unidad de longitud:
                           %
                           %       r(x) = -k*v(x)
                           %
                           % Para un elemento 2D, k debe entenderse como
                           % rigidez de fundacion por unidad de longitud.
%% ========================================================================
% RIGIDECES
% ========================================================================
EA = E*Ae;                 % [kN
EI = E*I;                  % [kN*m^2]
kappa = 5/6;               % coeficiente de correccion por cortante
Ks = kappa*G*Ae;           % [kN] rigidez efectiva a cortante
%% ========================================================================
% PARAMETROS DE LA SOLUCION HOMOGENEA
% EDO homogénea:
%       v'''' - a*v'' + b*v = 0
% con:
%       a = k/Ks
%       b = k/EI
% ========================================================================
a = k/Ks;
bb0 = k/EI;
%% Parámetros alfa y beta
alfa = 0.5*sqrt(a + 2*sqrt(bb0));
beta = 0.5*sqrt(2*sqrt(bb0) - a);
%% ========================================================================
% ADVERTENCIA SOBRE beta
% La forma:
% cosh(alfa*x)*cos(beta*x)
% corresponde al caso:
%       2*sqrt(k/EI) > k/Ks
% que es el caso de los datos utilizados.
% ========================================================================
if (2*sqrt(bb0)-a) <= 0
    warning(['Para estos datos beta no es real. ',...
             'Debe utilizarse la forma general basada en las raíces ',...
             'del polinomio característico.'])
end
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
%       v'''' - (k/Ks)*v'' + (k/EI)*v = 0
% ========================================================================
v = C1*cosh(alfa*x)*cos(beta*x) + ...
    C2*cosh(alfa*x)*sin(beta*x) + ...
    C3*sinh(alfa*x)*cos(beta*x) + ...
    C4*sinh(alfa*x)*sin(beta*x);
%% ========================================================================
% MOMENTO HOMOGENEO
% En general:
%       M = EI*(v'' - k/Ks*v + w/Ks)
% Para la solucion homogénea:
%       w = 0
% por tanto:
%       M = EI*(v'' - k/Ks*v)
% ========================================================================
M = EI*(diff(v,x,2) - k*v/Ks);
%% ========================================================================
% CORTANTE HOMOGENEO
%       M' = V
% ========================================================================
V = diff(M,x);
%% ========================================================================
% GIRO HOMOGENEO
% De:
%       v' = t - V/Ks
% resulta:
%       t = v' + V/Ks
% ========================================================================
t = diff(v,x) + V/Ks;
%% ========================================================================
% PARTE AXIAL HOMOGENEA
% Sin carga axial:
%       N' = 0
% ========================================================================
N = C5;
%% Desplazamiento axial
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
% COMPROBACION DE LA EDO HOMOGENEA DE WINKLER
% Debe cumplirse:
%       V' + k*v = 0
% ========================================================================
res_h = simplify(diff(V,x) + k*v);
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
    %
    %       u(x) = Nu*d
    %
    %       v(x) = Nv*d
    %
    %       t(x) = Nt*d
    %
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
    %
    % Convencion:
    %
    %       [-N1;
    %         V1;
    %        -M1;
    %         N2;
    %        -V2;
    %         M2]
    %
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
%% Eliminar pequeños errores de redondeo
K_TE2 = 0.5*(K_TE2 + K_TE2.');
%% ========================================================================
% COMPROBACION DE SIMETRIA
% ========================================================================
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
%
% Positiva en la direccion +u
% ========================================================================
b = ba + (bb-ba)*x/L;
%% ========================================================================
% CARGA TRANSVERSAL TRAPEZOIDAL
% Positiva en la direccion +v
% ========================================================================
w = wa + (wb-wa)*x/L;
%% ========================================================================
% SOLUCION PARTICULAR TRANSVERSAL
% EDO:
%       v'''' - k/Ks*v'' + k/EI*v
%                         w      w''
%                    = ------- - ----
%                         EI     Ks
% Para carga lineal:
%       w'' = 0
% Una solucion particular conveniente es:
%       vq = w/k
% porque:
%       k*vq = w
% ========================================================================
vq = w/k;
%% ========================================================================
% MOMENTO PARTICULAR
% IMPORTANTE:
% Se debe utilizar w y NO la variable qy=0 de la solucion homogénea.
%       Mq = EI*(vq'' - k/Ks*vq + w/Ks)
%
% ========================================================================
Mq = EI*(diff(vq,x,2) - k*vq/Ks + w/Ks);
Mq = simplify(Mq);
%% ========================================================================
% CORTANTE PARTICULAR
%       Vq = Mq'
% ========================================================================
Vq = diff(Mq,x);
Vq = simplify(Vq);
%% ========================================================================
% GIRO PARTICULAR
%       tq = vq' + Vq/Ks
% ========================================================================
tq = diff(vq,x) + Vq/Ks;
tq = simplify(tq);
%% Simplificar
vq = simplify(vq);
%% ========================================================================
% COMPROBACION DE LA EDO PARTICULAR
% Debe cumplirse:
%       Vq' + k*vq - w = 0
% ========================================================================
res_p = simplify(diff(Vq,x) + k*vq - w);
%% ========================================================================
% SOLUCION PARTICULAR AXIAL
% Para una carga b(x) positiva en +u:
%       N' + b = 0
% por tanto:
%       Nq' = -b
% =======================================================================
Nq = -int(b,x);
%% Desplazamiento axial particular
uq = int(Nq/EA,x);
%% Simplificar
Nq = simplify(Nq);
uq = simplify(uq);
%% ========================================================================
% COMPROBACION AXIAL
%       Nq' + b = 0
% ========================================================================
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
%
%       d = B*C + dp
%
%       f = S*C + fp
%
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
% Relacion entre las constantes C y los desplazamientos nodales.
%
%       d_h = B*C
%
% ========================================================================
Bd = equationsToMatrix([u;
                        v;
                        t],Un);
%% En x = 0
B0 = double(subs(Bd,x,0));
%% En x = L
BL = double(subs(Bd,x,L));
%% Matriz completa
B = [B0;
     BL];
%% ========================================================================
% MATRIZ S
%
% Estado homogéneo:
%
%       s = [N;
%            V;
%            M]
%
% ========================================================================
Sf = equationsToMatrix([N;
                        V;
                        M],Un);
%% Estado en x = 0
S0 = double(subs(Sf,x,0));
%% Estado en x = L
SL = double(subs(Sf,x,L));
%% ========================================================================
% MATRICES DE SIGNOS
% Nodo inicial:
%
%       [-N;
%         V;
%        -M]
% Nodo final:
%
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
%% Eliminar errores pequeños de redondeo
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
%
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
%% Aplicar convención nodal
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
% A partir de:
%       u = Nu*d
%       v = Nv*d
%       t = Nt*d
% tenemos:
% Deformacion axial:
%       epsilon = u'
% Curvatura:
%       chi = t'
% Deformacion por cortante:
%       gamma = t - v'
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
% Para Timoshenko + Winkler:
%       U =
%       1/2 integral[
%           EA*(u')^2
%         + EI*(t')^2
%         + Ks*(t-v')^2
%         + k*v^2 ] dx
% Por tanto:
%       Ktotal = Kviga + Kwinkler
% donde:
%       Kviga =
%          integral(Bdef'*D*Bdef dx)
% y:
%       Kwinkler =
%          integral(k*Nv'*Nv dx)
% ========================================================================
%% Rigidez de la viga
K_beam = int(Bdef.'*D*Bdef,x,0,L);
%% Rigidez de la fundación de Winkler
K_winkler = int(k*(Nv.'*Nv),x,0,L);
%% Rigidez total
K_B = K_beam + K_winkler;
%% Convertir a numérico
K_beam = double(K_beam);
K_winkler = double(K_winkler);
K_B = double(K_B);
%% Simetrización numérica
K_beam = 0.5*(K_beam+K_beam.');
K_winkler = 0.5*(K_winkler+K_winkler.');
K_B = 0.5*(K_B+K_B.');
%% ========================================================================
% COMPROBACION
% Debe cumplirse:
%       K_B = K
% ========================================================================
%% ========================================================================
% PARTE VII
% FUERZAS NODALES CONSISTENTES MEDIANTE FUNCIONES DE FORMA
% Carga externa:
%       b(x) -> axial
%       w(x) -> transversal
% IMPORTANTE:
% La reaccion de Winkler:
%       -k*v
% NO se introduce en el vector de cargas externas.
% Ya esta incluida dentro de:
%       K_winkler
% ========================================================================
%% ========================================================================
% CARGA AXIAL CONSISTENTE
%       feq_axial = integral(Nu'*b dx)
% ========================================================================
feq_axial = int(Nu.'*b,x,0,L);
%% ========================================================================
% CARGA TRANSVERSAL CONSISTENTE
%       feq_trans = integral(Nv'*w dx)
% ========================================================================
feq_trans = int(Nv.'*w,x,0,L);
%% Convertir
feq_axial = double(feq_axial);
feq_trans = double(feq_trans);
%% Vector total
feq_int = feq_axial + feq_trans;
%% ========================================================================
% COMPARACION
% En una formulacion variacional consistente debe cumplirse:
%       feq_int = feq
% =======================================================================
%% ========================================================================
% PARTE VIII
% COMPROBACIONES DE EQUILIBRIO GLOBAL
% ========================================================================
%% Carga axial total
Rb = double(int(b,x,0,L));
%% Carga transversal externa total
Rw = double(int(w,x,0,L));
%% Momento de la carga transversal respecto al nodo inicia
Mw = double(int(w*x,x,0,L));
%% ========================================================================
% PARTE IX
% COMPROBACION DE LAS FUNCIONES DE FORMA
%
% En x = 0 y x = L deben reproducir exactamente los GDL.
% ========================================================================
N0 = [subs(Nu,x,0);
      subs(Nv,x,0);
      subs(Nt,x,0)];
NL = [subs(Nu,x,L);
      subs(Nv,x,L);
      subs(Nt,x,L)];
disp(double(NL))
%% ========================================================================
% PARTE X
% COMPROBACION DIFERENCIAL DE LAS FUNCIONES DE FORMA
% Para cada columna Nv_j debe cumplirse:
%       Vj' + k*Nv_j = 0
% porque las funciones de forma corresponden al problema homogéneo.
% ========================================================================
res_N = sym(zeros(1,6));
for j = 1:6
    %% Momento correspondiente a la funcion de forma
    Mj = EI*(diff(Nv(j),x,2) - k*Nv(j)/Ks);
    %% Cortante
    Vj = diff(Mj,x);
    %% Residuo
    res_N(j) = simplify(diff(Vj,x) + k*Nv(j));
end