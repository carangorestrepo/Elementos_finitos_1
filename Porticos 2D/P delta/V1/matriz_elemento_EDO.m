function [Ke,feq,B,S,dp,fp] = ...
    matriz_elemento_EDO(Afun,Yp,x,L,D,R,Si,Sj)
%==========================================================================
% MATRIZ_ELEMENTO_EDO
%
% Obtiene la matriz de rigidez exacta y el vector de fuerzas nodales
% equivalentes a partir de la solucion general de un sistema de
% ecuaciones diferenciales lineales.
%
% La solucion del sistema debe expresarse como:
%
%       Y(x) = Afun(x)*C + Yp(x)
%
% donde:
%
%       Y     = vector de estado
%       Afun  = matriz fundamental de la solucion homogenea
%       C     = constantes de integracion
%       Yp    = solucion particular
%
% -------------------------------------------------------------------------
% ENTRADAS
%
% Afun : matriz fundamental simbolica
%
% Yp   : vector particular simbolico
%
% x    : variable simbolica
%
% L    : longitud del elemento
%
% D    : matriz que selecciona los desplazamientos del vector de estado
%
%       q = D*Y
%
% R    : matriz que selecciona las fuerzas internas
%
%       r = R*Y
%
% Si   : matriz de signos de fuerzas nodales en extremo inicial
%
% Sj   : matriz de signos de fuerzas nodales en extremo final
%
% -------------------------------------------------------------------------
% SALIDAS
%
% Ke   : matriz de rigidez del elemento
%
% feq  : vector de fuerzas nodales equivalentes
%
% B    : matriz desplazamiento-constantes
%
% S    : matriz fuerza-constantes
%
% dp   : desplazamientos nodales particulares
%
% fp   : fuerzas de extremo particulares
%
% Convencion:
%
%       fe = Ke*de - feq
%
%==========================================================================
%% ========================================================================
% 1. MATRIZ FUNDAMENTAL EN LOS EXTREMOS
% ========================================================================
A0 = subs(Afun,x,0);
AL = subs(Afun,x,L);
%% ========================================================================
% 2. SOLUCION PARTICULAR EN LOS EXTREMOS
% ========================================================================
P0 = subs(Yp,x,0);
PL = subs(Yp,x,L);
%% ========================================================================
% 3. RELACION DESPLAZAMIENTOS - CONSTANTES
%
% de = B*C + dp
% ========================================================================
B = [D*A0;
     D*AL];
%% Desplazamientos particulares
dp = [D*P0;
      D*PL];
%% ========================================================================
% 4. RELACION FUERZAS - CONSTANTES
%
% fe = S*C + fp
% ========================================================================
S = [Si*R*A0;
     Sj*R*AL];
%% Fuerzas particulares de extremo

fp = [Si*R*P0;
      Sj*R*PL];
%% ========================================================================
% 5. MATRIZ DE RIGIDEZ
%
% C = inv(B)*(de-dp)
%
% fe = S*inv(B)*de - S*inv(B)*dp + fp
%
% por tanto:
%
% Ke = S*inv(B)
%
% Numericamente:
%
% Ke = S/B
% ========================================================================
Ke = simplify(S/B);
%% ========================================================================
% 6. FUERZAS NODALES EQUIVALENTES
%
% fe = Ke*de - feq
%
% entonces:
%
% feq = Ke*dp - fp
% ========================================================================
feq = simplify(Ke*dp - fp);
end