%% ========================================================================
% VIGA DE TIMOSHENKO SOBRE FUNDACION ELASTICA TIPO WINKLER
%
% Sistema de ecuaciones diferenciales de primer orden:
%
% V' = -k*v + q        (1) Equilibrio de fuerzas cortantes
% M' =  V              (2) Relacion momento-cortante
% t' =  M/EI           (3) Relacion giro-momento
% v' =  t - V/Ac       (4) Relacion desplazamiento-giro
%
% donde:
%
%   v  = desplazamiento transversal
%   t  = giro de la seccion
%   M  = momento flector
%   V  = fuerza cortante
%   EI = rigidez a flexion
%   Ac = rigidez efectiva a cortante = kappa*G*A
%   k  = coeficiente de balasto
%   q  = carga transversal distribuida
%
% Objetivo:
%
% Reducir el sistema anterior a una unica ecuacion diferencial de
% cuarto orden en funcion de v(x), obtener su solucion real y deducir
% explicitamente los parametros alfa y beta.
%
%% ========================================================================


%% ========================================================================
% 1. ECUACION DIFERENCIAL EN TERMINOS DEL DESPLAZAMIENTO v
%% ========================================================================
%
% Partimos de la ecuacion cinematica:
%
%       v' = t - V/Ac                                      (4)
%
% Derivando respecto a x:
%
%       v'' = t' - V'/Ac
%
% Utilizando:
%
%       V' = -k*v + q                                      (1)
%
%       t' = M/EI                                          (3)
%
% resulta:
%
%       v'' = M/EI - (-k*v + q)/Ac
%
% por tanto:
%
%       v'' = M/EI + k*v/Ac - q/Ac
%
% Despejando M:
%
%       M/EI = v'' - k*v/Ac + q/Ac
%
% finalmente:
%
%       M = EI*(v'' - k*v/Ac + q/Ac)
%


%% ========================================================================
% 2. DERIVACION DEL MOMENTO Y ECUACION DE CUARTO ORDEN
%% ========================================================================
%
% Como:
%
%       M' = V                                             (2)
%
% derivando la expresion anterior:
%
%       M' = EI*(v''' - k*v'/Ac + q'/Ac)
%
% Por tanto:
%
%       V = EI*(v''' - k*v'/Ac + q'/Ac)
%
% Derivando nuevamente:
%
%       M'' = EI*(v'''' - k*v''/Ac + q''/Ac)
%
% Pero:
%
%       M' = V
%
% entonces:
%
%       M'' = V'
%
% y de la ecuacion de equilibrio:
%
%       V' = -k*v + q                                     (1)
%
% Por tanto:
%
% EI*(v'''' - k*v''/Ac + q''/Ac) = -k*v + q
%
% Desarrollando:
%
% EI*v'''' - EI*k/Ac*v'' + EI/Ac*q'' = -k*v + q
%
% Reorganizando:
%
%       EI*v'''' - EI*k/Ac*v'' + k*v
%                       = q - EI/Ac*q''
%
% Dividiendo entre EI:
%
%       v'''' - k/Ac*v'' + k/EI*v
%                       = q/EI - q''/Ac
%


%% ========================================================================
% 3. CARGA TRAPEZOIDAL
%% ========================================================================
%
% Una carga trapezoidal se puede escribir como:
%
%                       q2-q1
%       q(x) = q1 + -----------*x
%                         L
%
% Por tanto:
%
%       q'  = (q2-q1)/L
%
%       q'' = 0
%
% Debido a que q''=0, la ecuacion diferencial se simplifica a:
%
%       v'''' - k/Ac*v'' + k/EI*v = q(x)/EI
%
% o equivalentemente:
%
%       EI*v'''' - EI*k/Ac*v'' + k*v = q(x)
%


%% ========================================================================
% 4. ECUACION DIFERENCIAL HOMOGENEA
%% ========================================================================
%
% Para obtener la solucion complementaria hacemos q(x)=0:
%
%       v'''' - k/Ac*v'' + k/EI*v = 0
%
% Definimos:
%
%       a = k/Ac
%
%       b = k/EI
%
% Entonces:
%
%       v'''' - a*v'' + b*v = 0
%

a = k/Ac;
b = k/EI;


%% ========================================================================
% 5. ECUACION CARACTERISTICA
%% ========================================================================
%
% Se propone una solucion exponencial:
%
%       v = exp(r*x)
%
% cuyas derivadas son:
%
%       v''   = r^2*exp(r*x)
%
%       v'''' = r^4*exp(r*x)
%
% Sustituyendo:
%
%       r^4*exp(r*x) - a*r^2*exp(r*x)
%                    + b*exp(r*x) = 0
%
% Como exp(r*x) ~= 0:
%
%       r^4 - a*r^2 + b = 0
%
% Hacemos:
%
%       z = r^2
%
% obteniendo:
%
%       z^2 - a*z + b = 0
%
% Aplicando la formula cuadratica:
%
%                 a +/- sqrt(a^2 - 4*b)
%       z1,z2 = --------------------------
%                            2
%
% Como z = r^2:
%
%              a +/- sqrt(a^2 - 4*b)
%       r^2 = ----------------------------
%                         2
%


%% ========================================================================
% 6. CASO DE RAICES COMPLEJAS
%% ========================================================================
%
% Estamos interesados en el caso:
%
%       a^2 - 4*b < 0
%
% equivalente a:
%
%       a^2 < 4*b
%
% Entonces:
%
%       sqrt(a^2 - 4*b) = i*sqrt(4*b - a^2)
%
% y:
%
%              a       sqrt(4*b-a^2)
%       r^2 = --- +/- i---------------
%              2              2
%
% Consideremos una de las raices:
%
%              a       sqrt(4*b-a^2)
%       r^2 = --- + i ----------------
%              2              2
%
% Ahora escribimos la raiz compleja r en la forma:
%
%       r = alfa + i*beta
%
% El objetivo es encontrar alfa y beta.
%


%% ========================================================================
% 7. DEDUCCION DE alfa Y beta
%% ========================================================================
%
% Si:
%
%       r = alfa + i*beta
%
% entonces:
%
%       r^2 = (alfa + i*beta)^2
%
% Desarrollando:
%
%       r^2 = alfa^2 + 2*i*alfa*beta + (i*beta)^2
%
% como:
%
%       i^2 = -1
%
% resulta:
%
%       r^2 = alfa^2 - beta^2 + 2*i*alfa*beta
%
% Pero anteriormente obtuvimos:
%
%              a       sqrt(4*b-a^2)
%       r^2 = --- + i ----------------
%              2              2
%
% Igualando las partes reales:
%
%                       a
%       alfa^2-beta^2 = -
%                       2                       (A)
%
% Igualando las partes imaginarias:
%
%                           sqrt(4*b-a^2)
%       2*alfa*beta = ------------------------
%                                  2           (B)
%
% Necesitamos ahora encontrar otra relacion entre alfa y beta.
%


%% ========================================================================
% 8. MODULO DE r^2
%% ========================================================================
%
% Tenemos:
%
%              a       sqrt(4*b-a^2)
%       r^2 = --- + i ----------------
%              2              2
%
% El modulo es:
%
%                 _______________________________
%                / (a/2)^2 + (sqrt(4*b-a^2)/2)^2
%       |r^2| = V
%
% Entonces:
%
%                 __________________________
%                / a^2     4*b-a^2
%       |r^2| = V  ---- + -----------
%                   4          4
%
%                 __________
%                / 4*b
%       |r^2| = V  ----
%                   4
%
% Por tanto:
%
%       |r^2| = sqrt(b)
%
% Por otra parte, como:
%
%       r = alfa + i*beta
%
% tenemos:
%
%       |r|^2 = alfa^2 + beta^2
%
% y:
%
%       |r^2| = |r|^2
%
% Por consiguiente:
%
%       alfa^2 + beta^2 = sqrt(b)                 (C)
%


%% ========================================================================
% 9. CALCULO DE alfa
%% ========================================================================
%
% Tenemos las ecuaciones:
%
%       alfa^2 - beta^2 = a/2                     (A)
%
%       alfa^2 + beta^2 = sqrt(b)                 (C)
%
% Sumando (A)+(C):
%
%       2*alfa^2 = a/2 + sqrt(b)
%
% entonces:
%
%                       a
%       alfa^2 = sqrt(b)/2 + ---
%                       4
%
% Reorganizando:
%
%                 1
%       alfa^2 = ---*(a + 2*sqrt(b))
%                 4
%
% Finalmente:
%
%                   1
%       alfa = -----------*sqrt(a + 2*sqrt(b))
%                   2
%

alfa = 1/2*sqrt(a + 2*sqrt(b));


%% ========================================================================
% 10. CALCULO DE beta
%% ========================================================================
%
% Nuevamente:
%
%       alfa^2 + beta^2 = sqrt(b)                 (C)
%
%       alfa^2 - beta^2 = a/2                     (A)
%
% Restando (C)-(A):
%
%       2*beta^2 = sqrt(b) - a/2
%
% entonces:
%
%                        a
%       beta^2 = sqrt(b)/2 - ---
%                        4
%
% Reorganizando:
%
%                1
%       beta^2 = ---*(2*sqrt(b) - a)
%                4
%
% Finalmente:
%
%                  1
%       beta = -----------*sqrt(2*sqrt(b) - a)
%                  2
%

beta = 1/2*sqrt(2*sqrt(b) - a);


%% ========================================================================
% 11. EXPRESIONES EN FUNCION DE LAS PROPIEDADES FISICAS
%% ========================================================================
%
% Recordando:
%
%       a = k/Ac
%
%       b = k/EI
%
% se obtiene:
%
%                        ___________________________
%                 1     / k          / k
%       alfa = --------V  --- + 2* V ---
%                 2       Ac          EI
%
%
%                        ___________________________
%                 1     /     / k       k
%       beta = --------V  2* V ---  -  ---
%                 2           EI       Ac
%
% Es decir:
%
% alfa = 1/2*sqrt(k/Ac + 2*sqrt(k/EI))
%
% beta = 1/2*sqrt(2*sqrt(k/EI) - k/Ac)
%
% Para que beta sea real debe cumplirse:
%
%       2*sqrt(k/EI) > k/Ac
%
% o equivalentemente:
%
%       (k/Ac)^2 < 4*k/EI
%
% que es exactamente la condicion:
%
%       a^2 - 4*b < 0
%


%% ========================================================================
% 12. RAICES DE LA ECUACION CARACTERISTICA
%% ========================================================================
%
% Con alfa y beta definidos anteriormente, las cuatro raices son:
%
%       r1 =  alfa + i*beta
%       r2 =  alfa - i*beta
%       r3 = -alfa + i*beta
%       r4 = -alfa - i*beta
%
% Es decir:
%
%       r = +/- alfa +/- i*beta
%
% Las raices aparecen en pares conjugados debido a que la ecuacion
% diferencial tiene coeficientes reales.
%


%% ========================================================================
% 13. SOLUCION HOMOGENEA REAL
%% ========================================================================
%
% A partir de las cuatro raices anteriores puede escribirse:
%
% vh = exp(alfa*x)*(C1*cos(beta*x) + C2*sin(beta*x)) + ...
%      exp(-alfa*x)*(C3*cos(beta*x) + C4*sin(beta*x));
%
% Esta expresion es completamente real.
%
% Una forma equivalente y conveniente para desarrollar posteriormente
% la matriz fundamental es:
%

vh = C1*cosh(alfa*x)*cos(beta*x) + ...
     C2*cosh(alfa*x)*sin(beta*x) + ...
     C3*sinh(alfa*x)*cos(beta*x) + ...
     C4*sinh(alfa*x)*sin(beta*x);

% NOTA:
%
% Las cuatro funciones independientes son:
%
% f1 = cosh(alfa*x)*cos(beta*x)
% f2 = cosh(alfa*x)*sin(beta*x)
% f3 = sinh(alfa*x)*cos(beta*x)
% f4 = sinh(alfa*x)*sin(beta*x)
%
% Es importante que el ultimo termino sea:
%
%       sinh(alfa*x)*sin(beta*x)
%
% y NO:
%
%       cosh(alfa*x)*sin(beta*x)
%
% porque esta ultima ya aparece en f2 y se perderia independencia lineal.


%% ========================================================================
% 14. SOLUCION PARTICULAR PARA CARGA TRAPEZOIDAL
%% ========================================================================
%
% La carga es:
%
%       q(x) = q1 + (q2-q1)*x/L
%
% Buscamos una solucion particular:
%
%       vp = A + B*x
%
% Como vp es lineal:
%
%       vp''   = 0
%
%       vp'''' = 0
%
% Sustituyendo en:
%
%       EI*vp'''' - EI*k/Ac*vp'' + k*vp = q(x)
%
% resulta:
%
%       k*vp = q(x)
%
% Por tanto:
%
%       vp = q(x)/k
%

q = q1 + (q2-q1)*x/L;

vp = q/k;


%% ========================================================================
% 15. SOLUCION GENERAL DEL DESPLAZAMIENTO
%% ========================================================================
%
% La solucion completa es:
%
%       v = vh + vp
%

v = C1*cosh(alfa*x)*cos(beta*x) + ...
    C2*cosh(alfa*x)*sin(beta*x) + ...
    C3*sinh(alfa*x)*cos(beta*x) + ...
    C4*sinh(alfa*x)*sin(beta*x) + ...
    (q1 + (q2-q1)*x/L)/k;


%% ========================================================================
% 16. RECUPERACION DEL MOMENTO FLECTOR
%% ========================================================================
%
% De la deduccion inicial:
%
%       M = EI*(v'' - k*v/Ac + q/Ac)
%

M = EI*(diff(v,x,2) - k*v/Ac + q/Ac);


%% ========================================================================
% 17. RECUPERACION DEL CORTANTE
%% ========================================================================
%
% Como:
%
%       M' = V
%
% entonces:
%

V = diff(M,x);

% Equivalentemente:
%
%       V = EI*(v''' - k*v'/Ac + q'/Ac)
%
% Para la carga trapezoidal:
%
%       q' = (q2-q1)/L
%
% luego:
%
%       V = EI*(v''' - k*v'/Ac + (q2-q1)/(Ac*L))
%


%% ========================================================================
% 18. RECUPERACION DEL GIRO
%% ========================================================================
%
% La ecuacion cinematica original es:
%
%       v' = t - V/Ac
%
% Despejando t:
%
%       t = v' + V/Ac
%
% Por tanto:
%

t = diff(v,x) + V/Ac;


%% ========================================================================
% 19. COMPROBACION: LIMITE DE EULER-BERNOULLI
%% ========================================================================
%
% Si la rigidez a cortante tiende a infinito:
%
%       Ac -> infinito
%
% entonces:
%
%       a = k/Ac -> 0
%
% y:
%
%       alfa = 1/2*sqrt(2*sqrt(k/EI))
%
%       beta = 1/2*sqrt(2*sqrt(k/EI))
%
% Por tanto:
%
%       alfa = beta
%
% Simplificando:
%
%                         1/4
%                 1      / k
%       alfa = -------- *|---
%               sqrt(2)  \ EI
%
% Como:
%
%                     1/4
%                  / k
%       lambda =  |------
%                 \ 4*EI
%
% resulta:
%
%       alfa = beta = lambda
%
% recuperandose la solucion clasica de una viga de Euler-Bernoulli
% sobre fundacion elastica de Winkler.
%
%% ========================================================================




