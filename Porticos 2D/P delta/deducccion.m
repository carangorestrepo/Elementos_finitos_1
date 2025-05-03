%% Soluciones Exactas para Pandeo y Efectos de Segundo Orden
%% en Columnas de Viga Timoshenko Deformables por Cortante
%% Basado en Análisis Estructural Matricial

% Definición de variables:
% V = fuerza cortante
% M = momento flector
% t = ángulo de giro
% P = carga axial para cálculo P-Delta
% ksGAs = valor para considerar el esfuerzo cortante
% D = operador diferencial
% Q = fuerza transversal equivalente

%% Ecuaciones fundamentales
V == Q - P*D*v;          % (1) Relación fuerza cortante
D*M == Q;                % (2) Equilibrio de momentos
D*V == -w;               % (3) Equilibrio de fuerzas verticales
M == EI*D*t;             % (4) Relación momento-curvatura
Q == ksGAs*(t + D*v);    % (5) Relación fuerza transversal-deformación

%% Derivación de ecuaciones adicionales
% De (1):
V*D == Q*D - P*D^2*v;    % (1.1) Primera derivada de (1)
D*Q == P*D^2*v + V*D;    % (1.2) Reorganización de (1.1)

% De (2):
D^2*M == Q*D;            % (2.1) Segunda derivada del momento
D^2*M ==  Q*D;            % (2.2) Equivalente a (2.1)

% De (4):
D*M == EI*D^2*t;         % (4.1) Primera derivada de (4)
D^2*M == EI*D^3*t;       % (4.2) Segunda derivada de (4)

% De (5):
D*Q == ksGAs*(D*t + D^2*v);          % (5.1) Primera derivada de (5)
D*t == (D*Q - ksGAs*D^2*v)/ksGAs;    % (5.2) Despeje de D*t

%% Sustituciones y simplificaciones
% Sustituyendo (1.1) y (3) en (5.2):
D*t == (P*D^2*v - w - ksGAs*D^2*v)/ksGAs;

% Derivadas sucesivas del ángulo:
D^2*t == (P*D^3*v - w - ksGAs*D^3*v)/ksGAs;
D^3*t == (P*D^4*v - w - ksGAs*D^4*v)/ksGAs;

% Sustituyendo (2.1) y (4.1) en (1.1):
-w == EI*D^3*t - P*D^2*v;

% Despejando EI*D^3*t:
EI*D^3*t == P*D^2*v - w;

%% Solución para carga distribuida nula (w = 0)
EI*(P*D^4*v - ksGAs*D^4*v)/ksGAs == P*D^2*v;

% Simplificando:
D^4*v*EI*(P - ksGAs)/ksGAs - P*D^2*v == 0;

% Definición del parámetro de pandeo (lambda):
lambda = sqrt(P/((1 - P/ksGAs)*EI));

%% Solución general para la deformación vertical
y = C1 + C2*x + sin(x*lambda)*C3 + cos(x*lambda)*C4;

% Expresiones para momento, cortante y ángulo:
M = EI*(P - ksGAs)/ksGAs*diff(y, x, 2);
V = D*M - P*D*y;
Q = D*M;
t = Q/ksGAs - D*y;



