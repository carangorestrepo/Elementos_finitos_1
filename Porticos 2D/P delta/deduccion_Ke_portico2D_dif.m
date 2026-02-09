clc
clear

% # Definción de variables simbólicas
syms x L A(x) V(x) M(x) t(x) v(x) u(x) EI EA Ac b xi k bb ba wb wa Ae I g;

%%{
%% Se describen las propiedades de los materiales
g = 9.8066502;              % [m/s²] aceleracion de la gravedad
E = 24870062.3;             % [kPa] modulo de elasticidad del concreto (NSR 10 - art C.8.5)
G = 0.4*E;                  % [kPa] modulo de cortante (NSR 10 - art CR8.8.2)
L = 4;                      % [m] longitud de la columna
Ae = 0.4^2;                 % [m²] area de la sección transversal (columna de 0.4x0.4 m)
I = 0.4^4/12;               % [m?] momento de inercia de la sección
Ac = Ae*G*5/6;              % coeficiente de correccion del cortante para seccion rectangular
P = 1000;                   % [kN] carga axial aplicada (efecto de pandeo)
%}
% Propiedades de rigidez
EA = Ae*E;                  % rigidez axial
EI = E*I;                   % rigidez a flexión

%# Se calcula la matriz de rigidez considerando efectos de segundo orden

% Inicialización de matrices
K_TE2 = sym(zeros(6));      % Matriz de rigidez 6x6 (simbólica)
N_u2 = sym(zeros(1,6));     % Funciones de forma para desplazamiento axial
N_w2 = sym(zeros(1,6));     % Funciones de forma para desplazamiento vertical  
N_t2 = sym(zeros(1,6));     % Funciones de forma para rotación

% Bucle para calcular cada columna de la matriz de rigidez
for i = 1:6
    % Resuelve el sistema de ecuaciones diferenciales para el i-ésimo grado de libertad
    sol = dsolve(...       
           V == diff(M) - P*diff(v),...         % Ecuación de equilibrio con efecto P-delta
           diff(V) == 0,...                   % Relación cortante-carga distribuida
           M(x) == EI*diff(t),   ...            % Relación momento-rotación (ley constitutiva)
           diff(M) == Ac*(t + diff(v)), ...     % Ecuación que considera deformación por cortante
           diff(A) == 0,        ...             % Fuerza axial constante (no hay carga axial distribuida)
           diff(u) == A/(EA),     ...           % Relación deformación axial-desplazamiento axial
           u(0) == (i==1),      ...            % Desplazamiento axial en x=0 = 1 si i=1
           v(0) == (i==2),      ...            % Desplazamiento vertical en x=0 = 1 si i=2
           t(0) == (i==3),      ...            % Rotación en x=0 = 1 si i=3
           u(L) == (i==4),      ...            % Desplazamiento axial en x=L = 1 si i=4
           v(L) == (i==5),      ...            % Desplazamiento vertical en x=L = 1 si i=5
           t(L) == (i==6));                    % Rotación en x=L = 1 si i=6
    
    % Almacena las funciones de forma
    N_t2(i) = sol.t;        % Función de forma para rotación
    N_w2(i) = sol.v;        % Función de forma para desplazamiento vertical
    N_u2(i) = sol.u;        % Función de forma para desplazamiento axial
    
    % Calcula la matriz de rigidez
    K_TE2(:,i) = [
        -subs(sol.A,x,0);   % Reacción horizontal en apoyo izquierdo (X1)
        subs(sol.V,x,0);    % Reacción vertical en apoyo izquierdo (Y1)
        -subs(sol.M,x,0);   % Momento de reacción en apoyo izquierdo (M1)
        subs(sol.A,x,L);    % Reacción horizontal en apoyo derecho (X2)
        -subs(sol.V,x,L);   % Reacción vertical en apoyo derecho (Y2)
        subs(sol.M,x,L)     % Momento de reacción en apoyo derecho (M2)
    ];
end
% Convierte la matriz simbólica a numérica
K_TE2 = double(K_TE2);