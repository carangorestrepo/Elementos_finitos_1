clc
clear

% # Definción de variables simbólicas
syms x L A(x) V(x) M(x) t(x) v(x) u(x) EI EA Ac b xi k bb ba wb wa Ae I g;
syms  b1loc b2loc q1loc q2loc P Pv Mv
%{
%% Se describen las propiedades de los materiales
g = 9.8066502;              % [m/s²] aceleracion de la gravedad
E = 24870062.3;             % [kPa] modulo de elasticidad del concreto (NSR 10 - art C.8.5)
G = 0.4*E;                  % [kPa] modulo de cortante (NSR 10 - art CR8.8.2)
L = 4;                      % [m] longitud de la columna
Ae = 0.4^2;                 % [m²] area de la sección transversal (columna de 0.4x0.4 m)
I = 0.4^4/12;               % [m?] momento de inercia de la sección
Ac = Ae*G*5/6;              % coeficiente de correccion del cortante para seccion rectangular
P = 1000;                   % [kN] carga axial aplicada

% Cargas distribuidas variables a lo largo de la columna
b1loc = 25;                 % carga distribuida en el extremo inferior axial
b2loc = 35;                 % carga distribuida en el extremo superior axial
q1loc = 25;                 % otra carga distribuida en extremo inferior flexion
q2loc = 35;                 % otra carga distribuida en extremo superior
% Propiedades de rigidez
EA = Ae*E;                  % rigidez axial
EI = E*I;                   % rigidez a flexión
%}



%# Se calcula la matrix de rigidez
% Cargas distribuidas variables linealmente a lo largo de la columna
qxloc = b1loc - (x*(q1loc - q2loc))/L;      % carga distribuida en dirección x
qyloc = b1loc - (x*(b1loc - b2loc))/L;      % carga distribuida en dirección y
qxloc=0;


% Vector de fuerzas en los extremos (6 grados de libertad)
F_TE2 = (zeros(6,1));

% Resolución del sistema de ecuaciones diferenciales para pandeo
% El sistema incluye ecuaciones de equilibrio, compatibilidad y constitutivas
sol = dsolve(...       
       V == diff(M) - P*diff(v),...         % Ecuación de equilibrio vertical (incluye efecto P-delta)
       diff(V) == - qxloc,...               % Derivada de la fuerza cortante
       M(x) == EI*diff(t),   ...            % Relación momento-rotación (ley constitutiva)
       diff(M) == Ac*(t + diff(v)), ...     % Ecuación que relaciona derivada del momento con cortante
       diff(A) == qyloc,        ...         % Derivada de la fuerza axial
       diff(u) == A/(EA),     ...           % Relación deformación axial-fuerza axial
       u(0) == 0,      ...                  % Desplazamiento axial en x=0 es cero
       v(0) == 0,      ...                  % Desplazamiento vertical en x=0 es cero
       t(0) == 0,      ...                  % Rotación en x=0 es cero
       A(L) == 0,      ...                  % Desplazamiento axial en x=L es cero
       V(L) == 0,      ...                  % Desplazamiento vertical en x=L es cero   
       M(L) == Mv);                          % Rotación en x=L es cero
% # Se evaluan las reacciones en los apoyos
F_TE2(:,1) = double([
    -subs(sol.A,x,0);    % Reacción horizontal en apoyo izquierdo (X1)
    subs(sol.V,x,0);     % Reacción vertical en apoyo izquierdo (Y1)  
    -subs(sol.M,x,0);    % Momento de reacción en apoyo izquierdo (M1)
    subs(sol.A,x,L);     % Reacción horizontal en apoyo derecho (X2)
    -subs(sol.V,x,L);    % Reacción vertical en apoyo derecho (Y2)
    subs(sol.M,x,L)      % Momento de reacción en apoyo derecho (M2)
]);
