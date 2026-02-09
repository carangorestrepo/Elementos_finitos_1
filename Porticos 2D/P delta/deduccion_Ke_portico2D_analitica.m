clc
clear
syms EI P ksGA x k C1 C2 C3 C4 C5 C6 xi
%% PROPIEDADES DEL MATERIAL Y GEOMETRÍA
% Valores típicos para hormigón armado (ejemplo)
E = 24870062;      % Módulo de elasticidad [kN/m²]
G = 0.4*E;         % Módulo de cortante [kN/m²]
I = 0.4^4/12;      % Inercia de la sección [m?]
Ae = 0.4^2;        % Área efectiva [m²]
EI = E*I;          % Rigidez a flexión [kN·m²]
ksGA = Ae*G*5/6;    % Rigidez a cortante reducida [kN]
g = 9.8066502;     % Aceleración de gravedad
rho = 2.4;         % Densidad del concreto [Mg/m^3]

% Carga axial aplicada
P = 1000;          % Carga axial [kN]

% Longitud del elemento
L = 4;             % Longitud del elemento [m
AE = Ae*E;         % Rigidez axial [kN]

% Factor ?
chi = 1 - P/ksGA;

% Solución homogénea para y
k_val = sqrt(P/(EI*chi));  % = ?/L
y = C1 + C2*x + C3*cos(k_val*x) + C4*sin(k_val*x);

% Derivadas
dy = diff(y, x);
d2y = diff(y, x, 2);
d3y = diff(y, x, 3);
% Momento M
M = -EI*chi*d2y;
% Cortante Q
Q = -EI*chi*d3y;
% Cortante V
V = Q - P*dy;
% Rotación ?
t = Q/ksGA - dy;

%se definen las ecuaciones diferenciales a carga axial
b=0;
A=int(b,x)+C5;
u=int(A/AE,x)+C6;

%# Se calcula la matrix de rigidez
K_TE2 = zeros(6);   % Inicialización matriz de rigidez tangente
N_u2 = sym(zeros(1,6));
N_w2 = sym(zeros(1,6)); % Funciones de forma para desplazamiento
N_t2 = sym(zeros(1,6)); % Funciones de forma para rotación

for i = 1:6
    [c1,c2,c3,c4,c5,c6]=solve(subs(u,x,0)==(i==1),...
                              subs(y,x,0)==(i==2),...% con sus respectivas condiciones de frontera
                              subs(t,x,0)==(i==3),...
                              subs(u,x,L)==(i==4),...
                              subs(y,x,L)==(i==5),...
                              subs(t,x,L)==(i==6),...
                              [C1,C2,C3,C4,C5,C6]);
    % # se evaluan las reacciones horizontales y verticales y los momentos en los apoyos
	K_TE2(:,i)=[-subs(A,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0}); % X1
                 subs(V,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0}); % Y2
                -subs(M,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0}); % M2
                 subs(A,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L}); % X2 
                -subs(V,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L}); % Y2
                 subs(M,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L})];% M2
	N_t2(i) = subs(t,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2});
	N_w2(i) = subs(y,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2});
	N_u2(i) = subs(u,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2}); 
end

K_TE2(2,2)=-K_TE2(2,2);
K_TE2(3,2)=-K_TE2(3,2);
K_TE2(5,2)=-K_TE2(5,2);
K_TE2(6,2)=-K_TE2(6,2);

K_TE2(2,5)=-K_TE2(2,5);
K_TE2(3,5)=-K_TE2(3,5);
K_TE2(5,5)=-K_TE2(5,5);
K_TE2(6,5)=-K_TE2(6,5);
K_TE2=	double(K_TE2);