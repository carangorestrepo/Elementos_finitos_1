clc
clear all
close all

%% DEFINICI�N DE VARIABLES SIMB�LICAS
syms D EI Ac P C1 C2 C3 C4 C5 C6 x L AE landa xi
syms bb ba wa wb  % Para cargas distribuidas

%% PROPIEDADES DEL MATERIAL Y GEOMETR�A
% Valores para hormig�n armado (ejemplo)
E = 24870062;      % M�dulo de elasticidad [kN/m�]
G = 0.4*E;         % M�dulo de cortante [kN/m�]
I = 0.4^4/12;      % Inercia de la secci�n [m?]
Ae = 0.4^2;        % �rea efectiva [m�]
EI = E*I;          % Rigidez a flexi�n [kN�m�]
Ac = Ae*G*5/6;     % Rigidez a cortante reducida [kN]
g = 9.8066502;     % Aceleraci�n de gravedad [m/s�]
rho = 2.4;         % Densidad del concreto [Mg/m�]

% Carga axial aplicada
P = 1000;          % Carga axial [kN]

% Longitud del elemento
L = 4;             % Longitud del elemento [m]
AE = Ae*E;         % Rigidez axial [kN]

%% PAR�METRO DE PANDEO (?)
X = 1 - P/Ac;      % Factor de reducci�n por cortante
landa = sqrt(P/(X*EI)); % Par�metro de pandeo [1/m]

%% SOLUCI�N DE LA ECUACI�N DIFERENCIAL
% Soluci�n para desplazamientos transversales (v) y axiales (u)
v = C1 + C2*x + sin(x*landa)*C3 + cos(x*landa)*C4; % Transversal
b = 0;               % Carga axial distribuida (inicialmente cero)
A = int(b,x) + C5;   % Fuerza axial
u = int(A/AE,x) + C6;% Desplazamiento axial

%% DEFINICI�N DE ESFUERZOS INTERNOS
M = -(1 - P/Ac)*diff(v,x,2)*EI;  % Momento flector
V = diff(M,x) - P*diff(v,x);     % Cortante
Q = diff(M,x);                   % Flujo cortante
t = Q/Ac - diff(v,x);            % Rotaci�n por cortante

%% CONSTRUCCI�N DE LA MATRIZ DE RIGIDEZ (6x6)
K_TE2 = zeros(6);    % Inicializaci�n matriz de rigidez tangente
N_u2 = sym(zeros(1,6)); % Funciones de forma axiales
N_w2 = sym(zeros(1,6)); % Funciones de forma transversales
N_t2 = sym(zeros(1,6)); % Funciones de forma rotacionales

for i = 1:6
    % Resuelve constantes para condiciones de contorno
    [c1,c2,c3,c4,c5,c6] = solve(...
        subs(u,x,0) == (i==1),...  % Desplaz. axial nodo 1
        subs(v,x,0) == (i==2),...  % Desplaz. transv. nodo 1
        subs(t,x,0) == (i==3),...  % Rotaci�n nodo 1
        subs(u,x,L) == (i==4),...  % Desplaz. axial nodo 2
        subs(v,x,L) == (i==5),...  % Desplaz. transv. nodo 2
        subs(t,x,L) == (i==6),...  % Rotaci�n nodo 2
        [C1,C2,C3,C4,C5,C6]);
    
    % Ensambla matriz de rigidez
    K_TE2(:,i) = [
        -subs(A,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0});  % Fuerza axial nodo 1
         subs(V,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0});  % Cortante nodo 1
        -subs(M,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,0});  % Momento nodo 1
         subs(A,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L});  % Fuerza axial nodo 2
        -subs(V,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L});  % Cortante nodo 2
         subs(M,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L})   % Momento nodo 2
    ];
    
    % Funciones de forma en coordenadas naturales
    N_t2(i) = subs(t,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2});
    N_w2(i) = subs(v,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2});
    N_u2(i) = subs(u,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2}); 
end

%% CORRECCI�N DE SIGNOS PARA GRADOS DE LIBERTAD TRANSVERSALES
for i = [2,5]
    K_TE2(2,i) = -K_TE2(2,i);
    K_TE2(3,i) = -K_TE2(3,i);
    K_TE2(5,i) = -K_TE2(5,i);
    K_TE2(6,i) = -K_TE2(6,i);
end
K_TE2 = double(K_TE2);  % Conversi�n a num�rico

%% INTEGRACI�N NUM�RICA (CUADRATURA DE GAUSS-LEGENDRE)
dx_dxi = L/2;           % Jacobiano de transformaci�n

% Puntos y pesos de integraci�n (grado 15)
txi = 15; 
[xiv,wv] = gausslegendre_quad(txi);

% Conversi�n a funciones evaluables
N_t21 = matlabFunction(N_t2);
N_w21 = matlabFunction(N_w2);
N_u21 = matlabFunction(N_u2);

% Derivadas de las funciones de forma
dN_t21 = matlabFunction(diff(N_t2,xi)/dx_dxi,"Vars",{xi});
dN_w21 = matlabFunction(diff(N_w2,xi)/dx_dxi,"Vars",{xi});
dN_u21 = matlabFunction(diff(N_u2,xi)/dx_dxi,"Vars",{xi});

%% INICIALIZACI�N DE MATRICES GLOBALES
k = zeros(6);    % Matriz de rigidez el�stica
kG = zeros(6);   % Matriz de rigidez geom�trica
m = zeros(6);    % Matriz de masa
H = zeros(6);    % Matriz de Winkler
mq = zeros(6);   % Matriz de masa equivalente
MV = zeros(6,1); % Vector de momentos de empotramiento

%% DEFINICI�N DE CARGAS
% Carga transversal (polin�mica)
q1 = 25;         % Carga inicial [kN/m]
q2 = 25;         % Carga final [kN/m]
nq = 1;          % Exponente de variaci�n
qx = (q2 - q1)/L^nq * x^nq + q1;

% Carga axial (polin�mica)
b1a = 25;        % Carga inicial [kN/m]
b2a = 25;        % Carga final [kN/m]
naxi = 1;        % Exponente de variaci�n
bx = (b2a - b1a)/L^naxi * x^naxi + b1a;

% Coeficiente de balasto
b1 = 0.4;        % Ancho de viga [m]
kWinkler = 500;  % Coeficiente [kN/m�]

% Conversi�n a funciones evaluables
b = matlabFunction(subs(bx,x, L*(1+xi)/2),"Vars",{xi});
q = matlabFunction(subs(qx,x, L*(1+xi)/2),"Vars",{xi});

%% INTEGRACI�N NUM�RICA PARA MATRICES DE ELEMENTO
for i=1:txi
    %% Funciones de forma Lagrangianas
    
    Nu0v(:,i) = N_u21(xiv(i));%% Axial
    Nv0v(:,i) = N_w21(xiv(i));%% cortante
    Nt0v(:,i) = N_t21(xiv(i));%% flexion
    
    %% Defino las matrices de deformacion
    
    Nu0vd(:,i) = dN_u21(xiv(i));%% Axial
    Nv0vd(:,i) = dN_w21(xiv(i));%% cortante
    Nt0vd(:,i) = dN_t21(xiv(i));%% flexion
    qxi = q(xiv(i));
    bxi = b(xiv(i));
    
    %% Integro las matrices con una cuadratura de Gauss-Legendre de orden 4      
    % MATRIZ DE RIGIDEZ
    k = k + Nu0vd(:, i) * Nu0vd(:, i)' * wv(i) * AE * dx_dxi ...  % T�rmino axial
          + Nt0vd(:, i) * Nt0vd(:, i)' * wv(i) * EI * dx_dxi ...  % T�rmino de flexi�n
          + (Nt0v(:, i) - Nv0vd(:, i)) * (Nt0v(:, i) - Nv0vd(:, i))' * wv(i) * Ac * dx_dxi;% ... % T�rmino de cortante
          %+  Nt0vd(:, i) * Nt0vd(:, i)' * wv(i) *P * dx_dxi ...  % T�rmino axial
          %- (Nt0v(:, i) - Nv0vd(:, i)) * (Nt0v(:, i) - Nv0vd(:, i))' * wv(i) * (P/Ac) * EI * dx_dxi; % Correcci�n por interacci�n cortante
    % MATRIZ DE MASA ELEMENTO

    m = m + Nu0v(:, i) * Nu0v(:, i)' * wv(i) * rho * Ae * dx_dxi ...  % Masa axial
          + Nv0v(:, i) * Nv0v(:, i)' * wv(i) * rho * Ae * dx_dxi ...  % Masa vertical
          + Nt0v(:, i) * Nt0v(:, i)' * wv(i) * rho * I * dx_dxi;      % Masa rotacional
      
    % MATRIZ DE CARGA EXTERNA
    mq = mq + Nu0v(:,i)*Nu0v(:,i)'*wv(i)*bxi/g*dx_dxi + Nv0v(:,i)*Nv0v(:,i)'*wv(i)*qxi/g*dx_dxi  + Nt0v(:,i)*Nt0v(:,i)'*wv(i)*rho*I*dx_dxi;
    
    % 	MOMENTOS DE EMPOTRAMIENTO
    MV = MV + Nu0v(:,i)*wv(i)*bxi*dx_dxi + Nv0v(:,i)*wv(i)*qxi*dx_dxi;
    
    % MATRIZ DE GEOMETRICA
    kG = kG + Nu0vd(:,i)*Nu0vd(:,i)'*wv(i)*P*dx_dxi + Nt0vd(:,i)*Nt0vd(:,i)'*wv(i)*I/Ae*dx_dxi ;
    
    %% Matriz de rigidez asociada a la cimentaci�n el�stica de Winkler

    % Matriz de Winkler (H)
    H = H + Nu0v(:, i) * Nu0v(:, i)' * wv(i) * b1 * kWinkler * dx_dxi ...  % T�rmino axial
          + Nv0v(:, i) * Nv0v(:, i)' * wv(i) * b1 * kWinkler * dx_dxi;     % T�rmino vertical
end

a=1


