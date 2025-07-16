clc
clear all
close all

%% DEFINICIÓN DE VARIABLES SIMBÓLICAS
syms D EI Ac P C1 C2 C3 C4 C5 C6 x L AE landa xi
syms bb ba wa wb  % Para cargas distribuidas

%% PROPIEDADES DEL MATERIAL Y GEOMETRÍA
% Valores típicos para hormigón armado (ejemplo)
E = 24870062;      % Módulo de elasticidad [kN/m²]
G = 0.4*E;         % Módulo de cortante [kN/m²]
I = 0.4^4/12;      % Inercia de la sección [m?]
Ae = 0.4^2;        % Área efectiva [m²]
EI = E*I;          % Rigidez a flexión [kN·m²]
Ac = Ae*G*5/6;     % Rigidez a cortante reducida [kN]
g = 9.8066502;     % Aceleración de gravedad
rho = 2.4;         % Densidad del concreto [Mg/m^3]

% Carga axial aplicada
P = 1000;          % Carga axial [kN]

% Longitud del elemento
L = 4;             % Longitud del elemento [m]
AE = Ae*E;         % Rigidez axial [kN]

%% PARÁMETRO DE PANDEO (?)
% Calcula el parámetro que caracteriza el efecto P-Delta
X = 1 - P/Ac;      % Factor de reducción por cortante
landa = sqrt(P/(X*EI)); % Parámetro de pandeo [1/m]

%% SOLUCIÓN DE LA ECUACIÓN DIFERENCIAL
% Solución general para desplazamientos transversales
% La ecuación considera efectos de flexión, cortante y carga axial
v = C1 + C2*x + sin(x*landa)*C3 + cos(x*landa)*C4;

%% DEFINICIÓN DE ESFUERZOS INTERNOS
% Momento flector considerando efectos P-Delta
M = -(1 - P/Ac)*diff(v,x,2)*EI;

% Fuerza cortante (derivada del momento menos efecto de la carga axial)
V = diff(M,x) - P*diff(v,x);

% Flujo cortante
Q = diff(M,x);

% Rotación por cortante
t = Q/Ac - diff(v,x);

%% CONSTRUCCIÓN DE LA MATRIZ DE RIGIDEZ (6x6)
K_TE2 = zeros(6);   % Inicialización matriz de rigidez tangente
N_w2 = sym(zeros(1,6)); % Funciones de forma para desplazamiento
N_t2 = sym(zeros(1,6)); % Funciones de forma para rotación

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
                              subs(v,x,0)==(i==2),...% con sus respectivas condiciones de frontera
                              subs(t,x,0)==(i==3),...
                              subs(u,x,L)==(i==4),...
                              subs(v,x,L)==(i==5),...
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
	N_w2(i) = subs(v,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2});
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

%% cuadratura de Gauss-Legendre
dx_dxi = L/2;              % jacobiano de la transformacion isoparametrica
%wv=[ 0.347854845137454, 0.652145154862546,  0.652145154862546, 0.347854845137454];
%xiv= [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053];

txi = 15; % polinomio grado 4
[xiv,wv] = gausslegendre_quad(txi);

N_t21=matlabFunction(N_t2);
N_w21=matlabFunction(N_w2);
N_u21=matlabFunction(N_u2);

dN_t21=matlabFunction(diff(N_t2,xi)/dx_dxi,"Vars",{xi});
dN_w21=matlabFunction(diff(N_w2,xi)/dx_dxi,"Vars",{xi});
dN_u21=matlabFunction(diff(N_u2,xi)/dx_dxi,"Vars",{xi});

Nu0v = zeros(6,txi);
Nv0v = zeros(6,txi);
Nt0v = zeros(6,txi);
Nu0vd = zeros(6,txi);
Nv0vd = zeros(6,txi);
Nt0vd = zeros(6,txi);

qxi=zeros(1,txi);
bxi=zeros(1,txi);

k = zeros(6);
kG = zeros(6);
m = zeros(6);
H = zeros(6);
mq=zeros(6);

MV = zeros(6,1);

q1 = 25;       % Carga vertical inicial [kN/m]
q2 = 25;       % Carga vertical final [kN/m]
nq=1;       % exponente carga vertical final viga
qx = (q2 - q1)/L^nq * x^nq + q1;  % Carga vertical variable (polinómica)

b1=0.4;  %%% ancho de la viga
kWinkler=500; % coeficiente de balasto

naxi=1;     % exponenete carga axial final viga
b1a = 25;      % Carga axial inicial [kN/m]
b2a = 25;      % Carga axial final [kN/m]
bx = (b2a - b1a)/L^naxi * x^naxi + b1a;  % Carga axial variable  (polinómica)

b = matlabFunction(subs(bx,x, L*(1+xi)/2),"Vars",{xi});
q = matlabFunction(subs(qx,x, L*(1+xi)/2),"Vars",{xi});

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
    k = k + Nu0vd(:, i) * Nu0vd(:, i)' * wv(i) * AE * dx_dxi ...  % Término axial
          + Nt0vd(:, i) * Nt0vd(:, i)' * wv(i) * EI * dx_dxi ...  % Término de flexión
          + (Nt0v(:, i) - Nv0vd(:, i)) * (Nt0v(:, i) - Nv0vd(:, i))' * wv(i) * Ac * dx_dxi; % Término de cortante
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
    
    %% Matriz de rigidez asociada a la cimentación elástica de Winkler

    % Matriz de Winkler (H)
    H = H + Nu0v(:, i) * Nu0v(:, i)' * wv(i) * b1 * kWinkler * dx_dxi ...  % Término axial
          + Nv0v(:, i) * Nv0v(:, i)' * wv(i) * b1 * kWinkler * dx_dxi;     % Término vertical
end

a=1
