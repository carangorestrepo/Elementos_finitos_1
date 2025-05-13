clc
clear
% # Definción de variables
syms x L A(x) V(x) M(x) t(x) v(x) u(x) EI EA Ac b  xi k;
%% Se describen las propiedades de los materiales
E = 4700*sqrt(28)*1000;       % [kPa] modulo de elasticidad del concreto (NSR 10  - art C.8.5)
G = 0.4*E;                    % [kPa] modulo de cortante (NSR 10  - art CR8.8.2)
L = 4;                        % [m]   longitud
Ae = 0.4^2;                   % [m^2] area  
I = 0.4^4/12;                 % [m^4] inercia
Ac=Ae*G*5/6;                  % coeficiente de correccion del cortante para seccion rectangular


%% # Se define la carga: hace la ecuación diferencial homogénea
q =0;  
m = 0;
EA=Ae*E;
EI=E*I;
%# Se calcula la matrix de rigidez
K_TE2 = sym(zeros(6));
N_u2 = sym(zeros(1,6));
N_w2 = sym(zeros(1,6));
N_t2 = sym(zeros(1,6));
for i = 1:6
    sol = dsolve(...       
           diff(V)+k*v == q,        ... % diff(V,x)+k*v==q se definen las ecuaciones diferenciales
           diff(M) == V - m,    ... %-P*diff(v,x),    ...
           diff(t) == M/(EI),   ... 
           diff(v) == t - V/Ac, ...
           diff(A) == 0,        ...
           diff(u) == A/EA,     ...  
           u(0) == (i==1),      ...
           v(0) == (i==2),      ... % con sus respectivas condiciones de frontera
           t(0) == (i==3),      ...
           u(L) == (i==4),      ...
           v(L) == (i==5),      ...   
           t(L) == (i==6)); 
   N_t2(i) = sol.t;
   N_w2(i) = sol.v;
   N_u2(i) = sol.u; 
   % # se evaluan las reacciones horizontales y verticales y los momentos en los apoyos
       K_TE2(:,i)=[-subs(sol.A,x,0); % X1
                    subs(sol.V,x,0); % Y2
                   -subs(sol.M,x,0); % M2
                    subs(sol.A,x,L); % X2 
                   -subs(sol.V,x,L); % Y2
                    subs(sol.M,x,L)];% M2
end

%# Se convierten las funciones de forma a coordenadas naturales,
%% Funciones de forma
N_u2 = simplify(subs(N_u2,x, L*(1+xi)/2));
N_w2 = simplify(subs(N_w2,x, L*(1+xi)/2));  %# x = L*xi/2 + L/2
N_t2 = simplify(subs(N_t2,x, L*(1+xi)/2));
N1 = vertcat(N_u2,N_w2);

%% Derivadas funciones de forma 
dx_dxi = L/2;   % jacobiano de la transformacion isoparametrico
dN_u2 = diff(N_u2,1,xi)/dx_dxi;
dN_w2 = diff(N_w2,1,xi)/dx_dxi;
dN_t2 = diff(N_t2,1,xi)/dx_dxi;
% # Definción de variables
syms   rho  b w P;
rho = 2.4;            % [Mg/m^3] densidad del concreto
b = 9;                % [kN/m] carga uniforme distribuida axial
wa = 25;              % [kN/m] carga uniforme distribuida vertical inicial
wb = 40;              % [kN/m] carga uniforme distribuida vertical final
P=1000;               % [kN/m] carga axual de pandeo
%%  MATRIZ DE MASA

M_TE_rhoA = double(int(rho*Ae*(N1'*N1)*dx_dxi,xi, -1, 1));
N2= vertcat(N_t2);
M_TE_rhoI = double(int(rho*I*(N2'*N2)*dx_dxi,xi, -1, 1));
M=M_TE_rhoA+M_TE_rhoI;

%% MATRIZ DE RIGIDEZ
K = double(int(dN_u2'*dN_u2*EA*dx_dxi + dN_w2'*dN_t2*EI*dx_dxi + (N_t2-dN_w2)'*(N_t2-dN_w2)*Ac*dx_dxi,xi,-1,1));

%% MOMENTOS DE EMPOTRAMIENTO
w =(wb-wa)*x/L + wa;       % carga trapezoidal
wf = simplify(subs(w,x, L*(1+xi)/2));
MV =  double(int(N_u2'*b*dx_dxi + N_w2'*wf*dx_dxi,xi,-1,1));

%% MATRIZ GEOMETRICA
%A Unified Approach to the Timoshenko Geometric Stiffness Matrix
%Considering Higher-Order Terms in the Strain Tensor
kga=double(int(P*(dN_w2'*dN_w2)*dx_dxi+P*I/Ae*(dN_t2'*dN_t2)*dx_dxi,-1,1));
%% Matriz de rigidez asociada a la cimentación elástica de Winkler
% b es el ancho de la viga
% kWinkler es el coeficiente de balasto
syms b kWinkler
b=0.4;
kWinkler=500;
H = double(b*kWinkler*simplify(int(N1'*N1*dx_dxi,xi,-1,1)));
%disp('H = b*kWinkler*'); pretty(H/(b*kWinkler));


