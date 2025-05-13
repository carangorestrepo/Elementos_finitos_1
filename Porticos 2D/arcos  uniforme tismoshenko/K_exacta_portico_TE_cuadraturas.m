clc
clear
syms C1 C2 C3 C4 C5 C6 EI Ac AE x L xi ro A I
%se definen las ecuaciones diferenciales a flexion y cortante
%%{
E = 4700*sqrt(28)*1000;    % [kPa] modulo de elasticidad del concreto (NSR 10  - art C.8.5)
G = 0.4*E;                 % [kPa] modulo de cortante (NSR 10  - art CR8.8.2)
Ae = 0.4^2;                 % [m^2] area
I = 0.4^4/12;              % [m^4] inercia
L = 4;                     % [m]   longitud
ks = 5/6;                  % coeficiente de correccion del cortante para seccion rectangular
rho = 2.4;                 % [Mg/m^3] densidad del concreto
P=1;                       % [kN/m] carga axual de pandeo
g=9.8066502;      % aceleracion de la gravedad
Ac=G*ks*Ae;
EI=E*I;
AE=Ae*E;
%}
q=0;
b=0;
V=int(q,x)+C1;
M=int(V,x)+C2;
t=int(M/EI,x)+C3;
v=int(t-V/Ac,x)+C4;
%se definen las ecuaciones diferenciales a carga axial
A=int(b,x)+C5;
u=int(A/AE,x)+C6;

%# Se calcula la matrix de rigidez
K_TE2 = sym(zeros(6,6));
N_u2 = sym(zeros(1,6));
N_w2 = sym(zeros(1,6));
N_t2 = sym(zeros(1,6));
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
K_TE2=	double(K_TE2);

%% cuadratura de Gauss-Legendre
dx_dxi = L/2;              % jacobiano de la transformacion isoparametrica
wv=[ 0.347854845137454, 0.652145154862546,  0.652145154862546, 0.347854845137454];
xiv= [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053];
txi = 4; % polinomio grado 4

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


q1=25;      % carga vertical inicial viga
q2=30;      % carga vertical final viga
nq=1;       % exponente carga vertical final viga
qx=(q2-q1)/L^nq*x^nq+q1;
q = matlabFunction(subs(qx,x, L*(1+xi)/2));

b1=0.4;  %%% ancho de la viga
kWinkler=500; % coeficiente de balasto

naxi=1;     % exponenete carga axial final viga
b1a=25;     % carga axial inicial viga
b2a=30;     % carga axial final viga
bqx=(b2a-b1a)/L^naxi*x^naxi+b1a;
b = matlabFunction(subs(bqx,x, L*(1+xi)/2));

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
    k = k + Nu0vd(:,i)*Nu0vd(:,i)'*wv(i)*AE*dx_dxi + Nt0vd(:,i)*Nt0vd(:,i)'*wv(i)*EI*dx_dxi + (Nt0v(:,i)-Nv0vd(:,i))*(Nt0v(:,i)-Nv0vd(:,i))'*wv(i)*Ac*dx_dxi;
    
    % MATRIZ DE MASA ELEMENTO
    m = m + Nu0v(:,i)*Nu0v(:,i)'*wv(i)*rho*Ae*dx_dxi + Nv0v(:,i)*Nv0v(:,i)'*wv(i)*rho*Ae*dx_dxi  + Nt0v(:,i)*Nt0v(:,i)'*wv(i)*rho*I*dx_dxi;
    
    % MATRIZ DE CARGA EXTERNA
    mq = mq + Nu0v(:,i)*Nu0v(:,i)'*wv(i)*bxi/g*dx_dxi + Nv0v(:,i)*Nv0v(:,i)'*wv(i)*qxi/g*dx_dxi  + Nt0v(:,i)*Nt0v(:,i)'*wv(i)*rho*I*dx_dxi;
    
    % 	MOMENTOS DE EMPOTRAMIENTO
    MV = MV + Nu0v(:,i)*wv(i)*bxi*dx_dxi + Nv0v(:,i)*wv(i)*qxi*dx_dxi;
    
    % MATRIZ DE GEOMETRICA
    kG = kG + Nu0vd(:,i)*Nu0vd(:,i)'*wv(i)*P*dx_dxi + Nt0vd(:,i)*Nt0vd(:,i)'*wv(i)*I/Ae*dx_dxi ;
    
    %% Matriz de rigidez asociada a la cimentación elástica de Winkler

    H = H + Nu0v(:,i)*Nu0v(:,i)'*wv(i)*b1*kWinkler*dx_dxi+Nv0v(:,i)*Nv0v(:,i)'*wv(i)*b1*kWinkler*dx_dxi;
end





