clc
clear
syms C1 C2 C3 C4 C5 C6 x A B C D k Ac EI L M k1 re ima AE xi
E = 24870062.3;                 % Modulo de elasticidad del concreto [kPa]
G = 0.4*E;                      % Modulo de cortante [kPa]
L = 4;                          % Longitud [m]
Ae = 0.4^2;                     % Area [m^2]
I = 0.4^4/12;                   % Inercia [m^4]
Ac = Ae*G*5/6;                  % Coeficiente de correcci�n por cortante
k = 500;                        % Coeficiente de balasto
EI = E*I;                       % Rigidez a flexi�n
AE = Ae*E;                      % Rigidez axial
rho = 2.4;                      % Densidad del concreto [Mg/m^3]
g = 9.8066502;                  % Aceleraci�n de gravedad
P = 1;                          % Carga axial de pandeo [kN/m]
q = 0;                          % Carga vertical  o para matriz de rigidez[kN/m]
b = 0;                          % Carga axial o para matriz de rigidez[kN/m]
%A=((EI*k - (EI*k*(- 4*Ac^2 + EI*k))^(1/2))/(2*Ac*EI))^(1/2);
%B=((EI*k + (EI*k*(- 4*Ac^2 + EI*k))^(1/2))/(2*Ac*EI))^(1/2);
%C=(EI*k - (EI*k*(- 4*Ac^2 + EI*k))^(1/2))/(2*Ac) - (EI*k)/Ac;
%D=(EI*k + (EI*k*(- 4*Ac^2 + EI*k))^(1/2))/(2*Ac) - (EI*k)/Ac;
%M=C3*C*exp( x*A)...
%+ C4*C*exp(-x*A)... 
%+ C1*D*exp( x*B)...
%+ C2*D*exp(-x*B);
%syms Ac k y EI

%despejo todo en funcion de M
%V' == -k*v+q    (1) %Equilibrio de fuerzas cortantes
%M' == V       (2)  %%Relaci�n momento-cortante
%t' == M/EI    (3)  %Relaci�n giro-momento
%v' == t - V/Ac(4)  %Relaci�n desplazamiento-giro

%derivo (2)
%M''=V'
%(M''-q)/k=-v  (5)  (2) en (1)

%derivo de nuevo (5)

%(M^3-q')/k=-v' (7)

% despejo t en (4)

%t = v'+ V/Ac(6)
%remplazo  (7) y (4)

%t= (M^3-q')/k+M'/Ac

%derivo (6)

%t^1=(M^4-q'')/k+M^2/Ac

%M/EI=(M^4-q'')/k+M^2/Ac
%A=(solve(M^4/k-M^2/Ac+1/EI==0,M))

%A=(solve(M^4-k*M^2/Ac+k/EI==0,M))

%a=1
%b=-k/Ac
%c=k/EI

%m=(-b/(2*a)-(b^2-4*a*c)^(1/2)/(2*a))^(1/2)

%m=(k/(2*Ac) - (k^2/Ac^2 - (4*k)/EI)^(1/2)/2)^(1/2)

%m1=(k/(2*Ac)-(k^2/(Ac^2*4)-(k)/(EI))^(1/2))^(1/2)

re=k/(2*Ac);
ima=-(-k^2/(Ac^2*4)+(k)/(EI))^(1/2);

n=2;

r=(re^2+ima^2)^(1/2);
phi=atan2(ima,re);

k1=0;
wr=r^(1/n)*(cos((phi+2*pi*k1)/n));
wi=-r^(1/n)*(sin((phi+2*pi*k1)/n));
%syms wr wi
%M=cos(wi*x)*(C1*exp(wr*x)+C2*exp(-wr*x))+sin(wi*x)*(C3*exp(wr*x)+C4*exp(-wr*x));
M=cos(wi*x)*(C1*cosh(wr*x)+C2*sinh(wr*x))+sin(wi*x)*(C3*cosh(wr*x)+C4*sinh(wr*x));

n1=-(-(k*(((- 4*Ac^2 + EI*k)/(EI*k))^(1/2) - 1))/(2*Ac))^(1/2);
n2=-( (k*(((- 4*Ac^2 + EI*k)/(EI*k))^(1/2) + 1))/(2*Ac))^(1/2);
n3= (-(k*(((- 4*Ac^2 + EI*k)/(EI*k))^(1/2) - 1))/(2*Ac))^(1/2);
n4= ( (k*(((- 4*Ac^2 + EI*k)/(EI*k))^(1/2) + 1))/(2*Ac))^(1/2);

%a=-(k/(2*Ac)- ((k/(2*Ac))^2-k/(EI))^(1/2) )^(1/2)

%M=C1*exp(n1*x)+C2*exp(n2*x)+C3*exp(n3*x)+C4*exp(n4*x);
EA=E*Ae;
V=diff(M,x,1);
v=-diff(V,x,1)/k;
t=V/Ac+diff(v,x,1);

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
<<<<<<< HEAD
q2 = 25;       % Carga vertical final [kN/m]
=======
q2 = 30;       % Carga vertical final [kN/m]
>>>>>>> 3f995df81f4da25d2f1d935fb3b75b7b51361fd1
nq=1;       % exponente carga vertical final viga
qx = (q2 - q1)/L^nq * x^nq + q1;  % Carga vertical variable (polin�mica)

b1=0.4;  %%% ancho de la viga
kWinkler=500; % coeficiente de balasto

naxi=1;     % exponenete carga axial final viga
b1a = 25;      % Carga axial inicial [kN/m]
<<<<<<< HEAD
b2a = 25;      % Carga axial final [kN/m]
bx = (b2a - b1a)/L^naxi * x^naxi + b1a;  % Carga axial variable  (polin�mica)

b = matlabFunction(subs(bx,x, L*(1+xi)/2),"Vars",{xi});
q = matlabFunction(subs(qx,x, L*(1+xi)/2),"Vars",{xi});
=======
b2a = 30;      % Carga axial final [kN/m]
bx = (b2a - b1a)/L^naxi * x^naxi + b1a;  % Carga axial variable  (polin�mica)

b = matlabFunction(subs(bx,x, L*(1+xi)/2));
q = matlabFunction(subs(qx,x, L*(1+xi)/2));
>>>>>>> 3f995df81f4da25d2f1d935fb3b75b7b51361fd1

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
          + (Nt0v(:, i) - Nv0vd(:, i)) * (Nt0v(:, i) - Nv0vd(:, i))' * wv(i) * Ac * dx_dxi; % T�rmino de cortante
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



