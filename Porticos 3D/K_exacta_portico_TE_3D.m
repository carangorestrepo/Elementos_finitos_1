clc
clear
% # Definción de variables
syms x L  Ax(x) Vz(x) My(x) tz(x) vz(x) ux(x)      EIy EA Acz...
          T(x) Vy(x) Mz(x) ty(x) vy(x) uz(x) GJ   EIz    Acy;
%% Se describen las propiedades de los materiales
by=0.3;
hz=0.4;
E=24855578;
G=E*0.4;
A=by*hz;
Iy=by*hz^3/12;
Iz=hz*by^3/12;
Acy=5/6*G*A;
Acz=5/6*G*A;
EIy=E*Iy;
EIz=E*Iz;
EA=E*A;
L=6;   
bj=min([by,hz]);
hj=max([by,hz]);
J=(1/3-0.21*bj/hj*(1-1/12*(bj/hj)^4))*hj*bj^3;% coeficiente de correccion del cortante para seccion rectangular
GJ=G*J;
%# Se calcula la matrix de rigidez
K_TE2 = sym(zeros(12));
N_u2x = sym(zeros(1,12));
N_u2z = sym(zeros(1,12));
N_w2y = sym(zeros(1,12));
N_w2z = sym(zeros(1,12));
N_t2y = sym(zeros(1,12));
N_t2z = sym(zeros(1,12));
for i = 1:12
    sol = dsolve(...       
           diff(Vz,x)==0,            ...%diff(V,x)+k*v==0 se definen las ecuaciones diferenciales
           diff(My,x) == Vz,         ...%diff(M) == V-P*diff(v,x)
           diff(tz,x) == My/(EIy),   ... 
           diff(vz,x) == tz-Vz/(Acz),...%diff(v) == t-V/(Ac)
           diff(Vy,x)==0,            ...%diff(V,x)+k*v==0 se definen las ecuaciones diferenciales
           diff(Mz,x) == Vy,         ...%diff(M) == V-P*diff(v,x)
           diff(ty,x) == Mz/(EIz),   ... 
           diff(vy,x) == ty-Vy/(Acy),...%diff(v) == t-V/(Ac)
           diff(ux,x) == Ax/(EA),  ...
           diff(Ax,x) ==0,         ...  
           diff(uz,x) == T/(GJ),     ...
           diff(T,x) ==0,            ...
           ux(0) == (i==1),         ...%Axialx
           vy(0) == (i==2),         ...def y
           vz(0) == (i==3),         ...def z
           uz(0) == (i==4),         ...%torsion x
           tz(0) == (i==5),         ...giro z
           ty(0) == (i==6),         ...giro y
           ux(L) == (i==7),         ...%Axialx
           vy(L) == (i==8),         ...def y
           vz(L) == (i==9),         ...def z
           uz(L) == (i==10),        ...%torsion
           tz(L) == (i==11),        ...giro z
           ty(L) == (i==12)         ...giro y
           ); 
   N_t2z(i) = sol.tz;
   N_w2z(i) = sol.vz;
   N_u2x(i) = sol.ux;
   N_t2y(i) = sol.ty;
   N_w2y(i) = sol.vy;
   N_u2z(i) = sol.uz;
   
   % # se evaluan las reacciones horizontales y verticales y los momentos en los apoyos
       K_TE2(:,i)=[-subs(sol.Ax,x,0);
                   subs(sol.Vy,x,0);
                   -subs(sol.Vz,x,0);
                   -subs(sol.T,x,0);
                   subs(sol.My,x,0);
                   -subs(sol.Mz,x,0);
                   subs(sol.Ax,x,L);
                   -subs(sol.Vy,x,L);
                   subs(sol.Vz,x,L);
                   subs(sol.T,x,L);
                   -subs(sol.My,x,L);
                   subs(sol.Mz,x,L)];
end
K1=double(K_TE2)
%# Se convierten las funciones de forma a coordenadas naturales,
%% Funciones de forma
N_u2x = simplify(subs(N_u2x,x, L*(1+xi)/2));
N_w2z = simplify(subs(N_w2z,x, L*(1+xi)/2));  %# x = L*xi/2 + L/2
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



