
%function M_TE_rhow=matriz_de_masa(sec,q1,q2,nq,b1,b2,nb,h1,h2,nh,G,E,tf,tw,L)

clc
clear 
% # Definción de variables
syms b2 b1 nb x h1 h2 nh  E G  tw L tf tw C1 C2 C3 C4 C5 C6 xi
g=9.8066502;      % aceleracion de la gravedad
rho = 2.4;       % [Mg/m^3] densidad del concreto
E=24870062;      %Modulo elaticidad viga
v=0.25;          % coeficiente de poisson
G=E/(2*(1+v));   %Modulo cortante

b1=0.3;          % ancho inicial viga
b2=0.2;          % ancho final viga
nb=1;            % exponente ancho final viga

h1=0.5;          % altura inicial viga
h2=0.8;          % altura final viga
nh=1;            % exponente altura viga

L=4;        % longitud viga
q1=25;      % carga vertical inicial viga
q2=30;      % carga vertical final viga
nq=1;       % exponente carga vertical final viga

b1a=25;     % carga axial inicial viga
b2a=30;     % carga axial final viga
naxi=1;     % exponenete carga axial final viga
tf=0.09;    % espesor aleta sen secciones I, C , T L, Tubular hueca cirdular y rectangular 
tw=0.06;    % espesor alma sen secciones I, C , T L, Tubular hueca cirdular y rectangular
%}
datos=20;   % n es la enecima derivada  de la funcion f evaluada en le punto a

sec=2;      % seccion trasnversal a evaluar
[Ax,Ix,As2x]=secciones_3(sec,b1,b2,nb,h1,h2,nh,G,E,tf,tw,L);
%% # Se define la carga: hace la ecuación diferencial homogénea
q=0;
b=0;
%se definen las ecuaciones diferenciales
V=C1; 
M=int(V,x)+C2;

t=int(taylor(M/Ix,x,L/2,'Order',datos),x)+C3;
v=int(t-taylor(V/As2x,x,L/2,'Order',datos),x)+C4;
A=int(b,x)+C5;
u=int(taylor(A/Ax,x,L/2,'Order',datos),x)+C6;
%# Se calcula la matrix de rigidez
%K_TE2 = zeros(6);
N_u2 = sym(zeros(1,6));
N_w2 = sym(zeros(1,6));
%N_t2 = sym(zeros(1,6));

for i = 1:6
    [c1,c2,c3,c4,c5,c6]=solve(subs(u,x,0)==(i==1),...
                              subs(v,x,0)==(i==2),...% con sus respectivas condiciones de frontera
                              subs(t,x,0)==(i==3),...
                              subs(u,x,L)==(i==4),...
                              subs(v,x,L)==(i==5),...
                              subs(t,x,L)==(i==6),...
                              [C1,C2,C3,C4,C5,C6]);
    N_u2(i) = subs(u,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2}); 
    N_w2(i) = subs(v,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2});
	%N_t2(i) = subs(t,{C1,C2,C3,C4,C5,C6,x},{c1,c2,c3,c4,c5,c6,L*(1+xi)/2});
end         
%# Se convierten las funciones de forma a coordenadas naturales,
%% Funciones de forma
dx_dxi = L/2;   % jacobiano de la transformacion isoparametrico
N1 = vertcat(N_u2,N_w2);
qx=(q2-q1)/L^nq*x^nq+q1;
%Ae = simplify(subs(Ax/E,x, L*(1+xi)/2));
%I =  simplify(subs(Ix/E,x, L*(1+xi)/2));
q =  simplify(subs(qx/g,x, L*(1+xi)/2));
%%  MATRIZ DE MASA
%M_TE_rhoA = double(int(rho*Ae*(N1'*N1)*dx_dxi,xi, -1, 1));
%N2= vertcat(N_t2);
%M_TE_rhoI = double(int(rho*I*(N2'*N2)*dx_dxi,xi, -1, 1));
%M=M_TE_rhoA+M_TE_rhoI;
%%  MATRIZ DE MASA carga externa
M_TE_rhow = double(int(q*(N1'*N1)*dx_dxi,xi, -1, 1));


