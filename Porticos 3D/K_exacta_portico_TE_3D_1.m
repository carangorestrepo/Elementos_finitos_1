
% Programa para deducir la matriz de rigidez de un elemento de viga de
% tismoshenko a partir de la solucion de la ecuacion diferencial

clc
clear

%% Definición de variables
syms x L EIz EIy EA GJ Acy Acz b bt 
% Constantes de integración
syms C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12;
%%{
dy=0.4;          % [m]   altura elemento
bfx=0.55;        % [m]   ancho elemento
tw=0;            % [m]   ancho alma
tf=0;            % [m]   ancho aleta
L=4;             % [m]   longitud

%% secciones tranvesales a analizar
sec=2;

%Channel        1
%Rectangular    2
%Circle         3
%I/Wide Flange  4
%Pipe             Seccion circular hueca 5
%Box/Tube         Tubular hueca          6
%Angle            Perfil L               7
%Tee              Perfil T               8

[A,Ix,Iy,As2,As3,J]=propiedades_geometrica_perfiles_v1(dy,bfx,tw,tf);
%% 
E=24855578.06;   % [kPa] modulo de elasticidad del concreto (NSR 10  - art C.8.5)
G=10356490.86;   % [kPa] modulo de cortante (NSR 10  - art CR8.8.2)
EIz=Ix(2)*E;     % [m^4*kPa] inercia alrededor del eje x por modulo de eliticidad
EIy=Iy(2)*E;     % [m^4*kPa] inercia alrededor del eje y por modulo de eliticidad
Acy=As2(2)*G;    % [m^2*kPa]coeficiente de correccion del cortante para seccion por modulo de cortante
Acz=As3(2)*G;    % [m^2*kPa]coeficiente de correccion del cortante para seccion por modulo de cortante
EA=E*A(2);       % [m^2*kPa] area de seccion transvesal por modulo de eliticidad
GJ=G*J(2);       % [m^4*kPa] constante torsional por modulo de cortante
%}
b=0;
q=0;
bt=0;
%% Se calcula la matrix de rigidez

% se integran las ec. diferenciales
A = int(b,x)+C1;
u = int(A/EA,x)+C2;

T = int(bt,x)+C3;
tt = int(T/GJ,x)+C4;

V2=int(q,x)+C5;
M3=int(V2,x)+C6;
t3=int(M3/EIy,x)+C7;
v3=int(t3-V2/Acz,x)+C8;

V3=int(q,x)+C9;
M2=int(V3,x)+C10;
t2=int(M2/EIz,x)+C11;
v2=int(t2-V3/Acy,x)+C12;

K_EE = sym(zeros(12));

for i = 1:12
      [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12]=solve(subs(u ,x,0) == (i==1),...%%# Condiciones de frontera
                                                     subs(v2,x,0) == (i==2),...
                                                     subs(v3,x,0) == (i==3),...
                                                     subs(tt,x,0) == (i==4),...
                                                     subs(t3,x,0) == (i==5),...
                                                     subs(t2,x,0) == (i==6),...
                                                     subs(u ,x,L) == (i==7),...
                                                     subs(v2,x,L) == (i==8),...
                                                     subs(v3,x,L) == (i==9),...
                                                     subs(tt,x,L) == (i==10),...
                                                     subs(t3,x,L) == (i==11),...
                                                     subs(t2,x,L) == (i==12),...
                                                     [C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12]);                                         
        K_EE(:,i)=[-subs(A, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0}); 
                    subs(V3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0}); 
                    subs(V2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                   -subs(T, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                    subs(M3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                   -subs(M2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                    subs(A, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L}); 
                   -subs(V3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L}); 
                   -subs(V2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L});
                    subs(T, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L});
                   -subs(M3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L})
                    subs(M2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L})];
end
K_EE=simplify(K_EE);
K_EE(:,5)=-K_EE(:,5);
K_EE(:,11)=-K_EE(:,11);
%% matriz de rigidez de referencia
      %%Axial X                        %% vertical Z                           tranvesal Y 
Kloc=[[  EA/L,                                    0,                                    0,     0,                                                   0,                                                   0, -EA/L,                                    0,                                    0,     0,                                                   0,                                                   0];
      [     0,  (12*Acy*EIz)/(L*(Acy*L^2 + 12*EIz)),                                    0,     0,                                                   0,                      (6*Acy*EIz)/(Acy*L^2 + 12*EIz),     0, -(12*Acy*EIz)/(L*(Acy*L^2 + 12*EIz)),                                    0,     0,                                                   0,                      (6*Acy*EIz)/(Acy*L^2 + 12*EIz)];
      [     0,                                    0,  (12*Acz*EIy)/(L*(Acz*L^2 + 12*EIy)),     0,                     -(6*Acz*EIy)/(Acz*L^2 + 12*EIy),                                                   0,     0,                                    0, -(12*Acz*EIy)/(L*(Acz*L^2 + 12*EIy)),     0,                     -(6*Acz*EIy)/(Acz*L^2 + 12*EIy),                                                   0];
      [     0,                                    0,                                    0,  GJ/L,                                                   0,                                                   0,     0,                                    0,                                    0, -GJ/L,                                                   0,                                                   0];
      [     0,                                    0,      -(6*Acz*EIy)/(Acz*L^2 + 12*EIy),     0,    (4*EIy*(Acz*L^2 + 3*EIy))/(L*(Acz*L^2 + 12*EIy)),                                                   0,     0,                                    0,       (6*Acz*EIy)/(Acz*L^2 + 12*EIy),     0, -(2*EIy*(- Acz*L^2 + 6*EIy))/(L*(Acz*L^2 + 12*EIy)),                                                   0];
      [     0,       (6*Acy*EIz)/(Acy*L^2 + 12*EIz),                                    0,     0,                                                   0,    (4*EIz*(Acy*L^2 + 3*EIz))/(L*(Acy*L^2 + 12*EIz)),     0,      -(6*Acy*EIz)/(Acy*L^2 + 12*EIz),                                    0,     0,                                                   0, -(2*EIz*(- Acy*L^2 + 6*EIz))/(L*(Acy*L^2 + 12*EIz))];
      [ -EA/L,                                    0,                                    0,     0,                                                   0,                                                   0,  EA/L,                                    0,                                    0,     0,                                                   0,                                                   0];
      [     0, -(12*Acy*EIz)/(L*(Acy*L^2 + 12*EIz)),                                    0,     0,                                                   0,                     -(6*Acy*EIz)/(Acy*L^2 + 12*EIz),     0,  (12*Acy*EIz)/(L*(Acy*L^2 + 12*EIz)),                                    0,     0,                                                   0,                     -(6*Acy*EIz)/(Acy*L^2 + 12*EIz)];
      [     0,                                    0, -(12*Acz*EIy)/(L*(Acz*L^2 + 12*EIy)),     0,                      (6*Acz*EIy)/(Acz*L^2 + 12*EIy),                                                   0,     0,                                    0,  (12*Acz*EIy)/(L*(Acz*L^2 + 12*EIy)),     0,                      (6*Acz*EIy)/(Acz*L^2 + 12*EIy),                                                   0];
      [     0,                                    0,                                    0, -GJ/L,                                                   0,                                                   0,     0,                                    0,                                    0,  GJ/L,                                                   0,                                                   0];
      [     0,                                    0,      -(6*Acz*EIy)/(Acz*L^2 + 12*EIy),     0, -(2*EIy*(- Acz*L^2 + 6*EIy))/(L*(Acz*L^2 + 12*EIy)),                                                   0,     0,                                    0,       (6*Acz*EIy)/(Acz*L^2 + 12*EIy),     0,    (4*EIy*(Acz*L^2 + 3*EIy))/(L*(Acz*L^2 + 12*EIy)),                                                   0];
      [     0,       (6*Acy*EIz)/(Acy*L^2 + 12*EIz),                                    0,     0,                                                   0, -(2*EIz*(- Acy*L^2 + 6*EIz))/(L*(Acy*L^2 + 12*EIz)),     0,      -(6*Acy*EIz)/(Acy*L^2 + 12*EIz),                                    0,     0,                                                   0,    (4*EIz*(Acy*L^2 + 3*EIz))/(L*(Acy*L^2 + 12*EIz))]];

  
   AcyEIz12=(12*Acy*EIz)/(L*(Acy*L^2 + 12*EIz));
   AcyEIz6=(6*Acy*EIz)/(Acy*L^2 + 12*EIz);
   AcyEIz4=(4*EIz*(Acy*L^2 + 3*EIz))/(L*(Acy*L^2 + 12*EIz));
   AcyEIz2=(2*EIz*(- Acy*L^2 + 6*EIz))/(L*(Acy*L^2 + 12*EIz));
   
   AczEIy12=(12*Acz*EIy)/(L*(Acz*L^2 + 12*EIy));
   AczEIy6=(6*Acz*EIy)/(Acz*L^2 + 12*EIy);
   AczEIy4=(4*EIy*(Acz*L^2 + 3*EIy))/(L*(Acz*L^2 + 12*EIy));
   AczEIy2=(2*EIy*(- Acz*L^2 + 6*EIy))/(L*(Acz*L^2 + 12*EIy));
   
Kloc1=[[  EA/L,         0,         0,     0,       0,        0, -EA/L,        0,         0,     0,        0,       0];
       [     0,  AcyEIz12,         0,     0,       0,  AcyEIz6,     0,-AcyEIz12,         0,     0,        0, AcyEIz6];
       [     0,         0,  AczEIy12,     0,-AczEIy6,        0,     0,        0, -AczEIy12,     0, -AczEIy6,       0];
       [     0,         0,         0,  GJ/L,       0,        0,     0,        0,         0, -GJ/L,        0,       0];
       [     0,         0,  -AczEIy6,     0, AczEIy4,        0,     0,        0,   AczEIy6,     0, -AczEIy2,       0];
       [     0,   AcyEIz6,         0,     0,       0,  AcyEIz4,     0, -AcyEIz6,         0,     0,        0,-AcyEIz2];
       [ -EA/L,         0,         0,     0,       0,        0,  EA/L,        0,         0,     0,        0,       0];
       [     0, -AcyEIz12,         0,     0,       0, -AcyEIz6,     0, AcyEIz12,         0,     0,        0,-AcyEIz6];
       [     0,         0, -AczEIy12,     0, AczEIy6,        0,     0,        0,  AczEIy12,     0,  AczEIy6,       0];
       [     0,         0,         0, -GJ/L,       0,        0,     0,        0,         0,  GJ/L,        0,       0];
       [     0,         0,  -AczEIy6,     0,-AczEIy2,        0,     0,        0,   AczEIy6,     0,  AczEIy4,       0];
       [     0,   AcyEIz6,         0,     0,       0, -AcyEIz2,     0, -AcyEIz6,         0,     0,        0, AcyEIz4]];
      %% matriz de rigides viga apoyada empotrada RE
K_RE = sym(zeros(12));

for i = 1:12
      [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12]=solve(subs(u ,x,0) == (i==1),...%%# Condiciones de frontera
                                                     subs(v2,x,0) == (i==2),...
                                                     subs(v3,x,0) == (i==3),...
                                                     subs(tt,x,0) == (i==4),...
                                                     subs(M3,x,0) == 0     ,...
                                                     subs(M2,x,0) == 0     ,...
                                                     subs(u ,x,L) == (i==7),...
                                                     subs(v2,x,L) == (i==8),...
                                                     subs(v3,x,L) == (i==9),...
                                                     subs(tt,x,L) == (i==10),...
                                                     subs(t3,x,L) == (i==11),...
                                                     subs(t2,x,L) == (i==12),...
                                                     [C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12]);                                         
        K_RE(:,i)=[-subs(A, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0}); 
                    subs(V3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0}); 
                    subs(V2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                   -subs(T, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                    subs(M3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                   -subs(M2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                    subs(A, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L}); 
                   -subs(V3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L}); 
                   -subs(V2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L});
                    subs(T, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L});
                   -subs(M3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L})
                    subs(M2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L})];
end
K_RE=simplify(K_RE);
K_RE(:,5)=-K_RE(:,5);
K_RE(:,11)=-K_RE(:,11);

%% matriz de rigides viga empotrada apoyada  ER
K_ER = sym(zeros(12));

for i = 1:12
      [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12]=solve(subs(u ,x,0) == (i==1),...%%# Condiciones de frontera
                                                     subs(v2,x,0) == (i==2),...
                                                     subs(v3,x,0) == (i==3),...
                                                     subs(tt,x,0) == (i==4),...
                                                     subs(t3,x,0) == (i==5),...
                                                     subs(t2,x,0) == (i==6),...
                                                     subs(u ,x,L) == (i==7),...
                                                     subs(v2,x,L) == (i==8),...
                                                     subs(v3,x,L) == (i==9),...
                                                     subs(tt,x,L) == (i==10),...
                                                     subs(M3,x,L) == 0      ,...
                                                     subs(M2,x,L) == 0      ,...
                                                     [C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12]);                                         
        K_ER(:,i)=[-subs(A, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0}); 
                    subs(V3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0}); 
                    subs(V2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                   -subs(T, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                    subs(M3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                   -subs(M2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                    subs(A, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L}); 
                   -subs(V3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L}); 
                   -subs(V2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L});
                    subs(T, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L});
                   -subs(M3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L})
                    subs(M2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L})];
end
K_ER=simplify(K_ER);
K_ER(:,5)=-K_ER(:,5);
K_ER(:,11)=-K_ER(:,11);
%% matriz de rigides viga apoyada apoyada  RR
K_RR = sym(zeros(12));

for i = 1:12
      [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12]=solve(subs(u ,x,0) == (i==1),...%%# Condiciones de frontera
                                                     subs(v2,x,0) == (i==2),...
                                                     subs(v3,x,0) == (i==3),...
                                                     subs(tt,x,0) == (i==4),...
                                                     subs(M3,x,0) == 0     ,...
                                                     subs(M2,x,0) == 0     ,...
                                                     subs(u ,x,L) == (i==7),...
                                                     subs(v2,x,L) == (i==8),...
                                                     subs(v3,x,L) == (i==9),...
                                                     subs(tt,x,L) == (i==10),...
                                                     subs(M3,x,L) == 0     ,...
                                                     subs(M2,x,L) == 0     ,...
                                                     [C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12]);                                         
        K_RR(:,i)=[-subs(A, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0}); 
                    subs(V3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0}); 
                    subs(V2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                   -subs(T, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                    subs(M3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                   -subs(M2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,0});
                    subs(A, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L}); 
                   -subs(V3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L}); 
                   -subs(V2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L});
                    subs(T, {C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L});
                   -subs(M3,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L})
                    subs(M2,{C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,x},{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,L})];
end
K_RR=simplify(K_RR);
K_RR(:,5)=-K_RR(:,5);
K_RR(:,11)=-K_RR(:,11);

