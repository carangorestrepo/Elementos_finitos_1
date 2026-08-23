syms C1 C2 C3 C4 C5 C6 alfa beta k Ac q1 q2 x L EI AE b1 b2  real

bx = b1 + (b2-b1)*x/L;
%% Carga trapezoidal vertical 
qx = q1 + (q2-q1)*x/L;

%% Solucion particular
vp = qx/k;

%% Funciones base
f1 = cosh(alfa*x)*cos(beta*x);
f2 = cosh(alfa*x)*sin(beta*x);
f3 = sinh(alfa*x)*cos(beta*x);
f4 = sinh(alfa*x)*sin(beta*x);

%% Constantes de integracion
Un = [C1; C2; C3; C4; C5; C6];

%% Solucion homogenea
vh = f1*C1+f2*C2+f3*C3+f4*C4;

%% Desplazamiento total
v = vh + vp;

%% Momento
M = EI*(diff(v,x,2) - k*v/Ac + qx/Ac);

%% Cortante
V = diff(M,x);

%% Giro
t = diff(v,x) + V/Ac;

A=int(bx,x)+C5;
u=int(A/AE,x)+C6;

%% Matriz asociada a las constantes de integracion
An = equationsToMatrix([u; v; t; A; V; M],Un);
An = simplify(An);


An =[[                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                      0, x/AE, 1];
     [                                                                                                                                                                                                                                                                               cosh(alfa*x)*cos(beta*x),                                                                                                                                                                                                                                                                               cosh(alfa*x)*sin(beta*x),                                                                                                                                                                                                                                                                               cos(beta*x)*sinh(alfa*x),                                                                                                                                                                                                                                                                               sinh(alfa*x)*sin(beta*x),    0, 0];
     [ alfa*cos(beta*x)*sinh(alfa*x) - beta*cosh(alfa*x)*sin(beta*x) - (EI*((k*(alfa*cos(beta*x)*sinh(alfa*x) - beta*cosh(alfa*x)*sin(beta*x)))/Ac - beta^3*cosh(alfa*x)*sin(beta*x) - alfa^3*cos(beta*x)*sinh(alfa*x) + 3*alfa*beta^2*cos(beta*x)*sinh(alfa*x) + 3*alfa^2*beta*cosh(alfa*x)*sin(beta*x)))/Ac, alfa*sinh(alfa*x)*sin(beta*x) - (EI*(beta^3*cosh(alfa*x)*cos(beta*x) - alfa^3*sinh(alfa*x)*sin(beta*x) + (k*(alfa*sinh(alfa*x)*sin(beta*x) + beta*cosh(alfa*x)*cos(beta*x)))/Ac - 3*alfa^2*beta*cosh(alfa*x)*cos(beta*x) + 3*alfa*beta^2*sinh(alfa*x)*sin(beta*x)))/Ac + beta*cosh(alfa*x)*cos(beta*x), (EI*(alfa^3*cosh(alfa*x)*cos(beta*x) + beta^3*sinh(alfa*x)*sin(beta*x) + (k*(beta*sinh(alfa*x)*sin(beta*x) - alfa*cosh(alfa*x)*cos(beta*x)))/Ac - 3*alfa*beta^2*cosh(alfa*x)*cos(beta*x) - 3*alfa^2*beta*sinh(alfa*x)*sin(beta*x)))/Ac - beta*sinh(alfa*x)*sin(beta*x) + alfa*cosh(alfa*x)*cos(beta*x), alfa*cosh(alfa*x)*sin(beta*x) + beta*cos(beta*x)*sinh(alfa*x) - (EI*(beta^3*cos(beta*x)*sinh(alfa*x) - alfa^3*cosh(alfa*x)*sin(beta*x) + (k*(alfa*cosh(alfa*x)*sin(beta*x) + beta*cos(beta*x)*sinh(alfa*x)))/Ac + 3*alfa*beta^2*cosh(alfa*x)*sin(beta*x) - 3*alfa^2*beta*cos(beta*x)*sinh(alfa*x)))/Ac,    0, 0];
     [                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                      0,    1, 0];
     [                                                                     -EI*((k*(alfa*cos(beta*x)*sinh(alfa*x) - beta*cosh(alfa*x)*sin(beta*x)))/Ac - beta^3*cosh(alfa*x)*sin(beta*x) - alfa^3*cos(beta*x)*sinh(alfa*x) + 3*alfa*beta^2*cos(beta*x)*sinh(alfa*x) + 3*alfa^2*beta*cosh(alfa*x)*sin(beta*x)),                                                                     -EI*(beta^3*cosh(alfa*x)*cos(beta*x) - alfa^3*sinh(alfa*x)*sin(beta*x) + (k*(alfa*sinh(alfa*x)*sin(beta*x) + beta*cosh(alfa*x)*cos(beta*x)))/Ac - 3*alfa^2*beta*cosh(alfa*x)*cos(beta*x) + 3*alfa*beta^2*sinh(alfa*x)*sin(beta*x)),                                                                      EI*(alfa^3*cosh(alfa*x)*cos(beta*x) + beta^3*sinh(alfa*x)*sin(beta*x) + (k*(beta*sinh(alfa*x)*sin(beta*x) - alfa*cosh(alfa*x)*cos(beta*x)))/Ac - 3*alfa*beta^2*cosh(alfa*x)*cos(beta*x) - 3*alfa^2*beta*sinh(alfa*x)*sin(beta*x)),                                                                     -EI*(beta^3*cos(beta*x)*sinh(alfa*x) - alfa^3*cosh(alfa*x)*sin(beta*x) + (k*(alfa*cosh(alfa*x)*sin(beta*x) + beta*cos(beta*x)*sinh(alfa*x)))/Ac + 3*alfa*beta^2*cosh(alfa*x)*sin(beta*x) - 3*alfa^2*beta*cos(beta*x)*sinh(alfa*x)),    0, 0];
     [                                                                                                                                                       -EI*(beta^2*cosh(alfa*x)*cos(beta*x) - alfa^2*cosh(alfa*x)*cos(beta*x) + 2*alfa*beta*sinh(alfa*x)*sin(beta*x) + (k*cosh(alfa*x)*cos(beta*x))/Ac),                                                                                                                                                        EI*(alfa^2*cosh(alfa*x)*sin(beta*x) - beta^2*cosh(alfa*x)*sin(beta*x) + 2*alfa*beta*cos(beta*x)*sinh(alfa*x) - (k*cosh(alfa*x)*sin(beta*x))/Ac),                                                                                                                                                       -EI*(beta^2*cos(beta*x)*sinh(alfa*x) - alfa^2*cos(beta*x)*sinh(alfa*x) + 2*alfa*beta*cosh(alfa*x)*sin(beta*x) + (k*cos(beta*x)*sinh(alfa*x))/Ac),                                                                                                                                                        EI*(alfa^2*sinh(alfa*x)*sin(beta*x) - beta^2*sinh(alfa*x)*sin(beta*x) + 2*alfa*beta*cosh(alfa*x)*cos(beta*x) - (k*sinh(alfa*x)*sin(beta*x))/Ac),    0, 0]];
 