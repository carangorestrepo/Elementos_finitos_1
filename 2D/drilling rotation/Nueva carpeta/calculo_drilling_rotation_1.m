clc
clear
E=30000;
nu=0.25;
t=0.2;

%% Func. de forma y sus derivadas del EF rectangular serendÃ­pito de 8 nodos
% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa deduccion_funciones_forma/FF_serendipitos_Q4_Q8.m
Nforma = @(xi,eta)[1 - xi - eta
                   xi
                   eta];
% derivadas de las funciones de forma con respecto a xi
dN_dxi = @(xi,eta) [ ...
                    -1
                    1
                    0];



%% Derivadas de N con respecto a eta
dN_deta = @(xi,eta) [ ... 
                    -1;
                    0
                    1];

xe=[0;1.5;0];
ye=[0;0;1.5];

% matriz constitutiva para TENSION PLANA
De = (E/(1-nu^2)) * [ 1      nu  0
                             nu  1      0
                             0      0      (1-nu)/2 ];

%% Parametros de la cuadratura de Gauss-Legendre
x_gl= [0.5;
       0.5;
       0];
      
e_gl = [0;
        0.5;
        0.5];  
   
w_gl =[1/6;
       1/6;
       1/6];
n_gl = size(x_gl,1);  %# Número de puntos de Gauss.
x1=xe(1);
x2=xe(2);
x3=xe(3);
y1=ye(1);
y2=ye(2);
y3=ye(3);

x21 = xe(2) - xe(1);         y21 = ye(2) - ye(1); 
x32 = xe(3) - xe(2);         y32 = ye(3) - ye(2);
x31 = xe(3) - xe(1);         y31 = ye(3) - ye(1);

xji = [ x21 x32 x31];   yji = [y21 y32 y31]; 

Lk = hypot(xji, yji);

L12=Lk(1);
L23=Lk(2);
L31=Lk(3);

a1=-(4*(x1*x2*y3^2 + x3^2*y1*y2 - x1*x3*y2*y3 - x2*x3*y1*y3))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
a2 =-(4*(x1*y2*y3 - x2*y3^2 - x1*y3^2 + x2*y1*y3 - 2*x3*y1*y2 + x3*y1*y3 + x3*y2*y3))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
a3 =-(4*(x1*x3*y2 - x3^2*y2 - 2*x1*x2*y3 - x3^2*y1 + x2*x3*y1 + x1*x3*y3 + x2*x3*y3))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
a4 =-(4*(y1*y2 - y1*y3 - y2*y3 + y3^2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
a5 =-(4*(x1*y3 - x2*y1 - x1*y2 + x3*y1 + x2*y3 + x3*y2 - 2*x3*y3))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
a6 =-(4*(x1 - x3)*(x2 - x3))/(x1^2*y2^2 - 2*x1^2*y2*y3 + x1^2*y3^2 - 2*x1*x2*y1*y2 + 2*x1*x2*y1*y3 + 2*x1*x2*y2*y3 - 2*x1*x2*y3^2 + 2*x1*x3*y1*y2 - 2*x1*x3*y1*y3 - 2*x1*x3*y2^2 + 2*x1*x3*y2*y3 + x2^2*y1^2 - 2*x2^2*y1*y3 + x2^2*y3^2 - 2*x2*x3*y1^2 + 2*x2*x3*y1*y2 + 2*x2*x3*y1*y3 - 2*x2*x3*y2*y3 + x3^2*y1^2 - 2*x3^2*y1*y2 + x3^2*y2^2);

b1=-(4*(x2*x3*y1^2 + x1^2*y2*y3 - x1*x2*y1*y3 - x1*x3*y1*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
b2=-(4*(x1*y1*y2 - x3*y1^2 - x2*y1^2 + x1*y1*y3 - 2*x1*y2*y3 + x2*y1*y3 + x3*y1*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2; 
b3=-(4*(x1*x2*y1 - x1^2*y3 - x1^2*y2 + x1*x3*y1 + x1*x2*y3 + x1*x3*y2 - 2*x2*x3*y1))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
b4=(4*(y1*y2 + y1*y3 - y2*y3 - y1^2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;  
b5 =-(4*(x1*y2 - 2*x1*y1 + x2*y1 + x1*y3 + x3*y1 - x2*y3 - x3*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
b6 =-(4*(x1 - x2)*(x1 - x3))/(x1^2*y2^2 - 2*x1^2*y2*y3 + x1^2*y3^2 - 2*x1*x2*y1*y2 + 2*x1*x2*y1*y3 + 2*x1*x2*y2*y3 - 2*x1*x2*y3^2 + 2*x1*x3*y1*y2 - 2*x1*x3*y1*y3 - 2*x1*x3*y2^2 + 2*x1*x3*y2*y3 + x2^2*y1^2 - 2*x2^2*y1*y3 + x2^2*y3^2 - 2*x2*x3*y1^2 + 2*x2*x3*y1*y2 + 2*x2*x3*y1*y3 - 2*x2*x3*y2*y3 + x3^2*y1^2 - 2*x3^2*y1*y2 + x3^2*y2^2);

d1=-(4*(x1*x3*y2^2 + x2^2*y1*y3 - x1*x2*y2*y3 - x2*x3*y1*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
d2=-(4*(x2*y1*y2 - x3*y2^2 - x1*y2^2 + x1*y2*y3 - 2*x2*y1*y3 + x3*y1*y2 + x2*y2*y3))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
d3 =-(4*(x1*x2*y2 - x2^2*y3 - x2^2*y1 + x1*x2*y3 - 2*x1*x3*y2 + x2*x3*y1 + x2*x3*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2;
d4=(4*(y1*y2 - y1*y3 + y2*y3 - y2^2))/(x1^2*y2^2 - 2*x1^2*y2*y3 + x1^2*y3^2 - 2*x1*x2*y1*y2 + 2*x1*x2*y1*y3 + 2*x1*x2*y2*y3 - 2*x1*x2*y3^2 + 2*x1*x3*y1*y2 - 2*x1*x3*y1*y3 - 2*x1*x3*y2^2 + 2*x1*x3*y2*y3 + x2^2*y1^2 - 2*x2^2*y1*y3 + x2^2*y3^2 - 2*x2*x3*y1^2 + 2*x2*x3*y1*y2 + 2*x2*x3*y1*y3 - 2*x2*x3*y2*y3 + x3^2*y1^2 - 2*x3^2*y1*y2 + x3^2*y2^2);
d5 =-(4*(x1*y2 + x2*y1 - x1*y3 - 2*x2*y2 - x3*y1 + x2*y3 + x3*y2))/(x1^2*y2^2 - 2*x1^2*y2*y3 + x1^2*y3^2 - 2*x1*x2*y1*y2 + 2*x1*x2*y1*y3 + 2*x1*x2*y2*y3 - 2*x1*x2*y3^2 + 2*x1*x3*y1*y2 - 2*x1*x3*y1*y3 - 2*x1*x3*y2^2 + 2*x1*x3*y2*y3 + x2^2*y1^2 - 2*x2^2*y1*y3 + x2^2*y3^2 - 2*x2*x3*y1^2 + 2*x2*x3*y1*y2 + 2*x2*x3*y1*y3 - 2*x2*x3*y2*y3 + x3^2*y1^2 - 2*x3^2*y1*y2 + x3^2*y2^2);
d6 =(4*(x1 - x2)*(x2 - x3))/(x1^2*y2^2 - 2*x1^2*y2*y3 + x1^2*y3^2 - 2*x1*x2*y1*y2 + 2*x1*x2*y1*y3 + 2*x1*x2*y2*y3 - 2*x1*x2*y3^2 + 2*x1*x3*y1*y2 - 2*x1*x3*y1*y3 - 2*x1*x3*y2^2 + 2*x1*x3*y2*y3 + x2^2*y1^2 - 2*x2^2*y1*y3 + x2^2*y3^2 - 2*x2*x3*y1^2 + 2*x2*x3*y1*y2 + 2*x2*x3*y1*y3 - 2*x2*x3*y2*y3 + x3^2*y1^2 - 2*x3^2*y1*y2 + x3^2*y2^2);

x21 = xe(2) - xe(1);         y21 = ye(2) - ye(1); 
x32 = xe(3) - xe(2);         y32 = ye(3) - ye(2);
x31 = xe(3) - xe(1);         y31 = ye(3) - ye(1);

c12 =  (y2 - y1)/L12;    c23 =  (y3 - y2)/L23;    c31 =  (y1 - y3)/L31;
s12 = -(x2 - x1)/L12;    s23 = -(x3 - x2)/L23;    s31 = -(x1 - x3)/L31;

H= [[ (x2*y3 - x3*y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                               0,      0, -(x1*y3 - x3*y1)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                                0,      0, (x1*y2 - x2*y1)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                               0,      0];
    [       (y2 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                               0,      0,       -(y1 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                                0,      0,       (y1 - y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                               0,      0];
    [      -(x2 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                               0,      0,        (x1 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                                0,      0,      -(x1 - x2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),                                                               0,      0];
    [                                                               0, (x2*y3 - x3*y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),      0,                                                                0, -(x1*y3 - x3*y1)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),      0,                                                               0, (x1*y2 - x2*y1)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),      0];
    [                                                               0,       (y2 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),      0,                                                                0,       -(y1 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),      0,                                                               0,       (y1 - y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),      0];
    [                                                               0,      -(x2 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),      0,                                                                0,        (x1 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),      0,                                                               0,      -(x1 - x2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),      0];
    [                                                               0,                                                               0, -L12/6,                                                                0,                                                                0,  L12/6,                                                               0,                                                               0,      0];
    [                                                               0,                                                               0,      0,                                                                0,                                                                0, -L23/6,                                                               0,                                                               0,  L23/6];
    [                                                               0,                                                               0,  L31/6,                                                                0,                                                                0,      0,                                                               0,                                                               0, -L31/6]];

Ke = zeros(3*3);
fe = zeros(3*3,1);
Me=  zeros(3*2,1);
det_Je = zeros(n_gl,1); % matriz para almacenar los jacobianos
nef=1
B   = cell(nef,n_gl,n_gl); % contenedor para las matrices de deformacion
for e = 1:nef
    for p = 1:n_gl

        xi_gl  = x_gl(p);
        eta_gl = e_gl(p);
        % Se evaluan las funciones de forma y sus derivadas 
        % en los puntos de integracion de Gauss-Legendre
        NNforma  = Nforma (xi_gl, eta_gl);
        ddN_dxi  = dN_dxi (xi_gl, eta_gl);
        ddN_deta = dN_deta(xi_gl, eta_gl);

        dx_dxi  = sum(ddN_dxi  .* xe);   dy_dxi  = sum(ddN_dxi  .* ye);
        dx_deta = sum(ddN_deta .* xe);   dy_deta = sum(ddN_deta .* ye);
        % Se ensambla la matriz Jacbiana del elemento
        Je = [ dx_dxi   dy_dxi
                dx_deta  dy_deta ];      
        % Se calcula el determinante del Jacobiano
        det_Je(p) = det(Je);

        x  = sum(NNforma  .* xe);   y  = sum(NNforma  .* ye);

        B{e,p}=[[ 0, 1, 0, 0, 0, 0,                            c12*(a2 + 2*a4*x + a5*y),                                                                                                                                                            c23*(b5*y - (4*(x1*y1*y2 - x3*y1^2 - x2*y1^2 + x1*y1*y3 - 2*x1*y2*y3 + x2*y1*y3 + x3*y1*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2 + (8*x*(y1*y2 + y1*y3 - y2*y3 - y1^2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2),                            c31*(d5*y - (4*(x2*y1*y2 - x3*y2^2 - x1*y2^2 + x1*y2*y3 - 2*x2*y1*y3 + x3*y1*y2 + x2*y2*y3))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2 + (8*x*(y1*y2 - y1*y3 + y2*y3 - y2^2))/(x1^2*y2^2 - 2*x1^2*y2*y3 + x1^2*y3^2 - 2*x1*x2*y1*y2 + 2*x1*x2*y1*y3 + 2*x1*x2*y2*y3 - 2*x1*x2*y3^2 + 2*x1*x3*y1*y2 - 2*x1*x3*y1*y3 - 2*x1*x3*y2^2 + 2*x1*x3*y2*y3 + x2^2*y1^2 - 2*x2^2*y1*y3 + x2^2*y3^2 - 2*x2*x3*y1^2 + 2*x2*x3*y1*y2 + 2*x2*x3*y1*y3 - 2*x2*x3*y2*y3 + x3^2*y1^2 - 2*x3^2*y1*y2 + x3^2*y2^2))];
                [ 0, 0, 0, 0, 0, 1,                            s12*(a3 + a5*x + 2*a6*y),                                                                                                                                                                                                                                            s23*(b5*x + 2*b6*y - (4*(x1*x2*y1 - x1^2*y3 - x1^2*y2 + x1*x3*y1 + x1*x2*y3 + x1*x3*y2 - 2*x2*x3*y1))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            s31*(d3 + d5*x + 2*d6*y)];
                [ 0, 0, 1, 0, 1, 0, s12*(a2 + 2*a4*x + a5*y) + c12*(a3 + a5*x + 2*a6*y), s23*(b5*y - (4*(x1*y1*y2 - x3*y1^2 - x2*y1^2 + x1*y1*y3 - 2*x1*y2*y3 + x2*y1*y3 + x3*y1*y2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2 + (8*x*(y1*y2 + y1*y3 - y2*y3 - y1^2))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2) + c23*(b5*x + 2*b6*y - (4*(x1*x2*y1 - x1^2*y3 - x1^2*y2 + x1*x3*y1 + x1*x2*y3 + x1*x3*y2 - 2*x2*x3*y1))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2), s31*(d5*y - (4*(x2*y1*y2 - x3*y2^2 - x1*y2^2 + x1*y2*y3 - 2*x2*y1*y3 + x3*y1*y2 + x2*y2*y3))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)^2 + (8*x*(y1*y2 - y1*y3 + y2*y3 - y2^2))/(x1^2*y2^2 - 2*x1^2*y2*y3 + x1^2*y3^2 - 2*x1*x2*y1*y2 + 2*x1*x2*y1*y3 + 2*x1*x2*y2*y3 - 2*x1*x2*y3^2 + 2*x1*x3*y1*y2 - 2*x1*x3*y1*y3 - 2*x1*x3*y2^2 + 2*x1*x3*y2*y3 + x2^2*y1^2 - 2*x2^2*y1*y3 + x2^2*y3^2 - 2*x2*x3*y1^2 + 2*x2*x3*y1*y2 + 2*x2*x3*y1*y3 - 2*x2*x3*y2*y3 + x3^2*y1^2 - 2*x3^2*y1*y2 + x3^2*y2^2)) + c31*(d3 + d5*x + 2*d6*y)]];
            
        Ke = Ke + B{e,p}'*De*B{e,p}*det_Je(p)*w_gl(p);
    end
   % y se ensambla la matriz de rigidez del elemento y el vector de fuerzas
   % nodales del elemento en sus correspondientes GDL 
   %idx{e}           = reshape(gdl(LaG(e,:),:)', 1, 3*2);
   %idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2),:) ...
   %          gdl(LaG(e,3),:)  gdl(LaG(e,4),:) ...
   %          gdl(LaG(e,5),:)  gdl(LaG(e,6),:) ...
   %          gdl(LaG(e,7),:)  gdl(LaG(e,8),:) gdl(LaG(e,9),:) ];
   
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + H'*Ke*H;
end




