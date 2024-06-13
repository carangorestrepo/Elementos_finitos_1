
clc
clear

 %% Se define el elemento finito

EF = 10;

%# %%
%# Numeración local del EF triangular de 10 nodos:

%# Coordenadas de los nodos
%#                   xi  eta      # nodo   
xnod =          [[    0,   0 ];   %#  1
                  [   1,   0 ];   %#  2
                  [   0,   1 ];   %#  3
                  [ 1/3,   0 ];   %#  4
                  [ 2/3,   0 ];   %#  5
                  [ 2/3, 1/3 ];   %#  6
                  [ 1/3, 2/3 ];   %#  7
                  [   0, 2/3 ];   %#  8
                  [   0, 1/3 ];   %#  9
                  [ 1/3, 1/3 ]];  %#  10

plot(xnod(:,1),xnod(:,2),'*r')
              
nno = size(xnod, 1);

%# %% Se define la cuadratura de Gauss Legendre a utilizar
n = 4 ;%# Orden de la cuadratura (Lineal, cuadrática, cúbica, ...)
xw=TriGaussPoints(n);

x_gl = (xw(:,1));
e_gl = (xw(:,2));
w_gl =  xw(:,3);
n_gl = size(x_gl,1);


%% número de términos del polinomio interpolador
nterm = (n_gl); % 1   xi_gl   eta_gl    xi_gl²   xi_gl*eta_gl   eta_gl²   

 %%  Se define la matriz A1

A1 = (zeros(nterm));

for i =1 :n_gl
    xi = x_gl(i); 
    eta = e_gl(i);
    A1(i,:) =[1 xi eta xi^2 xi*eta eta^2];
end

 %% Se define la matriz A2
A2 = (zeros(nno, nterm));
for i = 1:nno
    xi_gl = xnod(i,1);
    eta_gl = xnod(i,2);
    A2(i,:) = [ 1, xi_gl, eta_gl, xi_gl^2, xi_gl*eta_gl, eta_gl^2 ];
end
 %% Se reporta la matriz
mensaje = 'EF triangular de 6 nodos';
fprintf('La matriz de interpolacion del %s es:\n', mensaje)
A = (A2/A1); %= A2*inv(A1)
disp(A) 
disp('o en forma numerica:')
format long
disp(double(A))
