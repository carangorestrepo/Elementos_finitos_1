
clc
clear

 %% Se define el elemento finito


%# %%
%# Numeración local del EF triangular de 6 nodos:

%# Coordenadas de los nodos
%#                   xi  eta      # nodo   
xnod =          [ [   0,   0 ];   %#  1
                  [   1,   0 ];   %#  2
                  [   0,   1 ];   %#  3
                  [ 1/4,   0 ];   %#  4
                  [ 1/2,   0 ];   %#  5
                  [ 3/4,   0 ];   %#  6
                  [ 3/4, 1/4 ];   %#  7
                  [ 1/2, 1/2 ];   %#  8
                  [ 1/4, 3/4 ];   %#  9
                  [   0, 3/4 ];   %#  10
                  [   0, 1/2 ];   %#  11
                  [   0, 1/4 ];   %#  12
                  [ 1/4, 1/4 ];   %#  13
                  [ 1/2, 1/4 ];   %#  14
                  [ 1/4, 1/2 ]];  %#  15

             
plot(xnod(:,1),xnod(:,2),'*r')
              
nno = size(xnod, 1);

%# %% Se define la cuadratura de Gauss Legendre a utilizar
n = 6 ;%# Orden de la cuadratura (Lineal, cuadrática, cúbica, ...)
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
    %A1(i,:) =[1 xi eta xi^2 xi*eta eta^2];
    A1(i,:) =[1 xi eta xi^2 xi*eta eta^2 xi^3 xi^2*eta xi*eta^2 eta^3 eta^3*xi eta*xi^3];
end

 %% Se define la matriz A2
A2 = (zeros(nno, nterm));
for i = 1:nno
    xi_gl = xnod(i,1);
    eta_gl = xnod(i,2);
    A2(i,:) = [ 1, xi_gl, eta_gl,...
                 xi_gl^2, xi_gl*eta_gl, eta_gl^2,...
                 xi_gl^3 xi_gl^2*eta_gl xi_gl*eta_gl^2 eta_gl^3,...
                eta_gl^3*xi_gl eta_gl*xi_gl^3];
end
 %% Se reporta la matriz
mensaje = 'EF triangular de 6 nodos';
fprintf('La matriz de interpolacion del %s es:\n', mensaje)
A = (A2/A1); %= A2*inv(A1)
disp(A) 
disp('o en forma numerica:')
format long
disp(double(A))
